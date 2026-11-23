#!/usr/bin/env python3
"""Fast optimizer for Fusion post-processor parameters on MITDB/SVDB + synthetic noise.

MinEl and TiTen remain 100% frozen; only the 6 post-processor parameters in
tenhz.fusion.FUSION (agree_s, w_single, lam, start_pen, win_s, thr) are tuned.
Records 118 and 119 are excluded so NSTDB base patients stay completely unseen.

Uses a Numba-compiled two-pointer closest-first matcher plus Latin-Hypercube
global search + multi-start bounded Powell refinement (~250 evals, ~25-35s total).

Usage:
    python scripts/optimize_fusion.py --dbs mitdb svdb --noise-recs 16
    python scripts/optimize_fusion.py --dbs mitdb svdb --noise-recs 16 --cv
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from numba import njit
from scipy.optimize import minimize
from scipy.signal import butter, filtfilt

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))

import tenhz  # noqa: E402
from tenhz._signal import filtfilt_matlab  # noqa: E402
from tenhz.detectors import minel_detect, minmax_stream, titen_detect  # noqa: E402
from tenhz.fusion import FUSION, adaptive, pair_flags, viterbi  # noqa: E402
from tenhz.params import PARAMS  # noqa: E402


def mround(x):
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


@njit(cache=True)
def _k_match_sorted(det, gt, tol):
    """Numba two-pointer closest-first one-to-one matcher (exact match to benchmark_three.match)."""
    nd = len(det)
    ng = len(gt)
    if nd == 0 or ng == 0:
        return 0, nd, ng
    tot = 0
    j_lo = 0
    for i in range(nd):
        d_val = det[i]
        while j_lo < ng and gt[j_lo] < d_val - tol:
            j_lo += 1
        j = j_lo
        while j < ng and gt[j] <= d_val + tol:
            tot += 1
            j += 1
    if tot == 0:
        return 0, nd, ng

    diffs = np.empty(tot, np.int64)
    ii = np.empty(tot, np.int32)
    jj = np.empty(tot, np.int32)
    idx = 0
    j_lo = 0
    for i in range(nd):
        d_val = det[i]
        while j_lo < ng and gt[j_lo] < d_val - tol:
            j_lo += 1
        j = j_lo
        while j < ng and gt[j] <= d_val + tol:
            diffs[idx] = abs(d_val - gt[j])
            ii[idx] = i
            jj[idx] = j
            idx += 1
            j += 1

    order = np.argsort(diffs, kind="mergesort")
    used_d = np.zeros(nd, np.uint8)
    used_g = np.zeros(ng, np.uint8)
    tp = 0
    for z in range(tot):
        k = order[z]
        i = ii[k]
        j = jj[k]
        if used_d[i] == 0 and used_g[j] == 0:
            used_d[i] = 1
            used_g[j] = 1
            tp += 1
    return tp, nd - tp, ng - tp


def load_records(data_dir, dbs):
    """Exact record loader and exclusions from benchmark_three.py."""
    out = []
    for db in dbs:
        d = os.path.join(data_dir, f"{db}-mat")
        mf = os.path.join(d, "_manifest.json")
        if not os.path.isdir(d) or not os.path.exists(mf):
            sys.exit(f"Missing {d} or {mf}")
        man = json.load(open(mf))
        names = sorted(
            f[:-4] for f in os.listdir(d)
            if f.endswith(".txt") and not f.endswith(("-ann.txt", "-fs.txt"))
        )
        if db == "mitdb":
            names = [n for n in names if n not in ("118", "119")]
        for n in names:
            if man.get(n, {}).get("status") != "OK":
                continue
            fs = float(man[n]["fs"])
            base = os.path.join(d, n)
            ecg = pd.read_csv(base + ".txt", header=None, sep=r"\s+", usecols=[0], engine="c").to_numpy(dtype=float, copy=True).ravel()
            ecg[~np.isfinite(ecg)] = 0.0
            ecg = ecg.astype(np.float32).astype(float)
            ann = pd.read_csv(base + "-ann.txt", header=None, sep=r"\s+", usecols=[0], engine="c").to_numpy(dtype=float, copy=True).ravel()
            ann = np.round(ann).astype(np.int64) - 1
            ann = np.sort(ann[(ann >= 0) & (ann < len(ecg))])
            out.append((f"{db}/{n}", ecg, ann, fs))
    return out


def synthesize_physiological_noise(n, fs, rng, mode="mixed"):
    """Non-stationary physiological ECG noise (electrode motion + muscle + baseline wander)."""
    nyq = 0.5 * fs
    t = np.arange(n) / fs

    bw = (0.6 * np.sin(2 * np.pi * rng.uniform(0.12, 0.35) * t + rng.uniform(0, 6.28)) +
          0.4 * np.sin(2 * np.pi * rng.uniform(0.40, 0.75) * t + rng.uniform(0, 6.28)))

    env = np.zeros(n)
    idx = 0
    while idx < n:
        dur = int(rng.uniform(4.0, 18.0) * fs)
        if rng.random() < 0.45:
            burst = np.hanning(min(dur, n - idx))
            env[idx:idx + len(burst)] = burst * rng.uniform(0.5, 1.5)
        idx += dur

    b_em, a_em = butter(2, [0.5 / nyq, min(15.0 / nyq, 0.95)], btype="bandpass")
    em = filtfilt(b_em, a_em, rng.standard_normal(n)) * (0.25 + env)

    hi_ma = min(120.0, 0.9 * nyq)
    b_ma, a_ma = butter(2, [20.0 / nyq, hi_ma / nyq], btype="bandpass")
    ma = filtfilt(b_ma, a_ma, rng.standard_normal(n)) * (0.20 + 0.8 * env)

    if mode == "em":
        return em + 0.3 * bw
    if mode == "ma":
        return ma + 0.2 * bw
    return 0.55 * em + 0.35 * ma + 0.25 * bw


def mix_at_snr(ecg, r_peaks, fs, snr_db, rng, mode="mixed"):
    """Scale noise relative to median R-wave peak-to-peak amplitude."""
    noise = synthesize_physiological_noise(len(ecg), fs, rng, mode=mode)
    hw = int(0.05 * fs)
    qrs_amps = [np.ptp(ecg[max(0, r - hw):min(len(ecg), r + hw + 1)]) for r in r_peaks]
    qrs_rms = (np.median(qrs_amps) if qrs_amps else np.std(ecg)) / (2.0 * np.sqrt(2.0))
    noise_rms = np.sqrt(np.mean(noise ** 2)) + 1e-12
    target_noise_rms = qrs_rms / (10.0 ** (snr_db / 20.0))
    return (ecg + noise * (target_noise_rms / noise_rms)).astype(np.float32).astype(float)


def precompute_front_end(records, band=(7, 32), label=""):
    """Run frozen MinEl and TiTen once per record/condition and cache (m, t)."""
    cache = []
    order = int(PARAMS["minel"]["order"])
    p_minel = PARAMS["minel"]["params"]
    p_titen = PARAMS["titen"]["params"]
    n_tot = len(records)
    for idx_r, (rec_id, ecg, r_true, fs) in enumerate(records, 1):
        b, a = butter(order, [band[0] / (fs / 2), band[1] / (fs / 2)], btype="bandpass")
        xf = filtfilt_matlab(b, a, ecg)
        val, idx = minmax_stream(xf, fs, 10)
        tsec = (idx + 1) / fs
        m = np.unique(idx[minel_detect(val, 10, p_minel, tsec)[0]])
        t = np.unique(idx[titen_detect(val, 10, p_titen, tsec)[0]])
        cache.append((rec_id, m, t, r_true, fs, int(round(21.0 / 1000.0 * fs)), mround(75.0 / 1000.0 * fs)))
        if label and (idx_r % 5 == 0 or idx_r == n_tot):
            print(f"\r  Precomputing {label:8s}: {idx_r:3d}/{n_tot}", end="", flush=True)
    if label:
        print("", flush=True)
    return cache


REGIME_WEIGHTS = {
    "clean": 0.35,
    "snr_12db": 0.25,
    "snr_6db": 0.20,
    "snr_0db": 0.12,
    "snr_m3db": 0.08,
}


def eval_cache(theta, cache_by_regime, tol_ms=75.0, shift_ms=21.0):
    """Evaluate pooled F1 across regimes."""
    agree_s, w_single, lam, start_pen, win_s, thr = theta
    use_cached_offsets = (tol_ms == 75.0 and shift_ms == 21.0)
    scores = {}
    w_sum = 0.0
    comp_sum = 0.0
    for regime, items in cache_by_regime.items():
        if not items:
            continue
        tp_sum = fp_sum = fn_sum = 0
        for _, m, t, r_true, fs, sh_samp, tol_samp in items:
            pm, pt = pair_flags(m, t, int(round(agree_s * fs)))
            V = viterbi(m, t, pm, pt, fs, w_single, lam, start_pen)
            out = adaptive(m, t, V, pm, pt, fs, win_s, thr)[0]
            if use_cached_offsets:
                out = np.maximum(out - sh_samp, 0)
                tp, fp, fn = _k_match_sorted(out, r_true, tol_samp)
            else:
                if shift_ms != 0.0:
                    out = np.maximum(out - int(round(shift_ms / 1000.0 * fs)), 0)
                tp, fp, fn = _k_match_sorted(out, r_true, mround(tol_ms / 1000.0 * fs))
            tp_sum += tp
            fp_sum += fp
            fn_sum += fn
        f1 = 100.0 * 2.0 * tp_sum / max(2 * tp_sum + fp_sum + fn_sum, 1)
        scores[regime] = f1
        w = REGIME_WEIGHTS[regime]
        comp_sum += w * f1
        w_sum += w

    composite = comp_sum / max(w_sum, 1e-12)
    return composite, scores


def build_regimes(clean_records, noise_recs=16, seed=42):
    rng = np.random.default_rng(seed)
    n_tot = len(clean_records)
    if noise_recs is None or noise_recs >= n_tot:
        noisy_idx = np.arange(n_tot)
    else:
        noisy_idx = np.linspace(0, n_tot - 1, int(noise_recs), dtype=int)

    datasets = {"clean": [], "snr_12db": [], "snr_6db": [], "snr_0db": [], "snr_m3db": []}
    for rec_id, ecg, r_true, fs in clean_records:
        datasets["clean"].append((rec_id, ecg, r_true, fs))

    for i in noisy_idx:
        rec_id, ecg, r_true, fs = clean_records[i]
        datasets["snr_12db"].append((rec_id, mix_at_snr(ecg, r_true, fs, 12.0, rng, "mixed"), r_true, fs))
        datasets["snr_6db"].append((rec_id, mix_at_snr(ecg, r_true, fs, 6.0, rng, "mixed"), r_true, fs))
        datasets["snr_0db"].append((rec_id, mix_at_snr(ecg, r_true, fs, 0.0, rng, "em"), r_true, fs))
        datasets["snr_m3db"].append((rec_id, mix_at_snr(ecg, r_true, fs, -3.0, rng, "em"), r_true, fs))

    print(f"Injecting 4 noise SNRs into {len(noisy_idx)}/{n_tot} records (keeping all {n_tot} clean).")
    return {k: precompute_front_end(v, label=k) for k, v in datasets.items()}


def latin_hypercube(n_samples, bounds, rng):
    dim = len(bounds)
    lo = np.array([b[0] for b in bounds], dtype=float)
    hi = np.array([b[1] for b in bounds], dtype=float)
    u = np.empty((n_samples, dim), dtype=float)
    for j in range(dim):
        perm = rng.permutation(n_samples)
        u[:, j] = (perm + rng.random(n_samples)) / n_samples
    return lo + u * (hi - lo)


def fit_fast(cache_by_regime, tol_ms=75.0, shift_ms=21.0, n_lhs=72, seed=42):
    """Latin-Hypercube global search + multi-start bounded Powell refinement (~250 evals)."""
    bounds = [
        (0.05, 0.14),   # agree_s
        (0.15, 0.75),   # w_single
        (0.25, 2.00),   # lam
        (0.30, 2.50),   # start_pen
        (5.00, 20.0),   # win_s
        (0.01, 0.25),   # thr
    ]
    lo = np.array([b[0] for b in bounds], dtype=float)
    hi = np.array([b[1] for b in bounds], dtype=float)
    default_x = np.array([FUSION["agree_s"], FUSION["w_single"], FUSION["lam"], FUSION["start_pen"], FUSION["win_s"], FUSION["thr"]], dtype=float)

    rng = np.random.default_rng(seed)
    candidates = np.vstack([default_x[None, :], latin_hypercube(n_lhs, bounds, rng)])
    n_tot = len(candidates)

    t_start = time.time()
    scored = []
    best_sc = -np.inf
    best_x = default_x.copy()

    for idx_c, x in enumerate(candidates, 1):
        sc = eval_cache(x, cache_by_regime, tol_ms=tol_ms, shift_ms=shift_ms)[0]
        scored.append((sc, x))
        if sc > best_sc:
            best_sc = sc
            best_x = x.copy()
        elapsed = time.time() - t_start
        eta = (elapsed / idx_c) * (n_tot - idx_c)
        width = 20
        filled = int(width * idx_c / n_tot)
        bar = "█" * filled + "░" * (width - filled)
        print(
            f"\r  Stage 1 LHS [{bar}] {idx_c:2d}/{n_tot} | best={best_sc:.3f} | "
            f"{elapsed:4.1f}s (ETA {eta:3.0f}s)",
            end="",
            flush=True,
        )
    print("", flush=True)

    # Pick top 2 diverse starting points from LHS + default
    scored.sort(key=lambda t: -t[0])
    starts = [scored[0][1]]
    span = hi - lo
    for _, x in scored[1:]:
        if np.linalg.norm((x - starts[0]) / span) > 0.25:
            starts.append(x)
            if len(starts) == 2:
                break

    eval_count = [0]
    t_pow = time.time()

    def obj(th):
        nonlocal best_sc, best_x
        th_c = np.clip(th, lo, hi)
        sc = eval_cache(th_c, cache_by_regime, tol_ms=tol_ms, shift_ms=shift_ms)[0]
        eval_count[0] += 1
        if sc > best_sc:
            best_sc = sc
            best_x = th_c.copy()
        if eval_count[0] % 5 == 0:
            print(
                f"\r  Stage 2 Powell refine | eval {eval_count[0]:3d} | best={best_sc:.3f} | "
                f"{time.time() - t_pow:4.1f}s",
                end="",
                flush=True,
            )
        return -sc

    for s_idx, x0 in enumerate(starts, 1):
        minimize(
            obj,
            x0,
            method="Powell",
            bounds=bounds,
            options={"maxfev": 140, "xtol": 2e-3, "ftol": 1e-4},
        )

    print(f"\n  Done in {time.time() - t_start:.1f}s ({n_tot + eval_count[0]} total evals).", flush=True)
    return best_x


def main():
    default_data = "../tenhz-qrs/data" if os.path.isdir("../tenhz-qrs/data") else "data"
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("TENHZ_DATA") or default_data)
    ap.add_argument("--dbs", nargs="+", default=["mitdb"], help="training DBs (default: mitdb; or mitdb svdb)")
    ap.add_argument("--noise-recs", type=int, default=16, help="number of records to inject noise into (default: 16; 0=all)")
    ap.add_argument("--tol-ms", type=float, default=75.0, help="target matching tolerance in ms (default: 75.0)")
    ap.add_argument("--shift-ms", type=float, default=21.0, help="fixed causal latency offset in ms (default: 21.0)")
    ap.add_argument("--n-lhs", type=int, default=72, help="Latin-Hypercube global search samples (default: 72)")
    ap.add_argument("--cv", action="store_true", help="run 5-fold record-grouped CV before full fit")
    args = ap.parse_args()

    tenhz.warmup()
    _k_match_sorted(np.array([10, 50], np.int64), np.array([10, 50], np.int64), 5)

    t0 = time.time()
    clean_records = load_records(args.data, args.dbs)
    print(f"Loaded {len(clean_records)} records from {args.dbs} (118 & 119 excluded).")

    print("Precomputing frozen MinEl & TiTen detections across clean + 4 noise SNRs...")
    cache = build_regimes(clean_records, noise_recs=(args.noise_recs or None), seed=42)
    print(f"Front-end precomputed in {time.time() - t0:.1f} s.\n")

    default_th = [FUSION["agree_s"], FUSION["w_single"], FUSION["lam"], FUSION["start_pen"], FUSION["win_s"], FUSION["thr"]]
    base_comp, base_scores = eval_cache(default_th, cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms)
    print(f"Default FUSION (+/-{args.tol_ms:g} ms, shift={args.shift_ms:g} ms):")
    print(f"  Composite = {base_comp:.3f} | " + " | ".join(f"{k}={v:.2f}" for k, v in base_scores.items()))

    if args.cv:
        n_rec = len(clean_records)
        perm = np.random.default_rng(42).permutation(n_rec)
        folds = np.array_split(perm, 5)
        cv_base, cv_opt = [], []
        print("\nRunning 5-fold record-grouped cross-validation...")
        for f_idx, val_idx in enumerate(folds, 1):
            print(f"-- Fold {f_idx}/5 --")
            val_ids = {clean_records[i][0] for i in val_idx}
            tr_cache = {k: [item for item in v if item[0] not in val_ids] for k, v in cache.items()}
            va_cache = {k: [item for item in v if item[0] in val_ids] for k, v in cache.items()}
            th_f = fit_fast(tr_cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms, n_lhs=args.n_lhs, seed=42 + f_idx)
            b_val, _ = eval_cache(default_th, va_cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms)
            o_val, _ = eval_cache(th_f, va_cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms)
            cv_base.append(b_val)
            cv_opt.append(o_val)
            print(f"  Fold {f_idx} val composite: {b_val:.3f} -> {o_val:.3f} ({o_val - b_val:+.3f})")
        print(f"Mean 5-fold val composite: {np.mean(cv_base):.3f} -> {np.mean(cv_opt):.3f} ({np.mean(cv_opt) - np.mean(cv_base):+.3f})")

    print("\nOptimizing Fusion post-processor parameters on full training set...")
    best_th = fit_fast(cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms, n_lhs=args.n_lhs, seed=42)
    opt_comp, opt_scores = eval_cache(best_th, cache, tol_ms=args.tol_ms, shift_ms=args.shift_ms)
    print(f"Optimized FUSION (+/-{args.tol_ms:g} ms, shift={args.shift_ms:g} ms):")
    print(f"  Composite = {opt_comp:.3f} | " + " | ".join(f"{k}={v:.2f}" for k, v in opt_scores.items()))

    keys = ["agree_s", "w_single", "lam", "start_pen", "win_s", "thr"]
    best_dict = dict(band=(7, 32), **{k: round(float(v), 4) for k, v in zip(keys, best_th)})
    print("\nDrop-in dictionary for src/tenhz/fusion.py:")
    print(f"FUSION = {best_dict}")


if __name__ == "__main__":
    main()