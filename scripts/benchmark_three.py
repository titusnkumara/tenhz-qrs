#!/usr/bin/env python3
"""Benchmark MinEl, TiTen and fusion (Numba versions) with the paper's evaluation rules.

    python scripts/benchmark_three.py --data ../tenhz-qrs/data --out results/three_v1

Rules (identical to tools/benchmark_neurokit.py; match() and the bootstrap are copied from it):
  * records:   <data>/<db>-mat/<rec>.txt (first lead, mV) and <rec>-ann.txt (1-based beats);
               the sampling rate comes from <db>-mat/_manifest.json
  * training:  MITDB (118 and 119 removed) and SVDB
  * test:      EDB, NSRDB, QTDB without the 28 excerpts of MITDB/SVDB recordings, TWADB,
               NSTDB (118e*/119e* only)
  * matching:  +/-150 ms (the paper's tolerance, default), one-to-one, closest pairs first;
               --tol-ms 150 75 40 (the default) scores several tolerances from the same detections
  * metrics:   TP/FP/FN summed over records before Se, PPV, F1 (pooled); 95% CI of the pooled
               test F1 from 2000 bootstrap resamples of records, stratified by database
  * signal:    stored as float32 and read back as float, as in the MATLAB evaluation

--offset-correct adds an ORACLE diagnostic: each record's detections are shifted by minus their
median offset to the reference (estimated within +/-150 ms) before scoring. It shows detection quality independent of which
point of the QRS is marked; it is not a result for the detectors as they would run.

Competitors: --methods ours nk2 wfdb baseline (see competitors.py). Speed: run with --workers 1 and
quote only that. Times are detection only
(file reading excluded), after the Numba kernels are compiled.

Outputs in --out (must be new or empty): per_record.csv, summary.csv, summary.md.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tenhz  # noqa: E402
import competitors  # noqa: E402

TRAIN = ["mitdb", "svdb"]
TEST = ["edb", "nsrdb", "qtdb", "twadb", "nstdb"]
HELDOUT = ["ltafdb"]       # held-out database: nothing was tuned on it; reported separately, never pooled with TEST
QTDB_OVERLAP = {  # QTDB excerpts of MITDB and SVDB recordings (QTDB documentation, Table 2)
    "sel100", "sel102", "sel103", "sel104", "sel114", "sel116", "sel117", "sel123", "sel213", "sel221",
    "sel223", "sel230", "sel231", "sel232", "sel233",
    "sel803", "sel808", "sel811", "sel820", "sel821", "sel840", "sel847", "sel853", "sel871", "sel872",
    "sel873", "sel883", "sel891"}
TOL_MS = [150.0, 75.0, 40.0]   # default: +/-150 ms is the paper's primary tolerance; 75 and 40 for the tighter columns
OFFSET_WINDOW_MS = 150.0  # window in which a record's median offset is estimated for --offset-correct
OURS = competitors.OURS
LABEL = competitors.LABEL


def mround(x):
    return int(np.sign(x) * np.floor(abs(x) + 0.5))


def records(data, dbs):
    """Same record lists as benchmark_neurokit.records, plus each record's sampling rate."""
    out = []
    for db in dbs:
        d = os.path.join(data, f"{db}-mat")
        mf = os.path.join(d, "_manifest.json")
        if not os.path.isdir(d) or not os.path.exists(mf):
            print(f"  {db}: {d} or its _manifest.json not found, skipped"); continue
        man = json.load(open(mf))
        names = sorted(f[:-4] for f in os.listdir(d)
                       if f.endswith(".txt") and not f.endswith(("-ann.txt", "-fs.txt")))
        if db == "nstdb":
            names = [n for n in names if n.startswith(("118e", "119e"))]
        if db == "mitdb":
            names = [n for n in names if n not in ("118", "119")]
        if db == "qtdb":
            names = [n for n in names if n not in QTDB_OVERLAP]
        for n in names:
            if man.get(n, {}).get("status") != "OK":
                print(f"  {db}/{n}: not OK in manifest, skipped"); continue
            out.append((db, n, float(man[n]["fs"])))
    return out


def load(data, db, rec):
    """ECG (float32, as the MATLAB evaluation stores it, then float) and 0-based reference beats."""
    base = os.path.join(data, f"{db}-mat", rec)
    ecg = pd.read_csv(base + ".txt", header=None, sep=r"\s+", usecols=[0], engine="c").to_numpy(dtype=float, copy=True).ravel()
    ecg[~np.isfinite(ecg)] = 0
    ann = pd.read_csv(base + "-ann.txt", header=None, sep=r"\s+", usecols=[0], engine="c").to_numpy(dtype=float, copy=True).ravel()
    ann = np.round(ann).astype(np.int64) - 1
    ann = np.sort(ann[(ann >= 0) & (ann < len(ecg))])
    return ecg.astype(np.float32).astype(float), ann


def match(det, gt, tol):
    """One-to-one matching within +/-tol samples, closest pairs first (as the paper)."""
    det = np.sort(np.asarray(det, np.int64)); gt = np.sort(np.asarray(gt, np.int64))
    lo = np.searchsorted(gt, det - tol, "left"); hi = np.searchsorted(gt, det + tol, "right")
    cnt = hi - lo
    I = np.repeat(np.arange(len(det)), cnt)
    J = np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    D = det[I] - gt[J]
    used_d = np.zeros(len(det), bool); used_g = np.zeros(len(gt), bool); tp = 0; off = []
    for k in np.argsort(np.abs(D), kind="stable"):
        i, j = I[k], J[k]
        if not used_d[i] and not used_g[j]:
            used_d[i] = used_g[j] = True; tp += 1; off.append(D[k])
    return tp, len(det) - tp, len(gt) - tp, (float(np.median(off)) if off else np.nan)


def match_pairs(det, gt, tol):
    """Same one-to-one closest-first matching as match(); returns the matched detections (sample values)
    and the matched indices into the sorted reference."""
    det = np.sort(np.asarray(det, np.int64)); gt = np.sort(np.asarray(gt, np.int64))
    lo = np.searchsorted(gt, det - tol, "left"); hi = np.searchsorted(gt, det + tol, "right")
    cnt = hi - lo
    I = np.repeat(np.arange(len(det)), cnt)
    J = np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    D = det[I] - gt[J]
    used_d = np.zeros(len(det), bool); used_g = np.zeros(len(gt), bool); md, mg = [], []
    for k in np.argsort(np.abs(D), kind="stable"):
        i, j = I[k], J[k]
        if not used_d[i] and not used_g[j]:
            used_d[i] = used_g[j] = True; md.append(det[i]); mg.append(j)
    return np.asarray(md, np.int64), np.asarray(mg, np.int64), gt


RR_BINS = 401          # |RR error| histogram: 1 ms bins 0..399, last bin = 400 ms or more


def rr_error_stats(det, ann, fs):
    """RR-interval error over consecutive matched beats (matching at +/-150 ms):
    dRR = (d[i+1] - d[i]) - (r[i+1] - r[i]). A constant offset between the detector's marker and the
    annotation cancels exactly. Returns n, sum|dRR| (ms), sum dRR^2 (ms^2) and a 1 ms histogram of |dRR|."""
    md, mg, gt = match_pairs(det, ann, mround(150 / 1000 * fs))
    o = np.argsort(mg)
    md, mg = md[o], mg[o]
    k = np.flatnonzero(np.diff(mg) == 1)                      # reference beats i and i+1 both matched
    if len(k) == 0:
        return 0, 0.0, 0.0, np.zeros(RR_BINS, np.int64)
    ms = 1000.0 * ((md[k + 1] - md[k]) - (gt[mg[k + 1]] - gt[mg[k]])) / fs
    a = np.abs(ms)
    return int(len(k)), float(a.sum()), float((ms ** 2).sum()), np.bincount(np.minimum(a.astype(int), RR_BINS - 1),
                                                                            minlength=RR_BINS)


def _init():
    warnings.filterwarnings("ignore")


def run_detector(method, ecg, fs):
    if method in OURS:
        return tenhz.qrs_locations(ecg, fs, method)
    return competitors.detect(method, ecg, fs)


def task(args):
    data, db, rec, fs, method, tols, correct, save_dir = args
    try:
        ecg, ann = load(data, db, rec)
    except Exception as e:
        return [dict(db=db, record=rec, method=method, variant=v, tol_ms=t, tp=0, fp=0, fn=0, seconds=0.0,
                     hours=0.0, offset_ms=np.nan, error=f"could not load record: {type(e).__name__}: {e}"[:200], note="")
                for v in (["raw", "offset_corrected"] if correct else ["raw"]) for t in tols]
    t0 = time.perf_counter()
    competitors.take_note()
    try:
        det = run_detector(method, ecg, fs); err = ""
    except Exception as e:  # a failing detector scores every beat as missed
        det = np.zeros(0, np.int64); err = f"{type(e).__name__}: {e}"[:200]
    note = competitors.take_note()
    if save_dir and not err:
        os.makedirs(os.path.join(save_dir, db), exist_ok=True)
        np.save(os.path.join(save_dir, db, f"{rec}__{method}.npy"), np.asarray(det, np.int32))
    rr = rr_error_stats(det, ann, fs) if not err else (0, 0.0, 0.0, np.zeros(RR_BINS, np.int64))
    sec = time.perf_counter() - t0
    out = []
    variants = [("raw", det)]
    if correct:
        # Oracle diagnostic: shift this record's detections by minus their median offset to the
        # reference (estimated at the widest tolerance), then score again. Uses the reference, so it
        # measures detection independent of which QRS point is marked; it is not a deployable detector.
        off = match(det, ann, mround(OFFSET_WINDOW_MS / 1000 * fs))[3]
        shift = 0 if np.isnan(off) else int(round(off))
        variants.append(("offset_corrected", np.clip(det - shift, 0, None)))
    for v, d in variants:
        for t in tols:
            tp, fp, fn, off = match(d, ann, mround(t / 1000 * fs))
            out.append(dict(db=db, record=rec, method=method, variant=v, tol_ms=t, tp=tp, fp=fp, fn=fn,
                            seconds=sec, hours=len(ecg) / fs / 3600, offset_ms=1000 * off / fs, error=err, note=note))
    out[0].update(rr_n=rr[0], rr_sumabs=rr[1], rr_sumsq=rr[2], rr_hist=",".join(map(str, rr[3])))   # once per record/detector
    return out


def pooled(g):
    tp, fp, fn = g.tp.sum(), g.fp.sum(), g.fn.sum()
    return 100 * tp / max(tp + fn, 1), 100 * tp / max(tp + fp, 1), 100 * 2 * tp / max(2 * tp + fp + fn, 1)


def boot_ci(g, B=2000, seed=1):
    rng = np.random.default_rng(seed)
    f = []
    by = [x for _, x in g.groupby("db")]
    for _ in range(B):
        s = pd.concat([x.iloc[rng.integers(0, len(x), len(x))] for x in by])
        f.append(pooled(s)[2])
    return np.percentile(f, [2.5, 97.5])


def holm(p):
    """Holm-Bonferroni adjusted p-values."""
    p = np.asarray(p, float); o = np.argsort(p); n = len(p)
    adj = np.empty(n); run = 0.0
    for rank, i in enumerate(o):
        run = max(run, (n - rank) * p[i]); adj[i] = min(1.0, run)
    return adj


def paired(R, out, variant, tol, B=2000, seed=1, dbs=None, label="test"):
    """paired_<variant>_<tol>[_<label>].md/csv: each of our detectors against every other detector on the
    same records: difference in pooled F1 with a paired bootstrap CI (the same record resamples for all
    detectors, stratified by database) and a Wilcoxon signed-rank test on per-record F1 with Holm correction.

    dbs: database list to pool over (default TEST). label: 'test' for the primary tables, or a database name
    (e.g. 'ltafdb') for a held-out table written to paired_<variant>_<tol>_<label>.md/csv."""
    from scipy.stats import wilcoxon
    if dbs is None:
        dbs = TEST
    g = R[(R.variant == variant) & (R.tol_ms == tol) & R.db.isin(dbs)]
    if g.empty:
        return
    keys = g[["db", "record"]].drop_duplicates().sort_values(["db", "record"]).reset_index(drop=True)
    pos = {(d, r): i for i, (d, r) in enumerate(zip(keys.db, keys.record))}
    methods = [m for m in g.method.unique() if g[g.method == m].drop_duplicates(["db", "record"]).shape[0] == len(keys)]
    dropped = [m for m in g.method.unique() if m not in methods]
    TP = np.zeros((len(keys), len(methods))); FP = TP.copy(); FN = TP.copy()
    col = {m: j for j, m in enumerate(methods)}
    for r in g[g.method.isin(methods)].itertuples():
        i, j = pos[(r.db, r.record)], col[r.method]
        TP[i, j], FP[i, j], FN[i, j] = r.tp, r.fp, r.fn
    f1 = lambda tp, fp, fn: 100 * 2 * tp / np.maximum(2 * tp + fp + fn, 1)
    point = f1(TP.sum(0), FP.sum(0), FN.sum(0))
    den = 2 * TP + FP + FN
    rec = np.where(den > 0, 2 * TP / np.maximum(den, 1), np.nan)
    strata = [np.flatnonzero((keys.db == d).to_numpy()) for d in keys.db.unique()]
    rng = np.random.default_rng(seed)
    boot = np.empty((B, len(methods)))
    for b in range(B):
        idx = np.concatenate([rng.choice(ix, len(ix)) for ix in strata])
        boot[b] = f1(TP[idx].sum(0), FP[idx].sum(0), FN[idx].sum(0))
    refs = [m for m in OURS if m in col]
    others = sorted(methods, key=lambda m: -point[col[m]])
    _kind = "offset-corrected (oracle)" if variant != "raw" else "raw"
    if label == "test":
        _hdr = f"# Paired comparison: {_kind} F1 at +/-{tol:g} ms"
        _rec = f"{len(keys)} test records"
    else:
        _hdr = f"# Paired comparison: {_kind} F1 at +/-{tol:g} ms [{label.upper()} — held out, nothing tuned on it]"
        _rec = f"{len(keys)} held-out {label.upper()} records"
    rows, md = [], [_hdr, "",
                    f"{_rec}. Difference = reference minus comparator, in F1 points, pooled over records; "
                    "CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified "
                    "by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators "
                    "of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so "
                    "long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is "
                    "marked only when both agree: CI excludes 0 and p(Holm) < 0.05.", ""]
    for ref in refs:
        j = col[ref]; comps = [m for m in others if m != ref]
        pv, cnt = [], []
        for m in comps:
            d = rec[:, j] - rec[:, col[m]]; d = d[~np.isnan(d)]
            pv.append(1.0 if (len(d) == 0 or np.all(d == 0)) else wilcoxon(d).pvalue)
            cnt.append((int((d > 1e-9).sum()), int((d < -1e-9).sum()), int((np.abs(d) <= 1e-9).sum())))
        pa = holm(pv)
        md += [f"## {LABEL.get(ref, ref)} against each other detector (reference F1 {point[j]:.2f})", "",
               "| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |",
               "|---|---|---|---|---|---|"]
        for m, p_raw, p_adj, (nb_, nw_, nt_) in zip(comps, pv, pa, cnt):
            dd = point[j] - point[col[m]]
            lo, hi = np.percentile(boot[:, j] - boot[:, col[m]], [2.5, 97.5])
            mark = ("reference better" if (lo > 0 and p_adj < 0.05) else
                    "reference worse" if (hi < 0 and p_adj < 0.05) else "")
            md.append(f"| {LABEL.get(m, m)} | {point[col[m]]:.2f} | {dd:+.2f} [{lo:+.2f}, {hi:+.2f}] | {p_adj:.3g} | {nb_} / {nw_} / {nt_} | {mark} |")
            rows.append(dict(variant=variant, tol_ms=tol, reference=ref, comparator=m, ref_F1=point[j],
                             comp_F1=point[col[m]], diff=dd, ci_lo=lo, ci_hi=hi, p_raw=p_raw, p_holm=p_adj,
                             rec_better=nb_, rec_worse=nw_, rec_tied=nt_, verdict=mark))
        md.append("")
    if dropped:
        md.append(f"Not compared (missing records): {', '.join(dropped)}")
    tag = f"{'oracle' if variant != 'raw' else 'raw'}_{tol:g}"
    if label != "test":
        tag += f"_{label}"
    pd.DataFrame(rows).to_csv(os.path.join(out, f"paired_{tag}.csv"), index=False)
    open(os.path.join(out, f"paired_{tag}.md"), "w").write("\n".join(md) + "\n")


def rr_error(R, out, dbs=None, label="test"):
    """rr_error[_<label>].md/csv: RR-interval error over consecutive matched beats, per detector."""
    if dbs is None:
        dbs = TEST
    if "rr_hist" not in R or not (R.rr_hist.astype(str) != "").any():
        return
    S = R[R.db.isin(dbs) & (R.rr_hist.astype(str) != "")].drop_duplicates(["db", "record", "method"])
    if S.empty:
        return
    def agg(g):
        h = np.sum([np.array(s.split(","), dtype=np.int64) for s in g.rr_hist], axis=0)
        n = int(h.sum())
        if n == 0:
            return dict(n=0, median=np.nan, mean=np.nan, rmse=np.nan, ge20=np.nan, ge40=np.nan, ge75=np.nan)
        cum = np.cumsum(h)
        med = float(np.searchsorted(cum, 0.5 * n)) + 0.5                  # 1 ms bins: bin centre
        sa, sq = pd.to_numeric(g.rr_sumabs).sum(), pd.to_numeric(g.rr_sumsq).sum()
        return dict(n=n, median=med, mean=sa / n, rmse=float(np.sqrt(sq / n)),
                    ge20=100 * h[20:].sum() / n, ge40=100 * h[40:].sum() / n, ge75=100 * h[75:].sum() / n)
    methods = list(S.method.unique())
    base = R[(R.variant == "raw") & (R.tol_ms == 150) & R.db.isin(dbs)]
    rows = []
    for m in methods:
        g = S[S.method == m]
        r = dict(detector=LABEL.get(m, m), **agg(g))
        for db in dbs:
            x = g[g.db == db]
            r[f"median_{db}"] = agg(x)["median"] if len(x) else np.nan
        b_ = base[base.method == m]
        r["Se"], r["PPV"] = (pooled(b_)[0], pooled(b_)[1]) if len(b_) else (np.nan, np.nan)
        rows.append(r)
    T = pd.DataFrame(rows).sort_values("median")
    tag = "rr_error" if label == "test" else f"rr_error_{label}"
    T.to_csv(os.path.join(out, f"{tag}.csv"), index=False)
    f = lambda v, d=1: "-" if pd.isna(v) else f"{v:.{d}f}"
    title = ("# RR-interval error over consecutive matched beats (test databases)" if label == "test"
             else f"# RR-interval error over consecutive matched beats (held out: {label.upper()})")
    lines = [title, "",
             "dRR = (detected interval) - (reference interval) for two consecutive reference beats that were both matched "
             "(+/-150 ms, one-to-one). A constant offset between a detector's marker and the annotation cancels exactly, "
             "so this measures interval timing, which is what rhythm and HRV analysis use; it is not an oracle. "
             "It does not cancel a change of offset from one beat to the next (morphology-dependent marker position), "
             "and it is conditional on detection: read it next to Se and PPV. Median is read from 1 ms bins.", "",
             "| Detector | median abs dRR (ms) | mean abs (ms) | RMSE (ms) | >=20 ms (%) | >=40 ms (%) | >=75 ms (%) | "
             + " | ".join(f"median {d.upper()}" for d in dbs) + " | pairs | Se | PPV |",
             "|---|---|---|---|---|---|---|" + "---|" * len(dbs) + "---|---|---|"]
    for _, r in T.iterrows():
        lines.append(f"| {r.detector} | {f(r['median'])} | {f(r['mean'])} | {f(r['rmse'])} | {f(r.ge20)} | {f(r.ge40)} | "
                     f"{f(r.ge75)} | " + " | ".join(f(r[f'median_{d}']) for d in dbs) + f" | {int(r.n)} | {f(r.Se, 2)} | {f(r.PPV, 2)} |")
    open(os.path.join(out, f"{tag}.md"), "w").write("\n".join(lines) + "\n")


def report(R, out, boot):
    all_rows, md = [], []
    for v in ["raw", "offset_corrected"]:
        for tol in sorted(R.tol_ms.unique(), reverse=True):
            g = R[(R.tol_ms == tol) & (R.variant == v)]
            if len(g):
                md += report_one(g, tol, boot, all_rows, v)
    pd.DataFrame(all_rows).to_csv(os.path.join(out, "summary.csv"), index=False)
    headline(all_rows, out, HEADLINE_TOL[0])
    rr_error(R, out)
    for db in HELDOUT:
        if R.db.eq(db).any():
            rr_error(R, out, dbs=[db], label=db)
    if PAIRED[0]:
        for v, t in (("offset_corrected", HEADLINE_TOL[0]), ("raw", 150.0)):
            paired(R, out, v, t)
            for db in HELDOUT:
                if R.db.eq(db).any():
                    paired(R, out, v, t, dbs=[db], label=db)
    md += ["", f"tenhz {tenhz.__version__}; numpy {np.__version__}."]
    open(os.path.join(out, "summary.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))

def heldout_table(R, db, methods, boot):
    """Markdown table for one held-out database: pooled over its records, CI from record resampling."""
    g = R[R.db == db]
    recs = g.drop_duplicates(["db", "record"])
    lines = ["", f"## Held-out: {db.upper()} ({len(recs)} records, {recs.hours.sum():.0f} h; nothing was tuned on it)", "",
             "| Detector | F1 [95% CI] | Se | PPV | marker offset (ms) | s per h of ECG | failed |", "|---|---|---|---|---|---|---|"]
    rows = []
    for m in methods:
        x = g[g.method == m]
        if x.empty:
            continue
        se, ppv, f1 = pooled(x)
        ci = boot_ci(x) if boot else (np.nan, np.nan)
        o = (x.offset_ms * x.tp).sum() / max(x.tp[x.offset_ms.notna()].sum(), 1) if x.offset_ms.notna().any() else np.nan
        rows.append((f1, m, se, ppv, ci, o, x.seconds.sum() / max(x.hours.sum(), 1e-9), int((x.error != "").sum())))
    ours = [r for r in rows if r[1] in OURS]
    rest = sorted([r for r in rows if r[1] not in OURS], key=lambda r: -r[0])
    for f1, m, se, ppv, ci, o, sph, nf in ours + rest:
        cs = "" if np.isnan(ci[0]) else f" [{ci[0]:.2f}, {ci[1]:.2f}]"
        lines.append(f"| {LABEL.get(m, m)} | {f1:.2f}{cs} | {se:.2f} | {ppv:.2f} | "
                     f"{'-' if np.isnan(o) else format(o, '+.0f')} | {sph:.4f} | {nf} |")
    return lines + [""]


HEADLINE_TOL = [75.0]
PAIRED = [True]


def headline(all_rows, out, tol):
    """headline.md/csv: one row per detector, ranked by offset-corrected F1 at +/-tol ms, next to raw +/-150 ms."""
    A = pd.DataFrame(all_rows)
    c = A[(A.variant == "offset_corrected") & (A.tol_ms == tol)]
    r = A[(A.variant == "raw") & (A.tol_ms == 150)]
    if c.empty or r.empty or c.test_F1.isna().all():
        return                                   # nothing to rank (e.g. a run with training databases only)
    H = c.merge(r[["detector", "test_F1", "offset_ms"]].rename(columns={"test_F1": "raw150_F1", "offset_ms": "raw_offset_ms"}),
                on="detector")
    ours = [LABEL[m] for m in OURS]
    HELD = [d for d in HELDOUT if d in H and H[d].notna().any()]
    H = H.assign(_a=[0 if d in ours else 1 for d in H.detector]).sort_values(
        ["_a", "test_F1"], ascending=[True, False]).drop(columns="_a")
    H.to_csv(os.path.join(out, "headline.csv"), index=False)
    f = lambda v: "-" if pd.isna(v) else f"{v:.2f}"
    lines = [f"# Headline: offset-corrected (oracle) F1 at +/-{tol:g} ms, with raw +/-150 ms alongside", "",
             f"Offset-corrected: each record's detections are shifted by minus their median offset to the reference "
             f"(estimated within +/-150 ms) before matching at +/-{tol:g} ms. It uses the reference, applies to every "
             "detector identically, and measures detection independent of which QRS point is marked; it is not a "
             "result for a detector as deployed. Raw +/-150 ms is the conventional score.", "",
             f"| Detector | F1 oracle +/-{tol:g} ms [95% CI] | Se | PPV | EDB | NSRDB | QTDB | TWADB | NSTDB | mean of test-DB F1 | "
             + "".join(f"{d.upper()} (held out) | " for d in HELD) + "raw +/-150 ms F1 | raw marker offset (ms) | s per h of ECG | failed |",
             "|---|---|---|---|---|---|---|---|---|---|" + "---|" * len(HELD) + "---|---|---|---|"]
    for _, x in H.iterrows():
        ci = "" if pd.isna(x.ci_lo) else f" [{x.ci_lo:.2f}, {x.ci_hi:.2f}]"
        lines.append(f"| {x.detector} | {x.test_F1:.2f}{ci} | {x.test_Se:.2f} | {x.test_PPV:.2f} | "
                     + " | ".join(f(x[d]) for d in TEST) + f" | {f(x.per_db_mean_F1)} | " + "".join(f"{f(x[d])} | " for d in HELD)
                     + f"{x.raw150_F1:.2f} | "
                     f"{x.raw_offset_ms:+.0f} | {x.s_per_h:.4f} | {int(x.failed)} |")
    lines += ["", "s per h comes from the run that produced these rows: indicative only unless run with --workers 1."]
    open(os.path.join(out, "headline.md"), "w").write("\n".join(lines) + "\n")


def report_one(R, tol, boot, all_rows, variant="raw"):
    rows = []
    methods = [m for m in R.method.unique()]
    def _f1(g):
        return pooled(g)[2] if len(g) else -1.0
    eval_dbs = TEST if R.db.isin(TEST).any() else HELDOUT
    f1_of = {m: _f1(R[(R.method == m) & R.db.isin(eval_dbs)]) for m in methods}
    methods = [m for m in OURS if m in methods] + sorted([m for m in methods if m not in OURS],
                                                         key=lambda m: -f1_of[m])
    for m in methods:
        g = R[R.method == m]
        if g.empty:
            continue
        te = g[g.db.isin(TEST)]
        se, ppv, f1 = pooled(te) if len(te) else (np.nan, np.nan, np.nan)   # empty test set: show "-"
        ci = boot_ci(te) if (boot and len(te)) else (np.nan, np.nan)
        row = dict(variant=variant, tol_ms=tol, detector=LABEL.get(m, m), test_F1=f1, ci_lo=ci[0], ci_hi=ci[1], test_Se=se, test_PPV=ppv)
        for db in TEST + TRAIN:
            x = g[g.db == db]
            row[db] = pooled(x)[2] if len(x) else np.nan
        for db in HELDOUT:
            x = g[g.db == db]
            row[db] = pooled(x)[2] if len(x) else np.nan
        row["per_db_mean_F1"] = np.nanmean([row[db] for db in TEST]) if any(
            len(g[g.db == db]) for db in TEST) else np.nan
        row["s_per_h"] = te.seconds.sum() / te.hours.sum() if len(te) and te.hours.sum() > 0 else np.nan
        o = (te.offset_ms * te.tp).sum() / max(te.tp[te.offset_ms.notna()].sum(), 1)
        row["offset_ms"] = o if te.offset_ms.notna().any() else np.nan
        row["failed"] = int((g.error != "").sum())
        row["records"] = len(g)
        rows.append(row)
    all_rows += rows
    n_te = R[R.db.isin(TEST)].drop_duplicates(["db", "record"])
    missing = [db.upper() for db in TEST if not (n_te.db == db).any()]
    title = ("" if variant == "raw" else " [OFFSET-CORRECTED, ORACLE: per-record median offset removed using the "
             "reference; diagnostic only]")
    lines = [f"# 10 Hz detectors and competitors, matching tolerance +/-{tol:g} ms{title}", "",
             f"Test set: {len(n_te)} records ({', '.join(f'{db.upper()} {sum(n_te.db == db)}' for db in TEST)}). "
             f"Matching +/-{tol:g} ms, one-to-one, closest first; pooled over records. "
             "MinEl and TiTen were tuned on MITDB and SVDB, so their training columns are in-sample. "]
    held = [db for db in HELDOUT if R.db.eq(db).any()]
    if not len(n_te) and held:
        lines += ["", f"**This run holds only the held-out database(s) {', '.join(h.upper() for h in held)}: the test-set columns "
                  "are empty by design; see the held-out table below.**"]
    elif not len(n_te):
        lines += ["", "**No test-set records in this run (only training databases, MITDB/SVDB): the test columns are empty "
                  "and the MITDB/SVDB columns are in-sample for MinEl, TiTen and fusion.**"]
    elif missing:
        lines += ["", f"**WARNING: no records for {', '.join(missing)}. This is not the full 274-record test set; "
                  "pooled numbers are not comparable with the paper's.**"]
    if len(n_te) or not held:                       # a run with only a held-out database has no test-set table
        lines += ["",
                 "| Detector | Test F1 [95% CI] | Se | PPV | EDB | NSRDB | QTDB | TWADB | NSTDB | "
                 "mean of test-DB F1 | MITDB | SVDB | marker offset (ms) | s per h of ECG | failed |",
                 "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows:
            ci = "" if np.isnan(r["ci_lo"]) else f" [{r['ci_lo']:.2f}, {r['ci_hi']:.2f}]"
            f = lambda v: "-" if np.isnan(v) else f"{v:.2f}"
            lines.append(f"| {r['detector']} | {f(r['test_F1'])}{ci} | {f(r['test_Se'])} | {f(r['test_PPV'])} | "
                         + " | ".join(f(r[db]) for db in TEST) + f" | {f(r['per_db_mean_F1'])} | "
                         + " | ".join(f(r[db]) for db in TRAIN)
                         + f" | {'-' if np.isnan(r['offset_ms']) else format(r['offset_ms'], '+.0f')} | "
                         + f"{'-' if np.isnan(r['s_per_h']) else format(r['s_per_h'], '.4f')} | {r['failed']} |")
        lines += ["", "Pooled test F1 is dominated by NSRDB (about two thirds of test beats); the mean of the "
                  "per-database F1 values weights every database equally.",
                  "Marker offset: median of (detection - reference beat) over matched beats, test set."]
    for db in held:
        lines += heldout_table(R, db, methods, boot)
    N = R[R.note != ""].drop_duplicates(["db", "record", "method"])
    if len(N):
        lines += ["", "Fallbacks used (the detector's own library raised an exception on these records):", ""]
        for (m, why), g in N.groupby(["method", "note"]):
            lines.append(f"- {m}: {len(g)} record(s), {why} ({', '.join(f'{r.db}/{r.record}' for r in g.itertuples())})")
    E = R[R.error != ""]
    if len(E):
        lines += ["", "Failures (scored as all beats missed):", ""]
        lines += [f"- {e.method} on {e.db}/{e.record}: {e.error}" for e in E.itertuples()]
    return lines + [""]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("TENHZ_DATA") or "data")
    ap.add_argument("--out", required=True, help="new folder for the results")
    ap.add_argument("--dbs", nargs="*", default=TRAIN + TEST)
    ap.add_argument("--tol-ms", nargs="+", type=float, default=TOL_MS,
                    help="matching tolerance(s) in ms (default 150 75 40; the paper's primary tolerance is 150)")
    ap.add_argument("--methods", nargs="+", default=["ours"],
                    help="detector and/or group names: ours (default), nk2, wfdb, pyecg, baseline, all; "
                         "e.g. --methods ours nk2 wfdb baseline")
    ap.add_argument("--workers", type=int, default=1, help="1 for clean timing; more only for a quick look")
    ap.add_argument("--offset-correct", action=argparse.BooleanOptionalAction, default=True,
                    help="also score after removing each record's median offset (oracle diagnostic, see "
                         "docstring); on by default, --no-offset-correct to skip")
    ap.add_argument("--exclude", nargs="+", default=[], help="detector names to leave out of --methods")
    ap.add_argument("--merge", nargs="+", default=None, metavar="DIR",
                    help="no detection: combine the per_record.csv of earlier runs (e.g. one per database) "
                         "into --out and write the summary tables")
    ap.add_argument("--prefer-last", action="store_true",
                    help="with --merge: where a record/detector appears in several folders keep the row from "
                         "the last folder (to replace a re-run detector)")
    ap.add_argument("--headline-tol", type=float, default=75.0,
                    help="tolerance (ms) of the offset-corrected column in headline.md (default 75)")
    ap.add_argument("--drop-dbs", nargs="+", default=[], choices=TRAIN + TEST + HELDOUT,
                    help="leave these databases out of the tables, e.g. --drop-dbs nstdb for a fusion-clean view")
    ap.add_argument("--save-detections", action=argparse.BooleanOptionalAction, default=True,
                    help="save each detector's beat indices to <out>/detections/<db>/<rec>__<method>.npy so later "
                         "analyses need no re-detection (on by default)")
    ap.add_argument("--no-paired", action="store_true", help="skip the paired comparison tables")
    ap.add_argument("--no-boot", action="store_true", help="skip the bootstrap CI")
    a = ap.parse_args()
    HEADLINE_TOL[0] = a.headline_tol
    PAIRED[0] = not (a.no_paired or a.no_boot)
    a.methods = [m for m in competitors.expand(a.methods) if m not in a.exclude]
    if os.path.exists(a.out) and os.listdir(a.out):
        sys.exit(f"{a.out} exists and is not empty: give a new --out folder")
    os.makedirs(a.out, exist_ok=True)
    if a.merge:
        parts = []
        for d in a.merge:
            f = os.path.join(d, "per_record.csv")
            if not os.path.exists(f):
                sys.exit(f"{f} not found")
            parts.append(pd.read_csv(f, dtype={"record": str}, keep_default_na=False))
        R = pd.concat(parts, ignore_index=True)
        for c_ in ("rr_n", "rr_sumabs", "rr_sumsq", "rr_hist"):
            if c_ not in R:
                R[c_] = ""
        if "note" not in R:
            R["note"] = ""
        R["note"] = R["note"].fillna("").astype(str)
        key = ["db", "record", "method", "variant", "tol_ms"]
        if a.prefer_last:
            n0 = len(R)
            R = R.drop_duplicates(key, keep="last").reset_index(drop=True)
            print(f"--prefer-last: replaced {n0 - len(R)} rows from earlier folders")
        dup = R[R.duplicated(key, keep=False)]
        if len(dup):
            sys.exit(f"{len(dup)} rows appear in more than one folder, e.g. "
                     f"{dup.iloc[0][['db', 'record', 'method']].tolist()}: merge disjoint runs only")
        for c in ["tp", "fp", "fn"]:
            R[c] = R[c].astype(int)
        for c in ["seconds", "hours", "offset_ms"]:
            R[c] = pd.to_numeric(R[c], errors="coerce")
        R["error"] = R["error"].fillna("").astype(str)
        R = R[~R.db.isin(a.drop_dbs)].reset_index(drop=True)
        print(f"merged {len(a.merge)} folders: {R.groupby(['db']).record.nunique().to_dict()} records per database; "
              f"{R.method.nunique()} detectors")
        R.to_csv(os.path.join(a.out, "per_record.csv"), index=False)
        report(R, a.out, not a.no_boot)
        return
    save_dir = os.path.join(a.out, "detections") if a.save_detections else None
    recs = records(a.data, a.dbs)
    if not recs:
        sys.exit("no records found")
    todo = [(a.data, db, r, fs, m, tuple(a.tol_ms), a.offset_correct, save_dir) for (db, r, fs) in recs for m in a.methods]
    # longest records first so the 24 h NSRDB files do not finish last
    todo.sort(key=lambda t: (t[1] != "nsrdb", t[1] != "edb"))
    print(f"{len(recs)} records x {len(a.methods)} detectors = {len(todo)} runs; {a.workers} worker(s)", flush=True)
    if a.workers > 1:
        print("note: with several workers the 's per h' column is indicative only; quote speed from "
              "benchmark_timing.py or a --workers 1 run", flush=True)
    rows = []; t0 = time.time()
    if a.workers == 1:
        _init()
        for i, t in enumerate(todo, 1):
            res = task(t); rows += res; r = res[0]
            print(f"[{i}/{len(todo)}] {r['method']:6s} {r['db']}/{r['record']}: "
                  f"F1 {100 * 2 * r['tp'] / max(2 * r['tp'] + r['fp'] + r['fn'], 1):6.2f}  {r['seconds']:7.3f} s"
                  + (f"  ERROR {r['error']}" if r["error"] else ""), flush=True)
    else:
        # Compile Numba kernels once in the parent process.
        warnings.filterwarnings("ignore")
        print("Warming up Numba kernels...", flush=True)
        tenhz.warmup()
        print("Numba warmup complete. Starting workers.", flush=True)

        import multiprocessing as mp

        ctx = mp.get_context("fork")

        with ProcessPoolExecutor(a.workers, initializer=_init, mp_context=ctx) as ex:
            futs = [ex.submit(task, t) for t in todo]
            for i, f in enumerate(as_completed(futs), 1):
                try:
                    res = f.result()
                except Exception:
                    print(traceback.format_exc()); continue
                rows += res
                if i % 25 == 0 or i == len(todo):
                    print(f"[{i}/{len(todo)}]", flush=True)
    print(f"\nRan {len(todo)} in {(time.time() - t0) / 60:.1f} min\n")
    R = pd.DataFrame(rows)
    R = R[~R.db.isin(a.drop_dbs)].reset_index(drop=True)
    R.to_csv(os.path.join(a.out, "per_record.csv"), index=False)
    report(R, a.out, not a.no_boot)


if __name__ == "__main__":
    main()