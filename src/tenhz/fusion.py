"""Fusion: MinEl + TiTen on one shared 7-32 Hz stream, an RR-plausibility
decoder (Viterbi) and an adaptive per-window switch.
"""
import numpy as np

from ._kernels import _k_local_median, _k_pair_flags, _k_viterbi
from ._signal import butter, filtfilt_matlab
from .detectors import RECOVER_S, minel_detect, minmax_stream, snap_to_r, titen_detect
from .params import PARAMS

# 2.1 decoder (agreed-beat median RR, no amplitude term): run_chain(..., params=FUSION_2_1) reproduces 2.1.x exactly.
FUSION_2_1 = {
    "band": (7, 32),
    "agree_s": 0.0832,
    "w_single": 0.4358,
    "lam": 0.8227,
    "start_pen": 1.1415,
    "win_s": 16.5281,
    "thr": 0.0183,
    "two_pass": False,
    "amp_lo": 0.0,
    "amp_hi": 0.0,
}

# 2.2 decoder (default). two_pass: the local median RR comes from a first decoding instead of from agreed MinEl beats
# only (in disagreement stretches the agreed-beat median is about twice the true RR, so true beats were skipped).
# amp_lo / amp_hi: score penalty per unit |log| for candidates smaller / larger than the local beat amplitude.
# All eight values tuned jointly on MITDB (without 118, 119) + SVDB, clean and with synthetic noise
# (scripts/optimize_fusion.py regimes, seed 42), pooled F1 at +/-150 ms; no test or held-out database was used.
FUSION = {
    "band": (7, 32),
    "agree_s": 0.07348387325273786,
    "w_single": 0.7298600541334451,
    "lam": 0.6519013192580632,
    "start_pen": 1.9657285125322124,
    "win_s": 16.420113268298174,
    "thr": 0.021988778315839785,
    "two_pass": True,
    "amp_lo": 0.7213747196579091,
    "amp_hi": 0.061795243350803886,
}


def pair_flags(a, b, w):
    """Which entries of a and b fusion_explore.pair matches (one-to-one within
    w samples, closest first, ties by a index then b index)."""
    a = np.asarray(a, np.int64); b = np.asarray(b, np.int64)
    lo = np.searchsorted(b, a - w, "left"); hi = np.searchsorted(b, a + w, "right")
    cnt = hi - lo
    I = np.repeat(np.arange(len(a)), cnt)
    J = np.repeat(lo, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    order = np.argsort(np.abs(a[I] - b[J]), kind="stable")
    return _k_pair_flags(I, J, order, len(a), len(b))


def viterbi(m, t, pm, pt, fs, w_single, lam, start_pen, ref=None, sig=None, amp_lo=0.0, amp_hi=0.0):
    """fusion_explore.viterbi; the logs come from one vectorised np.log call.

    ref    : beats the local median RR (and amplitude) is taken from; None = agreed MinEl beats (2.1 behaviour).
    sig    : polarity-corrected bandpassed ECG (native rate). With amp_lo/amp_hi > 0 a candidate's score is reduced by
             amp_lo * max(0, -log(a / a_ref)) + amp_hi * max(0, log(a / a_ref)), where a is its R amplitude and a_ref
             the local median amplitude of the reference beats."""
    cand = np.r_[m[pm], m[~pm], t[~pt]].astype(np.int64)
    sc = np.r_[np.ones(pm.sum()), np.full((~pm).sum(), w_single), np.full((~pt).sum(), w_single)]
    o = np.argsort(cand, kind="stable"); c, sc = cand[o], sc[o]
    if len(c) < 3:
        return c
    ag = np.sort(m[pm]) if ref is None else np.sort(np.asarray(ref, np.int64))
    rr = np.diff(ag) if len(ag) > 2 else np.array([0.8 * fs])
    k = np.searchsorted(ag, c).clip(1, max(1, len(rr)))
    med = _k_local_median(rr.astype(np.float64), k.astype(np.int64))
    if sig is not None and (amp_lo > 0 or amp_hi > 0) and len(ag) > 2:
        a_ref = _k_local_median(np.maximum(sig[ag], 1e-9).astype(np.float64), k.astype(np.int64))
        lr = np.log(np.maximum(sig[c], 1e-9) / np.maximum(a_ref, 1e-9))
        sc = sc - amp_lo * np.maximum(-lr, 0.0) - amp_hi * np.maximum(lr, 0.0)
    lo, hi = int(0.25 * fs), int(2.5 * fs)
    ptr = np.minimum(np.searchsorted(c, c - hi, "left"), np.arange(len(c)))
    cnt = np.arange(len(c)) - ptr
    I = np.repeat(np.arange(len(c)), cnt)
    J = np.repeat(ptr, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    dd = c[I] - c[J]
    ok = dd >= lo
    logv = np.zeros(len(dd))
    logv[ok] = np.log(dd[ok] / med[I[ok]])
    return _k_viterbi(c, sc, logv, lo, hi, float(lam), float(start_pen))


def decode(m, t, pm, pt, fs, P, sig=None):
    """Viterbi decoding with the options in P (FUSION): optional amplitude term, optional second pass whose local
    median RR comes from the first pass."""
    kw = dict(sig=sig, amp_lo=P.get("amp_lo", 0.0), amp_hi=P.get("amp_hi", 0.0))
    V = viterbi(m, t, pm, pt, fs, P["w_single"], P["lam"], P["start_pen"], **kw)
    if P.get("two_pass", False):
        V = viterbi(m, t, pm, pt, fs, P["w_single"], P["lam"], P["start_pen"], ref=V, **kw)
    return V


def adaptive(m, t, fused, pm, pt, fs, win_s, thr):
    """fusion_explore.adaptive, unchanged (already vectorised)."""
    W = win_s * fs
    nb = int(max(m.max() if len(m) else 0, t.max() if len(t) else 0) // W) + 1
    bm = (m // W).astype(int); bt = (t // W).astype(int)
    n_m = np.bincount(bm, minlength=nb)
    n_t = np.bincount(bt, minlength=nb)
    tot = n_m + n_t
    dis = np.bincount(bm[~pm], minlength=nb) + np.bincount(bt[~pt], minlength=nb)
    noisy = dis / np.maximum(tot, 1) > thr
    ti_fail = (tot >= 2) & (n_t < 0.7 * np.maximum(n_m, 1))
    noisy = noisy & ~ti_fail
    keep_m = m[~noisy[bm]]
    keep_f = fused[noisy[np.minimum((fused // W).astype(int), nb - 1)]]
    return np.sort(np.r_[keep_m, keep_f]).astype(np.int64), noisy.mean()


def run_chain(x, fs, snap=True, recover=True, params=None, startup="robust"):
    """Sorted, unique, 0-based native sample indices. With snap (default) MinEl and TiTen detections are moved to the
    R extremum with one polarity (decided on MinEl's first minute) before they are fused. params: FUSION overrides."""
    P = dict(FUSION, **(params or {}))
    lo, hi = P["band"]
    b, a = butter(int(PARAMS["minel"]["order"]), [lo / (fs / 2), hi / (fs / 2)], btype="bandpass")
    xf = filtfilt_matlab(b, a, x)
    val, idx = minmax_stream(xf, fs, 10)
    tsec = (idx + 1) / fs
    lm = minel_detect(val, 10, PARAMS["minel"]["params"], tsec)[0]
    lt = titen_detect(val, 10, dict(PARAMS["titen"]["params"], wd_s=RECOVER_S) if recover else PARAMS["titen"]["params"], tsec,
                      startup)[0]
    pol = 1
    if snap:
        lm, pol = snap_to_r(val, lm, tsec)
        nb = np.stack([np.maximum(lt - 1, 0), lt, np.minimum(lt + 1, len(val) - 1)], axis=1)
        lt = nb[np.arange(len(lt)), np.argmax(pol * val[nb], axis=1)] if len(lt) else lt
    m = np.unique(idx[lm])
    t = np.unique(idx[lt])
    pm, pt = pair_flags(m, t, int(round(P["agree_s"] * fs)))
    V = decode(m, t, pm, pt, fs, P, sig=pol * xf)
    return adaptive(m, t, V, pm, pt, fs, P["win_s"], P["thr"])[0]
