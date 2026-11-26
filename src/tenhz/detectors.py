"""MinEl and TiTen: bandpass at the native rate, 10 Hz min-max stream, detection.
"""
import math

import numpy as np

from ._kernels import (_k_amp_ref, _k_minel_main, _k_minel_searchback, _k_minmax_stream,
                       _k_peaks_dist, _k_peaks_time, _k_region_max, _k_titen_main)
from ._signal import EPS, block_edges, butter, conv_same, filtfilt_matlab, mround
from .params import PARAMS


def minmax_stream(xf, fs, ft=10):
    """Min and max of every 2/ft-second block, in sample order: values and
    0-based native sample indices."""
    edges = block_edges(len(xf), 2 * fs / ft)
    return _k_minmax_stream(np.ascontiguousarray(xf, dtype=np.float64), edges)


def minel_detect(x, fs, p, tsec):
    """minel_detect with compiled loops."""
    x = np.asarray(x, float)
    n = len(x)
    if n < p["win_beat"] + 10:
        return np.zeros(0, np.int64), np.zeros(0)
    xa = np.abs(x)
    ma_qrs = conv_same(xa, np.ones(int(p["win_qrs"])) / p["win_qrs"])
    ma_beat = conv_same(xa, np.ones(int(p["win_beat"])) / p["win_beat"])
    hybrid = p["alpha"] * xa + (1 - p["alpha"]) * ma_qrs
    thr1 = ma_beat + p["beta"] * np.mean(xa)
    d = np.diff(np.r_[0, (hybrid > thr1).astype(np.int8), 0])
    st = np.flatnonzero(d == 1)
    en = np.flatnonzero(d == -1) - 1
    if len(st) == 0:
        return np.zeros(0, np.int64), np.zeros(0)
    xs = xa if p["use_abs"] else x
    ca, cl = _k_region_max(xs, st, en)
    if p.get("thr_mode", "local") != "local":
        raise NotImplementedError("only thr_mode 'local' is ported")
    half = max(1, int(mround(p["win_len_sec"] * fs / 2)))
    base = p["height_thresh"] * _k_amp_ref(xs, cl, half, float(p["perc"]))
    okf = ca >= p.get("amp_floor", 0)
    cl, ca, base = cl[okf], ca[okf], base[okf]
    if p.get("time_rules", 0) and tsec is not None:
        tau = np.asarray(tsec, float); refr = p["refr_s"]
    else:
        tau = np.arange(1, n + 1, dtype=float); refr = p["min_dist"]
    cl_all, ca_all, base_all = cl, ca, base
    keep = ca > base
    locs, amps = _k_minel_main(cl[keep], ca[keep], base[keep], tau, float(refr),
                               float(p["rr_tol"]), float(p["low_factor"]), float(p["high_factor"]))
    if p["use_sb"] and len(locs) >= 3:
        rel = ca_all / np.maximum(base_all, EPS)
        locs, amps = _k_minel_searchback(locs, amps, tau, cl_all, ca_all, base_all, rel, tau[cl_all],
                                         float(refr), float(p["sb_mult"]), float(p["sb_frac"]))
    return locs.astype(np.int64), amps.astype(float)


def titen_detect(ecg, fs, p, tsec):
    """titen_detect with compiled loops."""
    ecg = np.asarray(ecg, float)
    n = len(ecg)
    empty = (np.zeros(0, np.int64), np.zeros(0))
    if n < 25:
        return empty
    w = int(p["win_len"])
    integ = np.convolve(ecg ** 2, np.ones(w) / w)
    ni = len(integ)
    use_time = bool(p.get("time_rules", 0)) and tsec is not None
    if ni < 3:
        return empty
    ismax = np.r_[False, (integ[1:-1] > integ[:-2]) & (integ[1:-1] > integ[2:]), False]
    pl = np.flatnonzero(ismax)
    pk = integ[pl]
    if len(pl) == 0:
        return empty
    order = np.argsort(-pk, kind="stable")
    if use_time:
        tau = np.asarray(tsec, float)[np.minimum(np.arange(ni), n - 1)]
        refr = p["refr_s"]
        keep = _k_peaks_time(ni, pl, order, tau, float(refr), int(math.ceil(refr * 20) + 2))
    else:
        tau = np.arange(1, ni + 1, dtype=float)
        keep = _k_peaks_dist(ni, pl, order, int(p["min_dist"]))
    pk, pl = pk[keep], pl[keep]
    if len(pk) == 0:
        return empty
    init_len = min(20, n)
    thr_sig = float(np.max(integ[:init_len])) * p["init_sig_frac"]
    thr_noise = float(np.mean(integ[:init_len])) * 0.5
    thr_sig_r = float(np.max(ecg[:init_len])) * p["init_sig_frac"]
    thr_noise_r = float(np.mean(ecg[:init_len])) * 0.5
    locs, amps = _k_titen_main(pk, pl.astype(np.int64), integ, ecg, tau, use_time,
                               float(p["sig_alpha"]), float(p["rr_start_beats"]),
                               float(p["rr_low_frac"]), float(p["rr_high_frac"]),
                               float(p["searchback_mult"]), float(p["thr_frac"]),
                               thr_sig, thr_noise, thr_sig_r, thr_noise_r)
    if p.get("amp_floor", 0) > 0:
        k = amps >= p["amp_floor"]
        locs, amps = locs[k], amps[k]
    return locs, amps


def snap_to_r(val, loc, tsec, minute_s=60.0, max_minutes=10, ratio=1.2, min_beats=20):
    """Move each detection to the R extremum of its QRS, with one QRS polarity per record.

    val, tsec : the 10 Hz min-max stream (values) and its times in s
    loc       : detections as indices into the stream
    In the stream every block contributes its minimum and maximum in sample order, so the R and S extrema of a QRS are
    neighbours: each detection is compared with the stream samples just before and after it (three gathers per beat).
    Polarity comes from the first decisive minute (scanning at most `max_minutes`): a minute with at least `min_beats`
    detections in which the median peak height and the median trough depth differ by at least `ratio`. This skips
    start-up artefact; if no minute is decisive (biphasic QRS) the record is treated as upright.
    Returns (snapped stream indices, polarity +1/-1).
    """
    loc = np.asarray(loc, np.int64)
    if len(loc) == 0:
        return loc, 1
    n = len(val)
    nb = np.stack([np.maximum(loc - 1, 0), loc, np.minimum(loc + 1, n - 1)], axis=1)
    v = np.asarray(val, float)[nb]
    t = np.asarray(tsec, float)[loc]
    edges = np.searchsorted(t, t[0] + minute_s * np.arange(max_minutes + 1), side="left")
    pol = 1
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < min_beats:
            continue
        pk = np.median(v[a:b].max(axis=1)); tr = np.median(-v[a:b].min(axis=1))
        if max(pk, tr) >= ratio * min(pk, tr):
            pol = 1 if pk >= tr else -1
            break
    j = np.argmax(pol * v, axis=1)
    return nb[np.arange(len(loc)), j], pol


def run_chain(x, fs, detector, snap=True):
    """Native-rate bandpass -> 10 Hz min-max stream -> detector -> (snap to R, default).
    Returns sorted, unique, 0-based native sample indices. snap=False gives the published detections."""
    det = PARAMS[detector]
    lo, hi = det["band"]
    b, a = butter(int(det["order"]), [lo / (fs / 2), hi / (fs / 2)], btype="bandpass")
    xf = filtfilt_matlab(b, a, x)
    val, idx = minmax_stream(xf, fs, 10)
    tsec = (idx + 1) / fs                       # MATLAB's 1-based times
    fn = minel_detect if detector == "minel" else titen_detect
    loc, _ = fn(val, 10, det["params"], tsec)
    if snap:
        loc, _ = snap_to_r(val, loc, tsec)
    return np.unique(idx[loc])
