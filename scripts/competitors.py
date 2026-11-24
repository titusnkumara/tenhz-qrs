"""Competing QRS detectors behind one calling convention: detect(name, ecg, fs) -> 0-based
sample indices (sorted, unique, int64). All see the same single lead at its native rate.

Groups (names for --methods in benchmark_three.py):
  ours      minel titen fusion                       (Numba, this package)
  nk2       the 17 NeuroKit2 R-peak methods except promac (ssf with a threshold set on the
            training databases, see SSF_THRESHOLD)
  wfdb      xqrs gqrs                                (WFDB-Python, the PhysioNet reference code)
  pyecg     pyecg_matched ...                        (py-ecg-detectors; only the detector that
                                                      NeuroKit2 does not already provide)
  baseline  energy                                   (band-pass, squared, moving average, find_peaks)
  all       everything above

NeuroKit2 methods run on ECG cleaned with their own ecg_clean method where NeuroKit2 has one,
otherwise 'neurokit', exactly as tools/benchmark_neurokit.py did.
"""
from __future__ import annotations

import numpy as np

NK2 = ["neurokit", "pantompkins1985", "hamilton2002", "christov2004", "engzeemod2012",
       "elgendi2010", "elgendi2013", "gamboa2008", "zong2003", "ssf", "manikandan2012",
       "kalidas2017", "khamis2016", "nabian2018", "martinez2004", "rodrigues2021", "emrich2023"]
NK_CLEAN = {"pantompkins1985": "pantompkins1985", "hamilton2002": "hamilton2002",
            "engzeemod2012": "engzeemod2012", "elgendi2010": "elgendi2010", "gamboa2008": "biosppy",
            "emrich2023": "vg"}
WFDB = ["xqrs", "gqrs"]
PYECG = ["pyecg_matched"]
BASELINE = ["energy"]
OURS = ["minel", "titen", "fusion"]
GROUPS = {"ours": OURS, "nk2": NK2, "wfdb": WFDB, "pyecg": PYECG, "baseline": BASELINE,
          "all": OURS + NK2 + WFDB + PYECG + BASELINE}

LABEL = {"minel": "MinEl (10 Hz)", "titen": "TiTen (10 Hz)", "fusion": "Fusion (10 Hz, exploratory)",
         "xqrs": "WFDB xqrs", "gqrs": "WFDB gqrs", "pyecg_matched": "Matched filter (py-ecg-detectors)",
         "energy": "Energy + find_peaks (baseline)"}


def expand(names):
    """Method names and group names -> list of method names, order kept, no duplicates."""
    out = []
    for n in names:
        for m in GROUPS.get(n, [n]):
            if m not in out:
                out.append(m)
    bad = [m for m in out if m not in GROUPS["all"]]
    if bad:
        raise SystemExit(f"unknown method(s) {bad}; choose from {sorted(GROUPS)} or {GROUPS['all']}")
    return out


def _idx(x):
    return np.unique(np.asarray(x, dtype=np.int64))


# A note about the last call (e.g. a fallback that was used), read once by take_note().
_NOTE = [""]


def take_note():
    n, _NOTE[0] = _NOTE[0], ""
    return n


def _manikandan_fixed(signal, sampling_rate=1000, **kwargs):
    """NeuroKit2 0.3.0 _ecg_findpeaks_manikandan with its final search-window loop removed. That
    loop builds `id_maxes`, which is never returned (the function returns the zero-crossing indices
    `idx`), but it raises IndexError / TypeError when a zero crossing lies within half a window of
    the end of the record (and, with NumPy >= 2.4, on a size-1 array argument to np.arange).
    Everything that determines the output is unchanged."""
    import scipy.signal

    def cheby1_bandpass_filter(data, lowcut, highcut, fs, order=5, rp=1):
        nyq = 0.5 * fs
        b, a = scipy.signal.cheby1(order, rp=rp, Wn=[lowcut / nyq, highcut / nyq], btype="bandpass")
        return scipy.signal.filtfilt(b, a, data)

    def running_mean(x, N):
        cumsum = np.cumsum(np.insert(x, 0, 0))
        return (cumsum[N:] - cumsum[:-N]) / float(N)

    filtered = cheby1_bandpass_filter(signal, lowcut=6, highcut=18, fs=sampling_rate, order=4)
    dn = np.append(filtered[1:], 0) - filtered
    if np.max(abs(dn)) == 0:
        return np.array([])
    dtn = dn / (np.max(abs(dn)))
    sn = -(dtn ** 2) * np.log10(dtn ** 2)
    window_len = int(0.15 * sampling_rate)
    window_len = window_len - 1 if window_len % 2 == 0 else window_len
    window = scipy.signal.windows.boxcar(window_len)
    see = scipy.signal.convolve(sn, window, mode="same")
    see = np.flip(see)
    see = scipy.signal.convolve(see, window, "same")
    see = np.flip(see)
    ht = np.imag(scipy.signal.hilbert(see))
    ma_len = int(2.5 * sampling_rate)
    ma_out = np.insert(running_mean(ht, ma_len), 0, [0] * (ma_len - 1))
    zn = ht - ma_out
    idx = np.argwhere(np.diff(np.sign(zn)) > 0).flatten().tolist()
    return np.array(idx)


# ------------------------------------------------------------------ NeuroKit2
def _vg_fixed(signal, sampling_rate=1000, window_seconds=2, window_overlap=0.5, accelerated=True, **kwargs):
    """NeuroKit2 0.3.0 _ecg_findpeaks_visibilitygraph (emrich2023) with one guard added. When a 2 s window
    has no local maximum above its median (a flat or saturated stretch of signal), the stock code builds a
    visibility graph from an empty segment and raises IndexError. Here such a window gets zero weights and
    the loop moves on. Every other line, and the result for every window that has candidate peaks, is the
    same as the installed code."""
    import scipy.signal
    import ts2vg
    from neurokit2.ecg.ecg_findpeaks import _ecg_findpeaks_visgraphthreshold

    N = len(signal)
    M = int(window_seconds * sampling_rate)
    L = 0
    R = M
    dM = int(np.ceil(window_overlap * M))
    n_segments = int(np.ceil((N - R) / (M - dM) + 1))
    weights = np.zeros(N)
    BETA = 0.55
    if np.max(signal) == np.min(signal):
        return np.array([])
    if N < M:
        M, R = N, N
        n_segments = 1
    for jj in range(n_segments):
        segment = signal[L:R]
        if accelerated:
            threshold = np.quantile(segment, 0.5)
            indices_maxima = scipy.signal.find_peaks(segment)[0]
            segment_maxima = segment[indices_maxima]
            greater = np.argwhere(segment_maxima > threshold).flatten()
            indices = indices_maxima[greater]
            segment = segment_maxima[greater]
        else:
            indices = np.arange(len(segment))
        if len(indices) >= 2:                                   # the added guard
            A = ts2vg.NaturalVG(directed="top_to_bottom").build(segment, indices).adjacency_matrix()
            size = len(indices)
            w = np.ones(size) / size
            while np.count_nonzero(w) / size >= BETA:
                _w = np.matmul(A, w)
                _w /= np.linalg.norm(_w)
                if np.any(np.isnan(_w)):
                    break
                w = _w
            if L == 0:
                weights[L + indices] = w
            elif N - dM + 1 <= L and L + 1 <= N:
                weights[L + indices] = 0.5 * (weights[L + indices] + w)
            else:
                weights[L + indices[indices <= dM]] = 0.5 * (w[indices <= dM] + weights[L + indices[indices <= dM]])
                weights[L + indices[indices > dM]] = w[indices > dM]
        if R - L < M:
            break
        L += M - dM
        if R + (M - dM) <= N:
            R += M - dM
        else:
            R = N
    rpeaks = _ecg_findpeaks_visgraphthreshold(signal * weights, sampling_rate)
    return np.array(rpeaks, dtype="int")


# NeuroKit2's ssf (slope sum function, from BioSPPy) marks a beat where the squared negative first
# difference exceeds a fixed threshold, default 20 in raw signal units per sample. For ECG in mV that is
# never reached, so with the default it finds no beats at all. The threshold below is in (mV/s)^2, so one
# value holds at every sampling rate (it is divided by fs^2 before the call). It was chosen on the training
# databases only (MITDB without 118/119, and SVDB) with scripts/tune_ssf.py, as MinEl and TiTen were.
# TENHZ_SSF_THRESHOLD overrides it (tune_ssf.py uses this).
SSF_THRESHOLD = 400.0


def _ssf_threshold():
    import os
    return float(os.environ.get("TENHZ_SSF_THRESHOLD", SSF_THRESHOLD))


def _nk2(method, ecg, fs):
    import neurokit2 as nk
    clean = nk.ecg_clean(ecg, sampling_rate=fs, method=NK_CLEAN.get(method, "neurokit"))
    kw = {"threshold": _ssf_threshold() / fs ** 2} if method == "ssf" else {}
    try:
        _, info = nk.ecg_peaks(clean, sampling_rate=fs, method=method, correct_artifacts=False, **kw)
    except Exception as e:
        # A few NeuroKit2 detectors raise on some records because of a bug in the library, not because the
        # method fails. Re-run only those records with a patched copy that fixes the bug; every other
        # record runs the installed NeuroKit2 code unchanged.
        patched = {"manikandan2012": ("_ecg_findpeaks_manikandan", _manikandan_fixed,
                                      "patched copy without its unused index block"),
                   "emrich2023": ("_ecg_findpeaks_visibilitygraph", _vg_fixed,
                                  "patched copy that skips windows with no candidate peaks")}
        if method not in patched:
            raise
        import importlib
        attr, fixed, why = patched[method]
        mod = importlib.import_module("neurokit2.ecg.ecg_findpeaks")
        orig = getattr(mod, attr)
        setattr(mod, attr, fixed)
        try:
            _, info = nk.ecg_peaks(clean, sampling_rate=fs, method=method, correct_artifacts=False, **kw)
        finally:
            setattr(mod, attr, orig)
        _NOTE[0] = f"NeuroKit2 raised {type(e).__name__}; {why}"
    return _idx(info["ECG_R_Peaks"])


# ----------------------------------------------------------------------- WFDB
def _wfdb(method, ecg, fs):
    from wfdb import processing
    if method == "xqrs":
        return _idx(processing.xqrs_detect(sig=ecg, fs=fs, verbose=False))
    return _idx(processing.gqrs_detect(sig=ecg, fs=fs))


# ------------------------------------------------------------ py-ecg-detectors
_TEMPLATE_FS = 250


def _pyecg_matched(ecg, fs):
    """py-ecg-detectors matched filter. Its stock templates exist only for 250 and 360 Hz, so
    for other rates the 250 Hz template is resampled linearly to fs (an adaptation, stated in
    the paper)."""
    from ecgdetectors import Detectors
    import ecgtemplates
    det = Detectors(fs)
    if fs in (250, 360):
        return _idx(det.matched_filter_detector(ecg))
    t0 = np.asarray(ecgtemplates.qrs_250Hz, float)
    n = max(3, int(round(len(t0) * fs / _TEMPLATE_FS)))
    t = np.interp(np.linspace(0, len(t0) - 1, n), np.arange(len(t0)), t0)
    import tempfile
    import os
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        np.savetxt(f, t)
        name = f.name
    try:
        return _idx(det.matched_filter_detector(ecg, template_file=name))
    finally:
        os.remove(name)


# ------------------------------------------------------------------- baseline
def _energy(ecg, fs):
    """Simple, fully specified baseline: 5-20 Hz zero-phase band-pass, squared, 100 ms moving
    average, local maxima at least 250 ms apart that exceed 0.35 x the running 8 s maximum of
    the energy; each peak moved to the largest |band-passed| sample within +/-75 ms."""
    from scipy.ndimage import maximum_filter1d, uniform_filter1d
    from scipy.signal import butter, filtfilt, find_peaks
    b, a = butter(2, [5 / (fs / 2), 20 / (fs / 2)], btype="bandpass")
    xf = filtfilt(b, a, ecg)
    e = uniform_filter1d(xf ** 2, max(1, int(round(0.1 * fs))))
    thr = 0.35 * maximum_filter1d(e, max(3, int(round(8 * fs))))
    pk, _ = find_peaks(e, height=thr, distance=max(1, int(round(0.25 * fs))))
    w = int(round(0.075 * fs))
    out = []
    for p in pk:
        lo, hi = max(0, p - w), min(len(xf), p + w + 1)
        out.append(lo + int(np.argmax(np.abs(xf[lo:hi]))))
    return _idx(out)


# ------------------------------------------------------------------- dispatch
def detect(name, ecg, fs):
    """0-based R-peak sample indices of detector `name`."""
    ecg = np.asarray(ecg, float)
    if name in NK2:
        return _nk2(name, ecg, fs)
    if name in WFDB:
        return _wfdb(name, ecg, fs)
    if name == "pyecg_matched":
        return _pyecg_matched(ecg, fs)
    if name == "energy":
        return _energy(ecg, fs)
    raise ValueError(f"unknown detector {name!r}")