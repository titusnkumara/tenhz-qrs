"""tenhz: MinEl and TiTen QRS detection from a 10 Hz min-max representation
of the ECG (Numba-compiled; detections identical to the published reference).

    >>> import tenhz
    >>> r = tenhz.detect(ecg, fs)                       # MinEl (default)
    >>> r.locations                                     # 0-based R-peak sample indices
    >>> r.peaks                                         # ECG values there
    >>> locs = tenhz.qrs_locations(ecg, fs, "titen")    # indices only

Offline processing (zero-phase bandpass). The first call in a fresh environment
compiles the kernels (a few seconds); later runs load them from Numba's cache.
Call tenhz.warmup() first if you are timing.
"""
from collections import namedtuple

import numpy as np

from . import fusion as _fusion
from .detectors import run_chain
from .params import PARAMS

__all__ = ["detect", "qrs_locations", "QRS", "parameters", "warmup", "DETECTORS", "EXPLORATORY", "ALL_DETECTORS"]
__version__ = "2.0.0"

DETECTORS = ("minel", "titen")          # the published detectors
EXPLORATORY = ("fusion",)                # opt-in: detect(..., "fusion")
ALL_DETECTORS = DETECTORS + EXPLORATORY

QRS = namedtuple("QRS", ["locations", "peaks"])
QRS.__doc__ = "R peaks: locations (0-based sample indices, int64) and peaks (ECG values there)."


def parameters():
    """Published parameters (band, filter order, detector settings), as a copy."""
    import copy
    return copy.deepcopy(PARAMS)


def detect(ecg, fs=None, detector="minel", *, sampling_rate=None):
    """Detect R peaks in a single-lead ECG.

    ecg      : 1-D array-like, one lead (mV) at its native sampling rate
    fs       : sampling rate in Hz (``sampling_rate=`` is accepted as an alias)
    detector : "minel" (default), "titen", or "fusion" (exploratory, not in the paper's
               published results: MinEl+TiTen on one shared stream with an RR decoder)

    Returns QRS(locations, peaks).
    """
    fs = fs if fs is not None else sampling_rate
    if fs is None:
        raise ValueError("give the sampling rate as fs or sampling_rate")
    fs = float(fs)
    x = np.asarray(ecg, dtype=float).squeeze()
    if x.ndim != 1:
        raise ValueError("ecg must be one lead (1-D)")
    if len(x) < 10 * fs:
        raise ValueError("ecg must be at least 10 s long")
    if not np.all(np.isfinite(x)):
        raise ValueError("ecg contains NaN or Inf (replace them, e.g. with 0, first)")
    key = str(detector).lower()
    if key not in ALL_DETECTORS:
        raise ValueError(f"detector must be one of {ALL_DETECTORS}, not {detector!r}")
    band = _fusion.FUSION["band"] if key == "fusion" else PARAMS[key]["band"]
    if fs <= 2 * band[1]:
        raise ValueError(f"fs = {fs:g} Hz is too low for the {band[0]}-{band[1]} Hz band")
    loc = (_fusion.run_chain(x, fs) if key == "fusion" else run_chain(x, fs, key)).astype(np.int64)
    return QRS(loc, x[loc])


def qrs_locations(ecg, fs, detector="minel"):
    """Only the 0-based R-peak sample indices."""
    return detect(ecg, fs, detector).locations


def warmup():
    """Compile (or load from Numba's cache) every kernel on a short synthetic signal."""
    fs = 250
    t = np.arange(30 * fs) / fs
    x = np.sin(2 * np.pi * 1.2 * t) ** 41 + 0.02 * np.random.default_rng(0).standard_normal(len(t))
    for d in ALL_DETECTORS:
        detect(x, fs, d)
