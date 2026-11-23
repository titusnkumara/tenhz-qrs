"""Vectorised helpers, copied unchanged from python/src/tenhz_qrs/core.py.

They follow MATLAB exactly (round half away from zero, conv(..., 'same')
alignment, filtfilt edge handling), so detections match the published results.
"""
import math

import numpy as np
from scipy.signal import butter, lfilter  # noqa: F401  (butter is re-exported)

EPS = np.finfo(float).eps


def mround(x):
    """MATLAB round: halves away from zero."""
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def conv_same(u, w):
    """MATLAB conv(u, w, 'same') for vectors (also for even-length w)."""
    full = np.convolve(u, w)
    start = len(w) // 2
    return full[start:start + len(u)]


def filtfilt_matlab(b, a, x):
    """MATLAB filtfilt for one channel: odd reflection of 3*(nfilt-1) samples
    at each end and MATLAB's steady-state initial conditions."""
    b = np.asarray(b, float); a = np.asarray(a, float)
    nfilt = max(len(a), len(b))
    if len(b) < nfilt: b = np.r_[b, np.zeros(nfilt - len(b))]
    if len(a) < nfilt: a = np.r_[a, np.zeros(nfilt - len(a))]
    b = b / a[0]; a = a / a[0]
    nfact = max(1, 3 * (nfilt - 1))
    if len(x) <= nfact:
        raise ValueError("signal too short for filtfilt")
    m = nfilt - 1
    A = np.eye(m)
    A[:, 0] += a[1:]
    A[:-1, 1:] -= np.eye(m - 1) if m > 1 else 0
    zi = np.linalg.solve(A, b[1:] - b[0] * a[1:])
    xt = 2 * x[0] - x[nfact:0:-1]
    _, zo = lfilter(b, a, xt, zi=zi * xt[0])
    yc2, zo = lfilter(b, a, x, zi=zo)
    xt = 2 * x[-1] - x[-2:-nfact - 2:-1]
    yc3, _ = lfilter(b, a, xt, zi=zo)
    _, zo = lfilter(b, a, yc3[::-1], zi=zi * yc3[-1])
    yc5, _ = lfilter(b, a, yc2[::-1], zi=zo)
    return yc5[::-1]


def block_edges(n, B):
    """Edges of consecutive blocks of B (possibly non-integer) samples."""
    nb = int(math.floor(n / B + 1e-9))
    return mround(np.arange(nb + 1) * B).astype(np.int64)
