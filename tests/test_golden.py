"""Regression test: detections on synthetic ECG must equal the reference.

tests/golden.npz holds detections from the pure-Python reference
(python/src/tenhz_qrs/core.py, itself checked beat for beat against MATLAB;
fusion: tools/tenhz_fast.detect_reference).

    python tests/test_golden.py          # or: pytest tests
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import tenhz  # noqa: E402

GOLDEN = os.path.join(os.path.dirname(__file__), "golden.npz")
FS = [128, 250, 360, 500]
NOISE = [0.01, 0.05, 0.2, 0.6]
N_SIGNALS = 16


def _synth(fs, seconds, seed, noise):
    """Synthetic ECG: irregular rhythm, ectopics, pauses, inverted beats, wander,
    noise bursts, flat lines and drift."""
    rng = np.random.default_rng(seed)
    n = int(seconds * fs); t = np.arange(n) / fs
    x = rng.uniform(0, 0.5) * np.sin(2 * np.pi * rng.uniform(0.1, 0.4) * t + rng.uniform(0, 6))
    hr = rng.uniform(45, 150); tb = rng.uniform(0.2, 1.0)
    while tb < seconds - 0.5:
        amp = rng.uniform(0.4, 1.8) * (-1 if rng.random() < 0.05 else 1)
        for c, wd, a in [(-0.2, 0.025, 0.12), (0, rng.uniform(0.008, 0.03), 1.0), (0.035, 0.012, -0.25),
                         (0.28, 0.05, 0.3)]:
            i0 = max(0, int((tb + c - 5 * wd) * fs)); i1 = min(n, int((tb + c + 5 * wd) * fs) + 1)
            x[i0:i1] += amp * a * np.exp(-(t[i0:i1] - tb - c) ** 2 / (2 * wd ** 2))
        rr = 60 / hr * (1 + 0.1 * rng.standard_normal())
        if rng.random() < 0.05: rr *= 0.6
        if rng.random() < 0.01: rr *= 2.2
        tb += max(0.25, rr)
    x += noise * rng.standard_normal(n)
    for _ in range(int(rng.integers(1, 6))):
        s = int(rng.integers(0, max(1, n - 5 * fs))); L = min(int(rng.uniform(0.3, 6) * fs), n - s)
        kind = rng.integers(0, 4)
        if kind == 0: x[s:s + L] += rng.standard_normal(L) * rng.uniform(0.3, 3)
        elif kind == 1: x[s:s + L] = x[s]
        elif kind == 2: x[s:s + L] += np.cumsum(rng.standard_normal(L)) * 0.05
        else: x[s:s + L] = np.round(x[s:s + L] * 4) / 4                    # coarse ADC: ties in extrema
    return x.astype(np.float32).astype(float)


def signal(s):
    fs = FS[s % 4]
    return _synth(fs, 300 * (1 + s % 3), seed=s, noise=NOISE[s // 4 % 4]), fs


def test_golden():
    g = np.load(GOLDEN)
    bad = []
    for s in range(N_SIGNALS):
        x, fs = signal(s)
        for d in tenhz.ALL_DETECTORS:
            if not np.array_equal(tenhz.qrs_locations(x, fs, d), g[f"{d}_{s}"]):
                bad.append(f"{d} signal {s} (fs {fs})")
    assert not bad, "differ from reference: " + ", ".join(bad)


if __name__ == "__main__":
    test_golden()
    print(f"{N_SIGNALS} synthetic signals x {len(tenhz.ALL_DETECTORS)} detectors: identical to the reference")
