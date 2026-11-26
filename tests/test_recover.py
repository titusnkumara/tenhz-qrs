"""recover (default True): TiTen must not stay silent after a large transient at the start of a record.

    python tests/test_recover.py
"""
import numpy as np

import tenhz


def _ecg(fs=250, minutes=10, seed=1):
    rng = np.random.default_rng(seed)
    n = minutes * 60 * fs
    t = np.arange(n) / fs
    beats = np.arange(1.0, minutes * 60 - 1, 0.85)
    beats = beats + rng.normal(0, 0.01, len(beats))
    x = np.zeros(n)
    for b in beats:
        x += np.exp(-0.5 * ((t - b) / 0.012) ** 2) - 0.25 * np.exp(-0.5 * ((t - b - 0.03) / 0.012) ** 2)
    return x + 0.02 * rng.standard_normal(n), beats, t


def _hit_rate(found, beats, fs):
    found = np.sort(found)
    b = (beats * fs).astype(int)
    j = np.clip(np.searchsorted(found, b), 1, len(found) - 1)
    return float(np.mean(np.minimum(np.abs(found[j] - b), np.abs(found[j - 1] - b)) <= 0.15 * fs))


def test_recover():
    fs = 250
    x, beats, t = _ecg(fs)
    y = x.copy()
    s = slice(int(0.8 * fs), int(2.5 * fs))
    y[s] += 25 * np.sin(2 * np.pi * 12 * t[s])                       # 25 mV burst at the start
    clean = tenhz.qrs_locations(x, fs, "titen")
    assert _hit_rate(clean, beats, fs) > 0.99
    assert np.array_equal(clean, tenhz.qrs_locations(x, fs, "titen", recover=False))  # no change without a failure
    assert _hit_rate(tenhz.qrs_locations(y, fs, "titen", recover=False), beats, fs) < 0.2  # recover=False stays stuck
    assert _hit_rate(tenhz.qrs_locations(y, fs, "titen"), beats, fs) > 0.9               # default recovers
    assert _hit_rate(tenhz.qrs_locations(y, fs, "fusion"), beats, fs) > 0.9


if __name__ == "__main__":
    test_recover()
    print("ok")
