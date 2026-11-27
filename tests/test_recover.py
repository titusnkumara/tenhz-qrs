"""recover (default True) and startup (default "robust"): TiTen must not stay silent after a large transient or a
rail-saturated stretch at the start of a record.

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
    old = dict(startup="2.2")                                                              # 2.2 start-up
    assert _hit_rate(tenhz.qrs_locations(y, fs, "titen", recover=False, **old), beats, fs) < 0.2  # no watchdog: stuck
    assert _hit_rate(tenhz.qrs_locations(y, fs, "titen", **old), beats, fs) > 0.9                 # 2.2 watchdog recovers
    assert _hit_rate(tenhz.qrs_locations(y, fs, "titen"), beats, fs) > 0.9                        # default (2.3) recovers
    assert _hit_rate(tenhz.qrs_locations(y, fs, "fusion"), beats, fs) > 0.9


def test_startup_rail():
    """2.3 start-up (default): a short record that starts with 8 s of rail saturation (5.5 mV, sign flips) is not lost;
    the 2.2 start-up loses it. Recovery takes a few seconds (5 s watchdog), so beats from 4 s after the clip are checked."""
    fs = 250
    x, beats, t = _ecg(fs, minutes=1, seed=2)
    x = x[:30 * fs]; beats = beats[beats < 29.5]
    y = x.copy(); e = 8 * fs
    y[:e] = 5.5 * np.where((np.arange(e) // (fs // 2)) % 2 == 0, 1.0, -1.0)
    y[e:] += y[e - 1] * np.exp(-np.arange(len(y) - e) / (0.5 * fs))
    after = beats[beats > 12.0]
    hit = lambda d: _hit_rate(d, after, fs)
    assert hit(tenhz.qrs_locations(y, fs, "titen")) > 0.9
    assert hit(tenhz.qrs_locations(y, fs, "titen", startup="2.2")) < 0.5
    assert hit(tenhz.qrs_locations(y, fs, "fusion")) > 0.9


if __name__ == "__main__":
    test_recover()
    test_startup_rail()
    print("ok")
