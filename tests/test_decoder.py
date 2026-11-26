"""2.2 Fusion decoder: the 2.1 decoder stays reproducible, and the new default still finds the beats.

    python tests/test_decoder.py
"""
import hashlib

import numpy as np

import tenhz
from tenhz import fusion
from test_recover import _ecg, _hit_rate

GOLDEN_2_1 = "df03eb514f69234805f0fbbe559d8907082d2c3b853bacc32e24f308355bb02a"   # tenhz 2.1.1 fusion, same signal


def _signal():
    x, beats, t = _ecg(fs=250, minutes=10, seed=3)
    rng = np.random.default_rng(7)
    return x + 0.15 * np.sin(2 * np.pi * 0.3 * t) + 0.08 * rng.standard_normal(len(x)), beats


def test_fusion_2_1_reproduced():
    x, _ = _signal()
    loc = fusion.run_chain(x, 250, params=fusion.FUSION_2_1).astype(np.int64)
    assert len(loc) == 704
    assert hashlib.sha256(loc.tobytes()).hexdigest() == GOLDEN_2_1


def test_fusion_2_2_default():
    x, beats = _signal()
    loc = tenhz.detect(x, 250, "fusion").locations
    assert _hit_rate(loc, beats, 250) > 0.99
    assert len(loc) < 1.02 * len(beats)


if __name__ == "__main__":
    test_fusion_2_1_reproduced(); test_fusion_2_2_default(); print("ok")
