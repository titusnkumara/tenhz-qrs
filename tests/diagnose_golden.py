"""Locate where fusion (or any detector) first departs from the golden values.

    python tests/diagnose_golden.py                      # numba as usual
    NUMBA_DISABLE_JIT=1 python tests/diagnose_golden.py  # same code as plain Python

If the plain-Python run also differs, the cause is NumPy/SciPy/platform, not Numba.
Stages (fusion): 10 Hz stream -> MinEl beats -> TiTen beats -> agreement flags
-> decoder output -> final beats. Prints the first stage that differs per signal.
"""
import os
import platform
import sys

import numpy as np
import scipy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import tenhz  # noqa: E402
import tenhz.fusion as F  # noqa: E402
from tenhz._signal import butter, filtfilt_matlab  # noqa: E402
from tenhz.detectors import minel_detect, minmax_stream, titen_detect  # noqa: E402
from tenhz.params import PARAMS  # noqa: E402
from tests.test_golden import GOLDEN, N_SIGNALS, signal  # noqa: E402

STAGES = os.path.join(os.path.dirname(__file__), "golden_stages.npz")


def stages(x, fs):
    lo, hi = F.FUSION["band"]
    b, a = butter(int(PARAMS["minel"]["order"]), [lo / (fs / 2), hi / (fs / 2)], btype="bandpass")
    xf = filtfilt_matlab(b, a, x)
    val, idx = minmax_stream(xf, fs, 10)
    tsec = (idx + 1) / fs
    m = np.unique(idx[minel_detect(val, 10, PARAMS["minel"]["params"], tsec)[0]])
    t = np.unique(idx[titen_detect(val, 10, PARAMS["titen"]["params"], tsec)[0]])
    pm, pt = F.pair_flags(m, t, int(round(F.FUSION["agree_s"] * fs)))
    V = F.viterbi(m, t, pm, pt, fs, F.FUSION["w_single"], F.FUSION["lam"], F.FUSION["start_pen"])
    fin = F.adaptive(m, t, V, pm, pt, fs, F.FUSION["win_s"], F.FUSION["thr"])[0]
    return dict(val=val, idx=idx, m=m, t=t, pm=pm, pt=pt, V=V, final=fin)


if __name__ == "__main__":
    print(f"python {platform.python_version()} on {platform.platform()} ({platform.machine()})")
    import numba
    print(f"numpy {np.__version__}, scipy {scipy.__version__}, numba {numba.__version__}, "
          f"jit {'OFF' if os.environ.get('NUMBA_DISABLE_JIT') else 'on'}")
    if "--make" in sys.argv:                                  # (re)create the stage file
        out = {}
        for s in range(N_SIGNALS):
            for k, v in stages(*signal(s)[:2]).items():
                out[f"{k}_{s}"] = v
        np.savez_compressed(STAGES, **out)
        print("wrote", STAGES); sys.exit(0)
    G = np.load(STAGES)
    gold = np.load(GOLDEN)
    nbad = 0
    for s in range(N_SIGNALS):
        x, fs = signal(s)
        got = stages(x, fs)
        final_ok = np.array_equal(got["final"], gold[f"fusion_{s}"])
        if final_ok:
            continue
        nbad += 1
        print(f"signal {s} (fs {fs}): final beats DIFFER "
              f"({np.setdiff1d(got['final'], gold[f'fusion_{s}']).size} only here, "
              f"{np.setdiff1d(gold[f'fusion_{s}'], got['final']).size} only golden)")
        dv = (np.max(np.abs(got["val"] - G[f"val_{s}"])) if len(got["val"]) == len(G[f"val_{s}"]) else np.inf)
        print(f"    stream: same index set {np.array_equal(got['idx'], G[f'idx_{s}'])}, max |value diff| {dv:.3g}")
        for k in ("m", "t", "pm", "pt", "V"):
            same = np.array_equal(got[k], G[f"{k}_{s}"])
            extra = "" if same else f"  ({len(got[k])} here vs {len(G[f'{k}_{s}'])} golden)"
            print(f"    {k:3s}: {'same' if same else 'DIFFERS'}{extra}")
            if not same:
                break
    print(f"{nbad} of {N_SIGNALS} fusion signals differ from golden")
