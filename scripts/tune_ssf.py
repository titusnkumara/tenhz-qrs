"""Choose the threshold of NeuroKit2's ssf detector on the training databases only.

    python scripts/tune_ssf.py --data data --workers 30

Runs ssf on MITDB (118 and 119 removed) and SVDB for a grid of thresholds in (mV/s)^2, scores each
with the benchmark's own matching (+/-150 ms, one-to-one, closest first, pooled over records) and prints
pooled F1 per threshold. Put the best value into SSF_THRESHOLD in scripts/competitors.py. No test or
held-out database is read.
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import benchmark_three as bt  # noqa: E402  (records, load, match, mround, TRAIN)

GRID = [25, 50, 100, 200, 300, 400, 600, 800, 1200, 1600, 2500, 4000, 6400, 10000]


def _one(args):
    data, db, rec, fs, grid = args
    import competitors
    ecg, ann = bt.load(data, db, rec)
    tol = bt.mround(150 / 1000 * fs)
    out = []
    for thr in grid:
        os.environ["TENHZ_SSF_THRESHOLD"] = str(thr)
        try:
            det = competitors.detect("ssf", ecg, fs)
        except Exception:
            det = np.zeros(0, np.int64)
        tp, fp, fn, _ = bt.match(det, ann, tol)
        out.append((thr, tp, fp, fn))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("TENHZ_DATA") or "data")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--grid", nargs="+", type=float, default=GRID, help="thresholds in (mV/s)^2")
    a = ap.parse_args()
    recs = bt.records(a.data, bt.TRAIN)
    if not recs:
        sys.exit("no training records found")
    jobs = [(a.data, db, r, fs, list(a.grid)) for db, r, fs in recs]
    if a.workers > 1:
        with ProcessPoolExecutor(a.workers) as ex:
            res = list(ex.map(_one, jobs))
    else:
        res = [_one(j) for j in jobs]
    tot = {t: np.zeros(3, np.int64) for t in a.grid}
    for rows in res:
        for thr, tp, fp, fn in rows:
            tot[thr] += (tp, fp, fn)
    print(f"ssf on {len(recs)} training records (MITDB + SVDB), +/-150 ms, pooled")
    print(f"{'threshold (mV/s)^2':>20} {'slope (mV/s)':>13} {'Se':>7} {'PPV':>7} {'F1':>7}")
    best = None
    for t in a.grid:
        tp, fp, fn = tot[t]
        se, ppv = 100 * tp / max(tp + fn, 1), 100 * tp / max(tp + fp, 1)
        f1 = 100 * 2 * tp / max(2 * tp + fp + fn, 1)
        print(f"{t:>20g} {np.sqrt(t):>13.1f} {se:7.2f} {ppv:7.2f} {f1:7.2f}")
        if best is None or f1 > best[1]:
            best = (t, f1)
    print(f"\nbest: SSF_THRESHOLD = {best[0]:g}  (F1 {best[1]:.2f}); set it in scripts/competitors.py")


if __name__ == "__main__":
    main()
