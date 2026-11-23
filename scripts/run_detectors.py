#!/usr/bin/env python3
"""Run MinEl and TiTen (and with --fusion the exploratory fusion detector) on exported records and write the detections.

    python scripts/run_detectors.py --data data --dbs nsrdb --out results/run_detectors_v1

Writes <out>/<db>/<rec>-<detector>.txt (0-based sample indices, one per line)
and <out>/summary.csv (record, detector, beats, reference beats, seconds).
Speed here includes file reading; quote speed only from benchmark_timing.py.
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import tenhz  # noqa: E402


def col(f):
    return pd.read_csv(f, header=None, sep=r"\s+", usecols=[0], engine="c").to_numpy(dtype=float, copy=True).ravel()


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default="data")
ap.add_argument("--dbs", nargs="+", required=True)
ap.add_argument("--records", nargs="*", default=None, help="only these records")
ap.add_argument("--fusion", action="store_true", help="also run the exploratory fusion detector")
ap.add_argument("--out", required=True, help="new folder for the results")
a = ap.parse_args()
if os.path.exists(a.out) and os.listdir(a.out):
    sys.exit(f"{a.out} exists and is not empty: give a new --out folder")
tenhz.warmup()
rows = []
for db in a.dbs:
    d = os.path.join(a.data, f"{db}-mat")
    names = sorted(f[:-4] for f in os.listdir(d) if f.endswith(".txt") and not f.endswith(("-ann.txt", "-fs.txt")))
    if a.records:
        names = [n for n in names if n in a.records]
    os.makedirs(os.path.join(a.out, db), exist_ok=True)
    FS = {db: {k: v["fs"] for k, v in json.load(open(os.path.join(d, "_manifest.json"))).items()
               if v.get("status") == "OK"}}      # sampling rates written by data/download_physionet.py
    names = [n for n in names if n in FS[db]]
    for rec in names:
        base = os.path.join(d, rec)
        x = col(base + ".txt"); x[~np.isfinite(x)] = 0
        fs = FS[db][rec]
        n_ref = len(col(base + "-ann.txt"))
        for det in (tenhz.ALL_DETECTORS if a.fusion else tenhz.DETECTORS):
            t0 = time.perf_counter()
            loc = tenhz.qrs_locations(x, fs, det)
            dt = time.perf_counter() - t0
            np.savetxt(os.path.join(a.out, db, f"{rec}-{det}.txt"), loc, fmt="%d")
            rows.append(dict(db=db, record=rec, detector=det, beats=len(loc), ref_beats=n_ref,
                             hours=len(x) / fs / 3600, seconds=dt))
            print(f"{db}/{rec} {det:5s} {len(loc):8d} beats (reference {n_ref}), {dt:.2f} s", flush=True)
with open(os.path.join(a.out, "summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
print(f"wrote {a.out}/summary.csv")
