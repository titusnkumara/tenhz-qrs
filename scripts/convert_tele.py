#!/usr/bin/env python3
"""Convert the TELE ECG database (Khamis et al. 2016, Harvard Dataverse doi:10.7910/DVN/QTG0EP) to the
<db>-mat layout used by benchmark_three.py.

    python scripts/convert_tele.py --src <folder with the *.dat files> --dst <data>/tele-mat

Each .dat line is: ECG (mV), QRS flag, visual artifact mask, software artifact mask (500 Hz).
Writes <rec>.txt, <rec>-ann.txt (1-based beats), <rec>-mask.txt (visual artifact mask as 1-based
inclusive 'start end' intervals) and _manifest.json. The reference QRS annotations exist only outside the
visual mask, so benchmark_three.py ignores detections inside it (see load_mask)."""
import argparse, json, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--src", default="../tenhz-qrs-data/data/dataverse_files")
ap.add_argument("--dst", default="../tenhz-qrs-data/data/tele-mat")
a = ap.parse_args()
os.makedirs(a.dst, exist_ok=True)
manifest = {}
for fname in sorted(os.listdir(a.src)):
    if not fname.endswith(".dat"):
        continue
    rec = fname[:-4]
    df = pd.read_csv(os.path.join(a.src, fname), header=None)
    if df.shape[1] < 3:
        print("skip (unexpected columns):", fname); continue
    ecg = df.iloc[:, 0].to_numpy(float)
    ann = np.flatnonzero(df.iloc[:, 1].to_numpy(float) > 0.5) + 1
    m = np.r_[0, (df.iloc[:, 2].to_numpy(float) > 0.5).astype(np.int8), 0]
    d = np.diff(m)
    iv = np.c_[np.flatnonzero(d == 1) + 1, np.flatnonzero(d == -1)]
    np.savetxt(os.path.join(a.dst, f"{rec}.txt"), ecg[:, None], fmt="%.6f")
    np.savetxt(os.path.join(a.dst, f"{rec}-ann.txt"), ann[:, None], fmt="%d")
    np.savetxt(os.path.join(a.dst, f"{rec}-mask.txt"), iv.reshape(-1, 2), fmt="%d")
    manifest[rec] = {"fs": 500.0, "status": "OK", "n_samples": int(len(ecg)), "n_beats": int(len(ann)),
                     "masked_fraction": round(float(m.mean()), 4)}
json.dump(manifest, open(os.path.join(a.dst, "_manifest.json"), "w"), indent=1)
print("wrote", len(manifest), "records to", a.dst)
