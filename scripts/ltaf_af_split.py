#!/usr/bin/env python3
"""Score saved detections on LTAFDB separately inside atrial-fibrillation and non-AF rhythm segments.

    python scripts/ltaf_af_split.py --data ../tenhz-qrs/data --run results/ltaf_v1 --out results/ltaf_af_v1

Needs
  * a benchmark_three.py run on LTAFDB that saved its detections (default; <run>/detections/ltafdb/),
  * <data>/ltafdb-mat/ (signal length, reference beats, _manifest.json) and <data>/ltafdb-meta/<rec>.csv
    (rhythm annotations written by scripts/download_physionet.py).

Rhythm segments come from the rhythm-change annotations ('(AFIB', '(N', '(B', '(T', ...): a label holds from its
sample to the next rhythm annotation or the end of the record. AF = the labels in --af-labels (default '(AFIB' only;
atrial flutter '(AFL' is counted as non-AF unless you add it). Reference beats and detections are restricted to the
segments of one class and matched there (one-to-one, closest first, same rule as the main benchmark), so a beat
detected within one tolerance of a segment boundary on the wrong side counts as one miss and one false alarm;
with segments of minutes this is negligible.

Offset-corrected variant: each record's detections are shifted by minus their median offset to the reference (estimated
within +/-150 ms over the whole record, not per class) before scoring, as in the main benchmark: an oracle diagnostic.

The reference is the PocketECG algorithm with manual verification by technicians, not a purely manual annotation.
Outputs in --out (must be new): ltaf_af.md, ltaf_af.csv (pooled), ltaf_af_per_record.csv.
"""
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import benchmark_three as B  # noqa: E402
import competitors  # noqa: E402

CLASSES = ["nonAF", "AF"]


def read_ann(data, rec, n):
    ann = pd.read_csv(os.path.join(data, "ltafdb-mat", rec + "-ann.txt"), header=None, sep=r"\s+", usecols=[0],
                      engine="c").to_numpy(dtype=float, copy=True).ravel()
    ann = np.round(ann).astype(np.int64) - 1
    return np.sort(ann[(ann >= 0) & (ann < n)])


def segments(data, rec, n, af_labels):
    """{'AF': (starts, ends), 'nonAF': (starts, ends)} in 0-based samples, from the rhythm annotations."""
    f = os.path.join(data, "ltafdb-meta", rec + ".csv")
    # Parsed by hand, not with a CSV reader: the exporter writes the annotation symbol and the note unquoted, and
    # WFDB symbols include the double quote (comment annotations), which a CSV reader takes as the start of a string.
    rows = []
    with open(f, encoding="utf-8", errors="replace") as fh:
        next(fh)                                            # header
        for line in fh:
            p = line.rstrip("\r\n").split(",", 2)         # sample, symbol, note (the note may contain commas)
            if len(p) == 3 and p[2].startswith("(") and p[0].isdigit():
                rows.append((int(p[0]) - 1, p[2].strip()))
    m = pd.DataFrame(rows, columns=["s", "aux_note"])
    m = m.sort_values("s").drop_duplicates("s", keep="last")
    starts = m.s.to_numpy(); labels = m.aux_note.str.strip().to_numpy()
    ends = np.r_[starts[1:], n]
    out = {c: ([], []) for c in CLASSES}
    for s, e, lab in zip(starts, ends, labels):
        c = "AF" if lab in af_labels else "nonAF"
        if out[c][0] and out[c][1][-1] == s:        # merge neighbours of the same class
            out[c][1][-1] = e
        else:
            out[c][0].append(s); out[c][1].append(e)
    return {c: (np.asarray(v[0], np.int64), np.asarray(v[1], np.int64)) for c, v in out.items()}


def inside(x, seg):
    starts, ends = seg
    if len(starts) == 0:
        return x[:0]
    i = np.searchsorted(starts, x, "right") - 1
    return x[(i >= 0) & (x < ends[np.maximum(i, 0)])]


def one_record(args):
    data, run, rec, fs, n, methods, tols, correct, af_labels = args
    ann = read_ann(data, rec, n)
    seg = segments(data, rec, n, af_labels)
    hours = {c: float((seg[c][1] - seg[c][0]).sum()) / fs / 3600 for c in CLASSES}
    nref = {c: int(len(inside(ann, seg[c]))) for c in CLASSES}
    rows = []
    for m in methods:
        f = os.path.join(run, "detections", "ltafdb", f"{rec}__{m}.npy")
        if not os.path.exists(f):
            continue
        det = np.sort(np.load(f).astype(np.int64))
        variants = [("raw", det)]
        if correct:
            off = B.match(det, ann, B.mround(B.OFFSET_WINDOW_MS / 1000 * fs))[3]
            variants.append(("offset_corrected", np.clip(det - (0 if np.isnan(off) else int(round(off))), 0, None)))
        for v, d in variants:
            for c in CLASSES:
                d_in, a_in = inside(d, seg[c]), inside(ann, seg[c])
                for t in tols:
                    tp, fp, fn, _ = B.match(d_in, a_in, B.mround(t / 1000 * fs))
                    rows.append(dict(record=rec, method=m, variant=v, tol_ms=t, cls=c, tp=tp, fp=fp, fn=fn,
                                     hours=hours[c], n_ref=nref[c]))
    return rows


def f1(tp, fp, fn):
    """F1 in %, NaN where there is nothing to score (e.g. a resample without AF beats)."""
    den = 2 * tp + fp + fn
    return np.where(den > 0, 100 * 2 * tp / np.maximum(den, 1), np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default=os.environ.get("TENHZ_DATA") or "data")
    ap.add_argument("--run", required=True, help="benchmark_three.py output folder of an LTAFDB run")
    ap.add_argument("--out", required=True, help="new folder")
    ap.add_argument("--af-labels", nargs="+", default=["(AFIB"], help="rhythm labels counted as AF (default '(AFIB')")
    ap.add_argument("--tol-ms", nargs="+", type=float, default=[150.0, 75.0])
    ap.add_argument("--offset-correct", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--methods", nargs="+", default=None, help="default: every detector found in <run>/detections")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--boot", type=int, default=1000, help="bootstrap resamples of records for the CIs (0 = none)")
    a = ap.parse_args()
    if os.path.exists(a.out) and os.listdir(a.out):
        sys.exit(f"{a.out} exists and is not empty: give a new --out folder")
    det_dir = os.path.join(a.run, "detections", "ltafdb")
    if not os.path.isdir(det_dir):
        sys.exit(f"{det_dir} not found: run benchmark_three.py on ltafdb with detections saved (the default) first")
    os.makedirs(a.out, exist_ok=True)
    man = json.load(open(os.path.join(a.data, "ltafdb-mat", "_manifest.json")))
    found = sorted({f.split("__")[1][:-4] for f in os.listdir(det_dir) if f.endswith(".npy")})
    methods = competitors.expand(a.methods) if a.methods else found
    methods = [m for m in methods if m in found]
    recs = [r for r in sorted(man) if man[r].get("status") == "OK"
            and os.path.exists(os.path.join(a.data, "ltafdb-meta", r + ".csv"))
            and os.path.exists(os.path.join(det_dir, f"{r}__{methods[0]}.npy"))]
    print(f"{len(recs)} records, {len(methods)} detectors, AF labels {a.af_labels}", flush=True)
    jobs = [(a.data, a.run, r, float(man[r]["fs"]), int(man[r]["n_samples"]), methods, tuple(a.tol_ms),
             a.offset_correct, tuple(a.af_labels)) for r in recs]
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        for i, rr in enumerate(ex.map(one_record, jobs), 1):
            rows += rr
            if i % 10 == 0 or i == len(jobs):
                print(f"[{i}/{len(jobs)}]", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(a.out, "ltaf_af_per_record.csv"), index=False)

    hrs = R.drop_duplicates(["record", "cls"]).groupby("cls").agg(hours=("hours", "sum"), beats=("n_ref", "sum"))
    n_af = int((R[R.cls == "AF"].drop_duplicates("record").hours > 0).sum())
    rng = np.random.default_rng(1)
    pooled_rows, md = [], ["# LTAFDB: AF versus non-AF rhythm segments", "",
        f"{len(recs)} records. AF = {', '.join(a.af_labels)}; non-AF = all other labelled rhythm segments. "
        f"Hours / reference beats: non-AF {hrs.loc['nonAF', 'hours']:.0f} h / {int(hrs.loc['nonAF', 'beats'])}, "
        f"AF {hrs.loc['AF', 'hours']:.0f} h / {int(hrs.loc['AF', 'beats'])}; {n_af} records contain AF. "
        "Reference: PocketECG algorithm with manual verification by technicians. "
        "Pooled over records; CI from resampling records.", ""]
    for v in ["raw", "offset_corrected"]:
        for t in sorted(a.tol_ms, reverse=True):
            g = R[(R.variant == v) & (R.tol_ms == t)]
            if g.empty:
                continue
            recs_ = sorted(g.record.unique()); ri = {r: i for i, r in enumerate(recs_)}
            ms = [m for m in methods if (g.method == m).any()]
            T = np.zeros((len(recs_), len(ms), 2, 3))
            for r in g.itertuples():
                T[ri[r.record], ms.index(r.method), CLASSES.index(r.cls)] = (r.tp, r.fp, r.fn)
            point = T.sum(0)                                          # methods x class x (tp, fp, fn)
            pf = f1(point[..., 0], point[..., 1], point[..., 2])
            boots = []
            for _ in range(a.boot):
                s = T[rng.integers(0, len(recs_), len(recs_))].sum(0)
                boots.append(f1(s[..., 0], s[..., 1], s[..., 2]))
            boots = np.asarray(boots) if a.boot else None
            title = f"{'Offset-corrected (oracle)' if v != 'raw' else 'Raw'}, matching tolerance +/-{t:g} ms"
            md += [f"## {title}", "",
                   "| Detector | non-AF F1 [CI] | AF F1 [CI] | AF minus non-AF [CI] | AF Se | AF PPV | non-AF Se | non-AF PPV |",
                   "|---|---|---|---|---|---|---|---|"]
            order = [m for m in competitors.OURS if m in ms] + sorted([m for m in ms if m not in competitors.OURS],
                                                                        key=lambda m: -pf[ms.index(m), 1])
            ci = lambda arr: "" if boots is None else f" [{np.nanpercentile(arr, 2.5):.2f}, {np.nanpercentile(arr, 97.5):.2f}]"
            for m in order:
                j = ms.index(m)
                se = 100 * point[j, :, 0] / np.maximum(point[j, :, 0] + point[j, :, 2], 1)
                pp = 100 * point[j, :, 0] / np.maximum(point[j, :, 0] + point[j, :, 1], 1)
                bn = boots[:, j, 0] if boots is not None else None
                ba = boots[:, j, 1] if boots is not None else None
                md.append(f"| {competitors.LABEL.get(m, m)} | {pf[j, 0]:.2f}{ci(bn) if bn is not None else ''} | "
                          f"{pf[j, 1]:.2f}{ci(ba) if ba is not None else ''} | {pf[j, 1] - pf[j, 0]:+.2f}"
                          f"{ci(ba - bn) if bn is not None else ''} | {se[1]:.2f} | {pp[1]:.2f} | {se[0]:.2f} | {pp[0]:.2f} |")
                pooled_rows.append(dict(variant=v, tol_ms=t, detector=m, nonAF_F1=pf[j, 0], AF_F1=pf[j, 1],
                                        AF_Se=se[1], AF_PPV=pp[1], nonAF_Se=se[0], nonAF_PPV=pp[0]))
            md.append("")
    pd.DataFrame(pooled_rows).to_csv(os.path.join(a.out, "ltaf_af.csv"), index=False)
    open(os.path.join(a.out, "ltaf_af.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()