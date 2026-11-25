#!/usr/bin/env python3
"""F1 against noise level on the MIT-BIH Noise Stress Test Database, for every detector in the benchmark.

    python scripts/plot_nstdb_snr.py --results results/benchmark_final --out figures

NSTDB records 118 and 119 with electrode-motion noise added at 24, 18, 12, 6, 0 and -6 dB SNR
(record names 118e24 ... 118e_6). F1 is pooled over the two base records at each noise level, with the
benchmark's matching (raw +/-150 ms by default). The three 10 Hz detectors are drawn in colour, the other
detectors in grey, and the dashed line is the best F1 any detector reaches at that noise level.

Writes <out>/fig_nstdb_snr.png and .pdf, and <out>/nstdb_snr.csv with the plotted numbers.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import competitors  # noqa: E402  (LABEL)

OURS = {"fusion": ("#0b5d2a", "o", 2.8), "minel": ("#2e9e4f", "s", 2.0), "titen": ("#7bc47f", "^", 2.0)}


def snr_of(rec):
    """118e24 -> 24, 118e06 -> 6, 118e00 -> 0, 118e_6 -> -6; None for any other record name."""
    m = re.fullmatch(r"(?:118|119)e(_?)(\d+)", rec)
    if not m:
        return None
    return -int(m.group(2)) if m.group(1) else int(m.group(2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="folder with per_record.csv")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--variant", default="raw", choices=["raw", "offset_corrected"])
    ap.add_argument("--tol", type=float, default=150.0, help="matching tolerance in ms (default 150)")
    a = ap.parse_args()

    R = pd.read_csv(os.path.join(a.results, "per_record.csv"), dtype={"record": str}, keep_default_na=False,
                    low_memory=False)
    g = R[(R.db == "nstdb") & (R.variant == a.variant) & (R.tol_ms == a.tol)].copy()
    g["snr"] = [snr_of(r) for r in g.record]
    g = g[g.snr.notna()]
    if g.empty:
        sys.exit("no NSTDB records 118e*/119e* in per_record.csv")
    s = g.groupby(["method", "snr"])[["tp", "fp", "fn"]].sum()
    s["f1"] = 200 * s.tp / (2 * s.tp + s.fp + s.fn)
    F = s.f1.unstack("snr").sort_index(axis=1, ascending=False)           # methods x SNR (24 ... -6)
    snrs = list(F.columns)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(a.out, exist_ok=True)
    F.rename(index=lambda m: competitors.LABEL.get(m, m)).round(2).to_csv(os.path.join(a.out, "nstdb_snr.csv"))
    best = F.max(axis=0)
    others = [m for m in F.index if m not in OURS]
    x = np.arange(len(snrs))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), sharex=True)
    for ax, (lo, hi), title in zip(axes, [(0, 101), (70, 100.6)], ["All detectors", "Zoom: 70-100%"]):
        for m in others:
            ax.plot(x, F.loc[m], color="0.7", lw=1.0, alpha=0.9, zorder=1)
        ax.plot(x, best, color="black", lw=1.4, ls="--", zorder=2, label="best of all detectors")
        for m, (c, mk, lw) in OURS.items():
            if m in F.index:
                ax.plot(x, F.loc[m], color=c, lw=lw, marker=mk, ms=6, zorder=3,
                        label=competitors.LABEL[m].replace(" (10 Hz, exploratory)", " (10 Hz)"))
        ax.plot([], [], color="0.7", lw=1.2, label=f"{len(others)} other detectors")
        ax.set_ylim(lo, hi)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{v:+d}" if v < 0 else f"{v}" for v in snrs])
        ax.set_xlabel("SNR of added electrode-motion noise (dB)  -  noisier to the right")
        ax.set_title(title, loc="left", fontsize=10.5)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(f"F1 (%), {a.variant.replace('_', ' ')} +/-{a.tol:g} ms, pooled over records 118 and 119")
    axes[1].legend(loc="lower left", fontsize=8.5, frameon=True)
    n_det = len(F.index)
    fig.suptitle(f"NSTDB: F1 against noise level, {n_det} detectors", fontsize=12.5, y=0.99)
    n_beats = int(g[g.method == "fusion"].groupby("snr").apply(lambda d: (d.tp + d.fn).sum(), include_groups=False).iloc[0])
    fig.text(0.01, 0.005, "MIT-BIH Noise Stress Test Database: records 118 and 119 with added electrode-motion noise. "
             f"Each point pools the two records at that SNR ({n_beats:,} beats).", fontsize=8, color="0.3")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"fig_nstdb_snr.{ext}"), dpi=200)
    print(f"saved fig_nstdb_snr.png / .pdf and nstdb_snr.csv in {a.out}")
    show = F.loc[[m for m in OURS if m in F.index]].rename(index=lambda m: competitors.LABEL[m])
    show.loc["best of all"] = best
    print(show.round(1).to_string())


if __name__ == "__main__":
    main()
