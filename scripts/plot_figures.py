#!/usr/bin/env python3
"""Generate the combined multi-panel publication figures from benchmark_three.py outputs.

Produces:
  1. fig_comprehensive_raw150.{png,pdf}      (3-panel: 8-DB F1 heatmap + RR error + Relative runtime)
  2. fig_comprehensive_offset75.{png,pdf}    (same 3-panel layout for offset-corrected +/-75 ms, if present)
  3. fig_speed_accuracy_cleansweep.{png,pdf} (2-panel: linear speed-accuracy + inset + paired clean sweep)

Usage:
    python scripts/plot_figures.py --results results/benchmark_v1 --out figures/
    python scripts/plot_figures.py --csv-dir results/benchmark_v1 --out-dir figures/

"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8.5,
    "axes.titlesize": 9.8,
    "axes.labelsize": 8.8,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.2,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

TRAIN_DBS = ["mitdb", "svdb"]
TEST_DBS = ["edb", "nsrdb", "qtdb", "twadb", "nstdb"]
HELDOUT_DBS = ["ltafdb", "tele"]
ALL_DBS = TRAIN_DBS + TEST_DBS + HELDOUT_DBS

DARK_GREEN = "#005a32"
MID_GREEN = "#1b7837"
MINT_TINT = "#e8f5ee"
COMP_GRAY = "#737a82"


def clean_name(label: str) -> str:
    s = str(label)
    s = s.replace(" (10 Hz, exploratory)", " (10 Hz)")
    s = s.replace("Matched filter (py-ecg-detectors)", "Matched filter (py-ecg)")
    s = s.replace("Energy + find_peaks (baseline)", "Energy + find_peaks")
    return s


def is_ours(label: str) -> bool:
    low = str(label).lower()
    return low.startswith(("fusion", "minel", "titen"))


def method_key_from_label(label: str, available_methods: set[str]) -> str:
    low = str(label).lower()
    for prefix in ("fusion", "minel", "titen"):
        if low.startswith(prefix) and prefix in available_methods:
            return prefix
    mapping = {
        "matched filter (py-ecg-detectors)": "pyecg_matched",
        "matched filter (py-ecg)": "pyecg_matched",
        "energy + find_peaks (baseline)": "energy_peaks",
        "energy + find_peaks": "energy_peaks",
        "wfdb xqrs": "xqrs",
        "wfdb gqrs": "gqrs",
    }
    if low in mapping and mapping[low] in available_methods:
        return mapping[low]
    tok = low.split()[0]
    if tok in available_methods:
        return tok
    return low


def find_fastest_reference(df: pd.DataFrame) -> tuple[int, float, str]:
    sph = df["s_per_h"].to_numpy(dtype=float)
    valid = np.isfinite(sph) & (sph > 0)
    if not np.any(valid):
        return 0, 1.0, clean_name(df.iloc[0]["detector"]).replace(" (10 Hz)", "")
    valid_indices = np.where(valid)[0]
    fastest_idx = int(valid_indices[np.argmin(sph[valid_indices])])
    fastest_sph = float(sph[fastest_idx])
    fastest_short = clean_name(df.iloc[fastest_idx]["detector"]).replace(" (10 Hz)", "")
    return fastest_idx, fastest_sph, fastest_short


def dataset_totals(per_record_path: str | None):
    """Records, hours and reference beats in the run, each record counted once (from per_record.csv)."""
    if not per_record_path or not os.path.exists(per_record_path):
        return None
    pr = pd.read_csv(per_record_path, dtype={"record": str}, low_memory=False)
    r = pr[(pr["variant"] == "raw")]
    r = r[r["tol_ms"] == r["tol_ms"].max()].drop_duplicates(["db", "record"])
    return len(r), float(r["hours"].sum()), int((r["tp"] + r["fn"]).sum())


def totals_text(tot) -> str:
    if tot is None:
        return ""
    n, h, b = tot
    return f"{n} Recordings, {h:,.0f} Hours, {b / 1e6:.1f} M Beats"


def holm_adjust(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.ones(m, dtype=float)
    running = 0.0
    for rank, idx in enumerate(order):
        val = pvals[idx] * (m - rank)
        running = max(running, val)
        adj[idx] = min(1.0, running)
    return adj.tolist()


def build_universal_f1_cmap():
    stops = [
        (0.00, "#0b0b0b"),
        (0.35, "#67000d"),
        (0.55, "#cb181d"),
        (0.72, "#ef6548"),
        (0.84, "#fec44f"),
        (0.91, "#a1d99b"),
        (0.96, "#31a354"),
        (1.00, "#00441b"),
    ]
    return mcolors.LinearSegmentedColormap.from_list("universal_f1", stops, N=512)


def draw_broken_slashes(ax_left, ax_right, y_positions=(0.0, 1.0), d=0.015):
    kwargs = dict(color="#222222", clip_on=False, linewidth=1.0)
    for y in y_positions:
        ax_left.plot((1 - d, 1 + d), (y - 2.2 * d, y + 2.2 * d), transform=ax_left.transAxes, **kwargs)
        ax_right.plot((-d, +d), (y - 2.2 * d, y + 2.2 * d), transform=ax_right.transAxes, **kwargs)


def compute_wilcoxon_stars(per_record_path: str, variant: str, tol_ms: float, detectors_ordered: list[str]):
    stars: dict[tuple[str, str], bool] = {}
    win_counts = {"Fusion": 0, "MinEl": 0, "TiTen": 0}
    total_pairs = 0
    if not per_record_path or not os.path.exists(per_record_path):
        print(f"  [wilcoxon] per_record not found, skipping stars: {per_record_path}")
        return stars, win_counts, total_pairs
    pr = pd.read_csv(per_record_path, dtype={"record": str}, low_memory=False)
    sub = pr[(pr["variant"] == variant) & (pr["tol_ms"] == tol_ms)].copy()
    if sub.empty:
        print(f"  [wilcoxon] no per_record rows for variant={variant} tol={tol_ms}")
        return stars, win_counts, total_pairs
    denom = 2 * sub["tp"] + sub["fp"] + sub["fn"]
    sub["rec_f1"] = np.where(denom > 0, 100.0 * 2 * sub["tp"] / denom.clip(lower=1), np.nan)
    avail = set(sub["method"].astype(str).str.lower())
    sub["m_low"] = sub["method"].astype(str).str.lower()
    piv = sub.pivot_table(index=["db", "record"], columns="m_low", values="rec_f1", aggfunc="first")
    if "fusion" not in piv.columns:
        print("  [wilcoxon] 'fusion' not in per_record methods, skipping stars")
        return stars, win_counts, total_pairs
    # Stars only where Fusion was not tuned (test + held-out); never on MITDB/SVDB (in-sample).
    # The pooled test column is starred from the paired table instead (see plot_three_panel_benchmark).
    db_groups = [(db, [db]) for db in TEST_DBS + HELDOUT_DBS]
    comp_dets = [d for d in detectors_ordered if not is_ours(d)]
    for col_key, dbs in db_groups:
        mask = piv.index.get_level_values(0).isin(dbs)
        block = piv[mask]
        if block.empty:
            continue
        p_raws = []
        valid_dets = []
        for det in comp_dets:
            mk = method_key_from_label(det, avail)
            if mk not in block.columns:
                continue
            pair = block[["fusion", mk]].dropna()
            diff = pair["fusion"].to_numpy() - pair[mk].to_numpy()
            nz = diff[np.abs(diff) > 1e-12]
            if len(nz) < 5:
                p = 1.0
            else:
                _, p = wilcoxon(nz, alternative="two-sided")
            p_raws.append(float(p))
            valid_dets.append((det, float(np.mean(diff)) if len(diff) else 0.0))
        if p_raws:
            p_adj = holm_adjust(p_raws)
            for (det, dm), pa in zip(valid_dets, p_adj):
                stars[(det, col_key)] = (pa < 0.05 and dm > 0)   # star = Fusion significantly BETTER
    unseen_dbs = TEST_DBS + HELDOUT_DBS
    total_pairs = len(unseen_dbs) * len(comp_dets)
    for our_lbl, our_key in [("Fusion", "fusion"), ("MinEl", "minel"), ("TiTen", "titen")]:
        if our_key not in piv.columns:
            continue
        wins = 0
        for db in unseen_dbs:
            block = piv[piv.index.get_level_values(0) == db]
            p_raws = []
            diffs_mean = []
            for det in comp_dets:
                mk = method_key_from_label(det, avail)
                if mk not in block.columns:
                    p_raws.append(1.0)
                    diffs_mean.append(0.0)
                    continue
                pair = block[[our_key, mk]].dropna()
                diff = pair[our_key].to_numpy() - pair[mk].to_numpy()
                diffs_mean.append(float(np.mean(diff)))
                nz = diff[np.abs(diff) > 1e-12]
                p = 1.0 if len(nz) < 5 else float(wilcoxon(nz, alternative="two-sided")[1])
                p_raws.append(p)
            p_adj = holm_adjust(p_raws)
            for pa, dm in zip(p_adj, diffs_mean):
                if pa < 0.05 and dm > 0:
                    wins += 1
        win_counts[our_lbl] = wins
    return stars, win_counts, total_pairs


# ===========================================================================
# FIGURE 1: 3-Panel Combined Benchmark Figure (Heatmap + RR Error + Runtime)
# ===========================================================================
def plot_three_panel_benchmark(
    summary_path: str,
    rr_path: str,
    per_record_path: str,
    out_dir: str,
    variant: str = "raw",
    tol_ms: float = 150.0,
    stem: str = "fig_comprehensive_raw150",
    paired_path: str | None = None,
):
    print(f"\n[three-panel] {stem}: variant={variant} tol={tol_ms}")
    print(f"  summary: {summary_path} exists={os.path.exists(summary_path) if summary_path else False}")
    print(f"  rr:      {rr_path} exists={os.path.exists(rr_path) if rr_path else False}")
    if not summary_path or not os.path.exists(summary_path):
        raise FileNotFoundError(f"summary/headline CSV not found: {summary_path}")
    if not rr_path or not os.path.exists(rr_path):
        raise FileNotFoundError(f"rr_error CSV not found: {rr_path} (need rr_error.csv in --results)")

    summ = pd.read_csv(summary_path)
    sub = summ[(summ["variant"] == variant) & (summ["tol_ms"] == tol_ms)].copy()
    sub = sub.dropna(subset=["test_F1"]).sort_values("test_F1", ascending=False).reset_index(drop=True)
    if sub.empty:
        print(f"  [three-panel] SKIP: no rows for variant={variant} tol={tol_ms} in {summary_path}")
        print(f"  available: {summ[['variant','tol_ms']].drop_duplicates().to_dict('records')}")
        return

    # FIX: normalise RR detector names so "Fusion (10 Hz, exploratory)" matches "fusion" / "Fusion"
    # Strip parentheticals and add short-key aliases (pyecg_matched, energy) on both sides.
    def _norm_det(s: str) -> str:
        v = clean_name(s).lower().replace(" (10 hz)", "").replace(" (10 hz, exploratory)", "").strip()
        v = v.replace("matched filter (py-ecg)", "matched filter").replace("energy + find_peaks", "energy")
        if v.startswith("wfdb "):
            v = v[5:]
        if v in ("pyecg_matched", "pyecg", "matched filter (py-ecg-detectors)"):
            v = "matched filter"
        if v in ("energy_peaks", "energy + find_peaks (baseline)"):
            v = "energy"
        return v
    rr_raw = pd.read_csv(rr_path)
    rr_by_clean = {}
    for _, r in rr_raw.iterrows():
        rr_by_clean[_norm_det(r["detector"])] = r
        rr_by_clean[clean_name(r["detector"]).lower()] = r
        rr_by_clean[str(r["detector"]).lower()] = r
        # aliases used in paired / per_record files
        if str(r["detector"]).lower() == "pyecg_matched":
            rr_by_clean["matched filter"] = r
        if str(r["detector"]).lower() == "energy":
            rr_by_clean["energy"] = r
    print(f"  detectors: {len(sub)}, rr rows: {len(rr_raw)}, rr matched: {sum(1 for d in sub['detector'] if _norm_det(d) in rr_by_clean)}/{len(sub)}")

    rr = rr_raw.set_index("detector")
    n_det = len(sub)
    detectors = sub["detector"].tolist()
    fastest_idx, fastest_sph, fastest_short = find_fastest_reference(sub)
    rel_rt = sub["s_per_h"].to_numpy(dtype=float) / max(fastest_sph, 1e-9)
    stars, win_counts, total_pairs = compute_wilcoxon_stars(per_record_path, variant, tol_ms, detectors)
    # Pooled test column: same rule as the paired table (bootstrap CI excludes 0 AND Holm p < 0.05).
    if paired_path and os.path.exists(paired_path):
        pp = pd.read_csv(paired_path)
        pp = pp[pp["reference"] == "fusion"]
        avail_keys = set(pp["comparator"].astype(str))
        for det in detectors:
            if is_ours(det):
                continue
            mk = method_key_from_label(det, avail_keys)
            row = pp[pp["comparator"] == mk]
            if len(row):
                stars[(det, "test_F1")] = str(row.iloc[0]["verdict"]) == "reference better"
    no_det = [d for d, se, ppv in zip(detectors, sub["test_Se"], sub["test_PPV"]) if se == 0 and ppv == 0]

    fig = plt.figure(figsize=(16.8, max(8.4, 0.43 * n_det)))
    gs = gridspec.GridSpec(1, 5, width_ratios=[4.3, 1.15, 0.44, 1.15, 0.48], wspace=0.08)
    ax_hm = fig.add_subplot(gs[0, 0])
    ax_rr_l = fig.add_subplot(gs[0, 1])
    ax_rr_r = fig.add_subplot(gs[0, 2])
    ax_rt_l = fig.add_subplot(gs[0, 3])
    ax_rt_r = fig.add_subplot(gs[0, 4])
    y_idx = np.arange(n_det)

    cols = [c for c in ALL_DBS if c in sub.columns] + ["test_F1"]
    col_headers = ["MITDB\n(train)", "SVDB\n(train)", "EDB\n(test)", "NSRDB\n(test)", "QTDB\n(test)", "TWADB\n(test)", "NSTDB\n(test)", "LTAFDB\n(held-out)", "TELE\n(held-out)", "POOLED\nTEST"]
    col_headers = [h for h, c in zip(col_headers, ALL_DBS + ["test_F1"]) if c in cols]
    mat = sub[cols].to_numpy(dtype=float)
    cmap = build_universal_f1_cmap()
    ax_hm.imshow(mat / 100.0, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
    for i, det in enumerate(detectors):
        ours = is_ours(det)
        for j, col in enumerate(cols):
            val = mat[i, j]
            if np.isnan(val):
                txt = "-"
                tcol = "#ffffff"
            else:
                star = "*" if stars.get((det, col), False) else ""
                txt = f"{val:.2f}{star}"
                tcol = "#ffffff" if (val >= 94.0 or val < 68.0) else "#111111"
            ax_hm.text(j, i, txt, ha="center", va="center", fontsize=7.3, fontweight="bold" if ours else "normal", color=tcol)
    if "ltafdb" in cols:
        ltaf_j = cols.index("ltafdb")
        ax_hm.axvline(ltaf_j - 0.5, color="white", linestyle="--", linewidth=1.4)
        last_j = cols.index("tele") if "tele" in cols else ltaf_j
        ax_hm.axvline(last_j + 0.5, color="white", linestyle="-", linewidth=2.2)
    ax_hm.set_xticks(np.arange(len(cols)))
    ax_hm.set_xticklabels(col_headers, fontsize=7.5, fontweight="bold", color="#222222")
    ax_hm.xaxis.tick_top()
    ax_hm.set_yticks(y_idx)
    ranked_labels = [f"{i + 1}. {clean_name(d)}" + (" \u2020" if d in no_det else "") for i, d in enumerate(detectors)]
    ax_hm.set_yticklabels(ranked_labels, fontsize=8.3)
    for tick, det in zip(ax_hm.get_yticklabels(), detectors):
        if is_ours(det):
            tick.set_color(MID_GREEN)
            tick.set_fontweight("bold")
        else:
            tick.set_color("#222222")
    title_tag = "(* = Fusion significantly better)" if stars else ""
    ax_hm.set_title(f"A. F1 Across All 9 Databases & Pooled Test {title_tag}", fontsize=9.8, pad=26, loc="left")
    cbar_ax = ax_hm.inset_axes([0.02, -0.055, 0.96, 0.018])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mcolors.Normalize(0, 100), cmap=cmap), cax=cbar_ax, orientation="horizontal")
    cb.set_ticks([0, 35, 55, 72, 84, 91, 96, 100])
    cb.set_ticklabels(["0%", "35%", "55%", "72%", "84%", "91%", "96%", "100%"])
    cb.ax.tick_params(labelsize=7.2)
    cb.set_label("Universal F1 Color Scale (0% = Black, 100% = Dark Green)", fontsize=7.6)

    def prep_broken_pair(ax_l, ax_r, xlim_l, xlim_r, xticks_l, xticks_r):
        for ax in (ax_l, ax_r):
            ax.set_ylim(n_det - 0.5, -0.5)
            ax.set_yticks(y_idx)
            ax.set_yticklabels([])
            ax.spines["top"].set_visible(True)
            ax.grid(True, axis="x", linestyle="--", linewidth=0.5, alpha=0.45)
            ax.axhspan(-0.5, 2.5, color=MINT_TINT, zorder=0)
            for y_sep in (2.5, 6.5):
                if y_sep < n_det - 0.5:
                    ax.axhline(y_sep, color="#bbbbbb", linestyle="--", linewidth=0.7, zorder=1)
        ax_l.set_xlim(*xlim_l)
        ax_r.set_xlim(*xlim_r)
        ax_l.set_xticks(xticks_l)
        ax_r.set_xticks(xticks_r)
        ax_l.spines["right"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        ax_r.tick_params(axis="y", length=0)
        draw_broken_slashes(ax_l, ax_r)

    # Panel B: RR Timing Error, using normalised lookup
    rr_means = []
    rr_rmses = []
    for d in detectors:
        row = rr_by_clean.get(_norm_det(d))
        if row is None:
            rr_means.append(np.nan)
            rr_rmses.append(np.nan)
        else:
            rr_means.append(float(row["mean"]))
            rr_rmses.append(float(row["rmse"]))
    rr_means = np.array(rr_means, dtype=float)
    rr_rmses = np.array(rr_rmses, dtype=float)
    split_rr = 46.0
    max_rr = max(165.0, np.nanmax(rr_rmses) * 1.38) if np.any(np.isfinite(rr_rmses)) else 165.0
    prep_broken_pair(ax_rr_l, ax_rr_r, (0, split_rr), (68, max_rr), [0, 15, 30], [80, 120])
    for i, det in enumerate(detectors):
        mae, rmse = rr_means[i], rr_rmses[i]
        if np.isnan(mae):
            continue
        col = DARK_GREEN if is_ours(det) else COMP_GRAY
        dcol = DARK_GREEN if is_ours(det) else "#111111"
        for ax in (ax_rr_l, ax_rr_r):
            ax.barh(i, mae, height=0.56, color=col, alpha=0.92, zorder=3)
            ax.scatter(rmse, i, marker="D", s=26, color=dcol, zorder=4)
        lbl_txt = f"{mae:.1f} / {rmse:.1f}"
        target_ax = ax_rr_r if rmse > split_rr else ax_rr_l
        target_ax.text(rmse + (2.0 if rmse > split_rr else 1.0), i, lbl_txt, va="center", ha="left", fontsize=7.0, color=DARK_GREEN if is_ours(det) else "#333333", fontweight="bold" if is_ours(det) else "normal")
    ax_rr_l.set_title("B. RR Timing Error (ms)", fontsize=9.8, pad=26, loc="left")
    ax_rr_l.set_xlabel("Test RR Error: MAE / RMSE (ms)", fontsize=8.3, x=0.72)
    leg_handles = [mpl.lines.Line2D([], [], color="#111111", marker="D", linestyle="None", markersize=5, label="RMSE (ms)"), patches.Patch(facecolor=DARK_GREEN, label="MAE (ms)")]
    ax_rr_l.legend(handles=leg_handles, loc="lower right", bbox_to_anchor=(0.98, 0.05), fontsize=7.0, frameon=True)

    # Panel C: Relative Runtime
    split_rt = 95.0
    max_rt_val = float(np.nanmax(rel_rt)) if np.any(np.isfinite(rel_rt)) else 1800.0
    lo_r = split_rt + 15.0
    max_rt = max(lo_r + 50.0, max_rt_val * 1.25)
    r_ticks = [t for t in mpl.ticker.MaxNLocator(3).tick_values(lo_r, max_rt) if lo_r <= t <= max_rt]
    prep_broken_pair(ax_rt_l, ax_rt_r, (0, split_rt), (lo_r, max_rt), [0, 25, 50, 75], r_ticks)
    ax_rt_l.axvline(1.0, color=MID_GREEN, linestyle="--", linewidth=0.9, zorder=2)
    for i, (det, rt) in enumerate(zip(detectors, rel_rt)):
        if np.isnan(rt):
            continue
        ours = is_ours(det)
        col = DARK_GREEN if ours else COMP_GRAY
        for ax in (ax_rt_l, ax_rt_r):
            ax.barh(i, rt, height=0.56, color=col, alpha=0.92, zorder=3)
        if i == fastest_idx:
            lbl = "1.0x (fastest ref)"
        elif rt < 10.0:
            lbl = f"{rt:.1f}x"
        else:
            lbl = f"{rt:,.0f}x"
        target_ax = ax_rt_r if rt > split_rt else ax_rt_l
        target_ax.text(rt + (max_rt * 0.018 if rt > split_rt else 1.8), i, lbl, va="center", ha="left", fontsize=7.1, color=MID_GREEN if ours else "#333333", fontweight="bold" if ours else "normal")
    ax_rt_l.set_title(f"C. Relative Runtime (1.0x = {fastest_short})", fontsize=9.8, pad=26, loc="left")
    ax_rt_l.set_xlabel(f"Relative Runtime (x {fastest_short}, linear)", fontsize=8.3, x=0.72)

    mode_str = f"Raw +/-{tol_ms:g} ms" if variant == "raw" else f"Offset-Adjusted +/-{tol_ms:g} ms"
    tot = totals_text(dataset_totals(per_record_path))
    fig.suptitle(f"QRS Detector Benchmark: {n_det} Detectors, {tot + ', ' if tot else ''}{mode_str}", fontsize=12.2, y=0.985)
    notes = ["Per-database stars: per-record Wilcoxon, Holm-corrected, Fusion better. "
             "Pooled-test stars: paired bootstrap CI excludes 0 and Holm p < 0.05 (as in the paired table)."]
    if no_det:
        notes.append("\u2020 produced no detections on any test record with its library defaults.")
    fig.text(0.01, -0.045, "  ".join(notes), ha="left", va="top", fontsize=7.0, color="#444444")
    if total_pairs > 0:
        badge = (f"Unseen-DB Win Rates vs. {n_det - 3} Competing Detectors ({len(TEST_DBS + HELDOUT_DBS)} Unseen DBs = {total_pairs} Pairs, Paired Wilcoxon $p_{{\\mathrm{{Holm}}}} < 0.05$):   Fusion: {win_counts['Fusion']}/{total_pairs} ({100*win_counts['Fusion']/total_pairs:.1f}%)  |  MinEl: {win_counts['MinEl']}/{total_pairs} ({100*win_counts['MinEl']/total_pairs:.1f}%)  |  TiTen: {win_counts['TiTen']}/{total_pairs} ({100*win_counts['TiTen']/total_pairs:.1f}%)")
        fig.text(0.01, -0.012, badge, ha="left", va="top", fontsize=8.2, color=DARK_GREEN, bbox=dict(boxstyle="round,pad=0.35", facecolor=MINT_TINT, edgecolor=MID_GREEN, linewidth=0.9))
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"{stem}.png"), dpi=300)
    fig.savefig(os.path.join(out_dir, f"{stem}.pdf"))
    plt.close(fig)
    print(f"  Saved {stem}.png and {stem}.pdf -> {out_dir}/")


# ===========================================================================
# FIGURE 2: 2-Panel Speed-Accuracy Trade-off + Record-Level Clean Sweep
# ===========================================================================
def plot_speed_accuracy_cleansweep(
    summary_path: str,
    paired_test_path: str,
    paired_ltaf_path: str,
    per_record_path: str,
    out_dir: str,
    stem: str = "fig_speed_accuracy_cleansweep",
):
    print(f"\n[clean-sweep] {stem}")
    print(f"  summary: {summary_path} exists={os.path.exists(summary_path) if summary_path else False}")
    print(f"  paired test: {paired_test_path} exists={os.path.exists(paired_test_path) if paired_test_path else False}")
    print(f"  paired ltaf: {paired_ltaf_path} exists={os.path.exists(paired_ltaf_path) if paired_ltaf_path else False}")
    if not summary_path or not os.path.exists(summary_path):
        raise FileNotFoundError(f"summary/headline CSV not found: {summary_path}")
    if not paired_test_path or not os.path.exists(paired_test_path):
        raise FileNotFoundError(f"paired test CSV not found: {paired_test_path} (need paired_raw_150.csv in --results)")

    summ = pd.read_csv(summary_path)
    sub = summ[(summ["variant"] == "raw") & (summ["tol_ms"] == 150.0)].dropna(subset=["test_F1", "s_per_h"]).copy()
    sub = sub.sort_values("test_F1", ascending=False).reset_index(drop=True)
    if sub.empty:
        print(f"  [clean-sweep] SKIP: no raw/150 rows in {summary_path}")
        return
    fastest_idx, fastest_sph, fastest_short = find_fastest_reference(sub)
    sub["rel_rt"] = sub["s_per_h"] / max(fastest_sph, 1e-9)

    fig = plt.figure(figsize=(16.8, 9.2))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.22, 1.02, 0.34], wspace=0.40)
    ax_sc = fig.add_subplot(gs[0, 0])
    ax_sw_l = fig.add_subplot(gs[0, 1])
    ax_sw_r = fig.add_subplot(gs[0, 2])

    in_range = sub[(sub["rel_rt"] <= 152.0) & (sub["test_F1"] >= 82.0)].copy()
    off_scale = sub[(sub["rel_rt"] > 152.0) | (sub["test_F1"] < 82.0)].copy()
    ours_df = in_range[in_range["detector"].apply(is_ours)].sort_values("rel_rt")
    comp_df = in_range[~in_range["detector"].apply(is_ours)]
    ax_sc.axvline(1.0, color="#333333", linestyle=":", linewidth=1.0, zorder=1)
    ax_sc.plot(ours_df["rel_rt"], ours_df["test_F1"], color="#111111", linewidth=1.5, zorder=3, label="Pareto frontier (10 Hz family)")
    for idx, r in ours_df.iterrows():
        is_fastest = (idx == fastest_idx)
        ax_sc.scatter(r["rel_rt"], r["test_F1"], s=115 if is_fastest else 85, color=DARK_GREEN if is_fastest else MID_GREEN, edgecolors="#111111", linewidths=1.0, zorder=5, label=f"{fastest_short} (10 Hz, 1.0x fastest ref)" if is_fastest else None)
    if len(ours_df) > 1:
        ax_sc.scatter([], [], s=85, color=MID_GREEN, edgecolors="#111111", linewidths=1.0, label="Other 10 Hz detectors")
    ax_sc.scatter(comp_df["rel_rt"], comp_df["test_F1"], s=36, color=COMP_GRAY, alpha=0.85, zorder=3, label=f"External baselines ({len(sub) - len(ours_df)} methods)")
    for _, r in comp_df.iterrows():
        lbl = clean_name(r["detector"])
        if r["test_F1"] < 99.15:
            dx, dy = (5, -2)
            if "elgendi2013" in lbl.lower():
                dy = -7
            elif "matched filter" in lbl.lower():
                dy = -8
            elif "gqrs" in lbl.lower():
                dx, dy = (-44, -9)
            ax_sc.annotate(lbl, (r["rel_rt"], r["test_F1"]), xytext=(dx, dy), textcoords="offset points", fontsize=7.1, color="#333333")
    slowest_idx = int(sub["rel_rt"].idxmax())
    sr = sub.loc[slowest_idx]
    if sr["rel_rt"] > 152.0:
        ax_sc.annotate(f"{clean_name(sr['detector'])}: {sr['test_F1']:.2f}% @ {sr['rel_rt']:,.0f}x slower ->", xy=(149.0, sr["test_F1"] + 0.15), xytext=(74.0, 98.0), fontsize=7.5, color="#333333", bbox=dict(boxstyle="round,pad=0.3", facecolor="#f7f7f7", edgecolor="#bbbbbb"), arrowprops=dict(arrowstyle="->", color="#222222", lw=1.0))
    ax_sc.legend(loc="center right", bbox_to_anchor=(0.985, 0.72), fontsize=7.3, frameon=True)
    ax_ins = ax_sc.inset_axes([0.21, 0.14, 0.52, 0.36])
    elite = sub[(sub["rel_rt"] <= 48.0) & (sub["test_F1"] >= 99.28)].copy()
    elite_ours = elite[elite["detector"].apply(is_ours)].sort_values("rel_rt")
    elite_comp = elite[~elite["detector"].apply(is_ours)]
    ax_ins.plot(elite_ours["rel_rt"], elite_ours["test_F1"], color="#111111", linewidth=1.2, zorder=3)
    for idx, r in elite_ours.iterrows():
        is_fastest = (idx == fastest_idx)
        is_fusion = str(r["detector"]).lower().startswith("fusion")
        ax_ins.scatter(r["rel_rt"], r["test_F1"], s=72 if (is_fastest or is_fusion) else 52, color=DARK_GREEN if is_fastest else MID_GREEN, edgecolors="#111111", zorder=5)
        dy = 2 if is_fusion else (-7 if "MinEl" in r["detector"] else -9)
        ax_ins.annotate(f"{clean_name(r['detector'])} ({r['test_F1']:.2f}%, {r['rel_rt']:.1f}x)", (r["rel_rt"], r["test_F1"]), xytext=(6, dy), textcoords="offset points", fontsize=6.5, color=MID_GREEN, fontweight="bold")
    for _, r in elite_comp.iterrows():
        ax_ins.scatter(r["rel_rt"], r["test_F1"], s=34, color="#444444", edgecolors="#111111", zorder=4)
        ax_ins.annotate(f"{clean_name(r['detector'])} ({r['test_F1']:.2f}%)", (r["rel_rt"], r["test_F1"]), xytext=(0, -10), textcoords="offset points", ha="center", fontsize=6.4, color="#222222")
    ax_ins.set_xlim(-2.0, 46.0)
    ax_ins.set_ylim(99.28, 99.67)
    ax_ins.set_title(f"Linear Zoom: Elite Regime (<45x {fastest_short}, >99.3% F1)", fontsize=7.2, pad=4)
    ax_ins.tick_params(labelsize=6.8)
    ax_ins.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    if not off_scale.empty:
        parts = [f"{clean_name(r['detector'])} ({r['test_F1']:.2f}%, {r['rel_rt']:,.0f}x)" for _, r in off_scale.iterrows()]
        mid = (len(parts) + 1) // 2
        line1 = "Off-scale slower/lower-F1:  " + "  |  ".join(parts[:mid])
        line2 = "  |  ".join(parts[mid:]) if len(parts) > mid else ""
        banner_txt = f"{line1}\n{line2}" if line2 else line1
        ax_sc.text(0.02, 0.018, banner_txt, transform=ax_sc.transAxes, fontsize=6.3, color="#8b0000", linespacing=1.35, bbox=dict(boxstyle="round,pad=0.28", facecolor="#fdecea", edgecolor="#e5989b", linewidth=0.8))
    ax_sc.set_xlim(-5.0, 152.0)
    ax_sc.set_ylim(82.0, 100.5)
    ax_sc.set_xlabel(f"Relative Compute Time (x {fastest_short}, linear scale, lower is faster)")
    ax_sc.set_ylabel("Pooled Test Set F1 (%, raw +/-150 ms)")
    ax_sc.set_title(f"A. Accuracy vs. Relative Speed (All {len(sub)} Detectors, Linear)", fontsize=10.2, loc="left", pad=10)
    ax_sc.grid(True, linestyle="--", linewidth=0.5, alpha=0.45)

    pt = pd.read_csv(paired_test_path)
    pt = pt[pt["reference"] == "fusion"].copy()
    pl = pd.read_csv(paired_ltaf_path) if paired_ltaf_path and os.path.exists(paired_ltaf_path) else pd.DataFrame()
    if not pl.empty:
        pl = pl[pl["reference"] == "fusion"].set_index("comparator")
    # FIX: paired files use short keys (fusion, minel, xqrs, pyecg_matched)
    # while summary uses long names (Fusion (10 Hz, exploratory), WFDB xqrs).
    # Build short->display map so ranks and labels line up.
    _short_to_display = {}
    for _d in sub["detector"]:
        _c = clean_name(_d)
        _short_to_display[_c.lower()] = _c
        _short_to_display[_c.lower().split()[0]] = _c
        _short_to_display[_c.lower().replace("wfdb ", "")] = _c
    _short_to_display.update({"pyecg_matched": "Matched filter (py-ecg)", "energy": "Energy + find_peaks", "energy_peaks": "Energy + find_peaks"})
    pt["comparator"] = pt["comparator"].apply(lambda c: _short_to_display.get(str(c).lower(), str(c)))
    if not pl.empty:
        pl.index = [_short_to_display.get(str(c).lower(), str(c)) for c in pl.index]
    rank_map = {clean_name(d): idx + 1 for idx, d in enumerate(sub["detector"])}
    pt["comp_clean"] = pt["comparator"].apply(clean_name)
    pt["rank"] = pt["comp_clean"].map(lambda c: rank_map.get(c, 99))
    pt = pt.sort_values("rank", ascending=True).reset_index(drop=True)

    wtl_map: dict[str, str] = {}
    if per_record_path and os.path.exists(per_record_path):
        pr = pd.read_csv(per_record_path, dtype={"record": str}, low_memory=False)
        pr_t = pr[(pr["variant"] == "raw") & (pr["tol_ms"] == 150.0) & (pr["db"].isin(TEST_DBS))].copy()
        denom = 2 * pr_t["tp"] + pr_t["fp"] + pr_t["fn"]
        pr_t["rec_f1"] = np.where(denom > 0, 100.0 * 2 * pr_t["tp"] / denom.clip(lower=1), np.nan)
        avail = set(pr_t["method"].astype(str).str.lower())
        piv = pr_t.pivot_table(index=["db", "record"], columns="method", values="rec_f1", aggfunc="first")
        if "fusion" in piv.columns:
            for comp_raw in pt["comparator"]:
                mk = method_key_from_label(comp_raw, avail)
                if mk in piv.columns:
                    d = (piv["fusion"] - piv[mk]).dropna().to_numpy()
                    w = int(np.sum(d > 1e-6))
                    l = int(np.sum(d < -1e-6))
                    t = int(len(d) - w - l)
                    wtl_map[comp_raw] = f"({w}/{t}/{l})"
    n_comp = len(pt)
    y_c = np.arange(n_comp)
    split_x = 18.0
    max_x = max(55.0, pt["ci_hi"].max() * 1.08)
    for ax in (ax_sw_l, ax_sw_r):
        ax.set_ylim(n_comp - 0.5, -0.5)
        ax.set_yticks(y_c)
        ax.set_yticklabels([])
        ax.grid(True, axis="x", linestyle="--", linewidth=0.5, alpha=0.45)
        ax.axhspan(-0.5, 1.5, color=MINT_TINT, zorder=0)
        ax.axhline(1.5, color="#bbbbbb", linestyle="--", linewidth=0.7)
    ax_sw_l.set_xlim(-0.6, split_x)
    ax_sw_r.set_xlim(21.0, max_x)
    ax_sw_l.set_xticks([0, 3, 6, 9, 12, 15])
    ax_sw_r.set_xticks([t for t in mpl.ticker.MaxNLocator(3).tick_values(21.0, max_x) if 21.0 <= t <= max_x])
    ax_sw_l.spines["right"].set_visible(False)
    ax_sw_r.spines["left"].set_visible(False)
    ax_sw_r.tick_params(axis="y", length=0)
    ax_sw_l.axvline(0.0, color="#111111", linestyle="--", linewidth=1.0, zorder=2)
    draw_broken_slashes(ax_sw_l, ax_sw_r)
    for i, r in pt.iterrows():
        comp_raw = r["comparator"]
        for ax in (ax_sw_l, ax_sw_r):
            ax.errorbar(r["diff"], i - 0.14, xerr=[[max(0.0, r["diff"] - r["ci_lo"])], [max(0.0, r["ci_hi"] - r["diff"])]], fmt="o", color="#111111", ecolor="#111111", elinewidth=1.2, capsize=2.0, markersize=4.5, zorder=4, label="Test Set F1 gain (274 recs)" if (i == 0 and ax is ax_sw_l) else None)
        if not pl.empty and comp_raw in pl.index:
            rl = pl.loc[comp_raw]
            for ax in (ax_sw_l, ax_sw_r):
                ax.errorbar(rl["diff"], i + 0.16, xerr=[[max(0.0, rl["diff"] - rl["ci_lo"])], [max(0.0, rl["ci_hi"] - rl["diff"])]], fmt="s", color=MID_GREEN, ecolor=MID_GREEN, elinewidth=1.2, capsize=2.0, markersize=4.3, zorder=4, label="Held-out LTAFDB gain (84 recs)" if (i == 0 and ax is ax_sw_l) else None)
        wtl = wtl_map.get(comp_raw, "")
        ax_sw_r.text(max_x + 1.5, i, f"{r['diff']:+.2f} pts  {wtl}", va="center", ha="left", fontsize=7.1, color="#222222")
    ax_sw_l.set_yticks(y_c)
    ylabels = [f"#{r['rank']} {r['comp_clean']}" for _, r in pt.iterrows()]
    ax_sw_l.set_yticklabels(ylabels, fontsize=8.0)
    for tick, (_, r) in zip(ax_sw_l.get_yticklabels(), pt.iterrows()):
        if is_ours(r["comparator"]):
            tick.set_color(MID_GREEN)
            tick.set_fontweight("bold")
    ax_sw_r.text(max_x + 1.5, -0.85, "Test dF1 (records W/T/L)", ha="left", va="bottom", fontsize=7.6, fontweight="bold")
    sig = (pt["verdict"].astype(str) == "reference better")
    if sig.all():
        sw_title = f"B. Fusion Significantly Better Than All {n_comp} Other Detectors (Test Set)"
    else:
        sw_title = f"B. Fusion vs. {n_comp} Other Detectors: Significantly Better Than {int(sig.sum())} (Test Set)"
    ax_sw_l.set_title(sw_title, fontsize=10.2, loc="left", pad=10)
    ax_sw_l.set_xlabel("Paired dF1 = Fusion - Comparator (F1 points, 95% bootstrap CI, linear)", fontsize=8.8, x=0.72)
    ax_sw_l.legend(loc="lower right", bbox_to_anchor=(0.98, 0.16), fontsize=7.5, frameon=True)
    tot = totals_text(dataset_totals(per_record_path))
    fig.suptitle(f"Speed-Accuracy Trade-off and Paired Comparison, {len(sub)} Detectors (Raw +/-150 ms{', ' + tot if tot else ''})", fontsize=12.4, y=0.965)
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(os.path.join(out_dir, f"{stem}.png"), dpi=300)
    fig.savefig(os.path.join(out_dir, f"{stem}.pdf"))
    plt.close(fig)
    print(f"  Saved {stem}.png and {stem}.pdf -> {out_dir}/")


# ===========================================================================
# MAIN
# ===========================================================================
def find_first(results_dir: Path, candidates: list[str]) -> str | None:
    for name in candidates:
        p = results_dir / name
        if p.exists():
            return str(p)
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Generate publication figures from benchmark_three.py outputs.")
    ap.add_argument("--results", default="results/benchmark_v1", help="Folder with benchmark CSVs")
    ap.add_argument("--csv-dir", default=None, help="Alias for --results")
    ap.add_argument("--out", default="figures", help="Output folder for figures")
    ap.add_argument("--out-dir", default=None, help="Alias for --out")
    args = ap.parse_args()

    results_dir = Path(args.csv_dir or args.results)
    out_dir = Path(args.out_dir or args.out)
    print(f"Results dir: {results_dir}  exists={results_dir.exists()}")
    if results_dir.exists():
        print("  Files in results dir:")
        for p in sorted(results_dir.iterdir()):
            print(f"    {p.name}  ({p.stat().st_size} bytes)")
    else:
        print(f"ERROR: results dir does not exist: {results_dir}", file=sys.stderr)
        print("  Did you run benchmark_three.py first? Or is --results pointing at the wrong folder?", file=sys.stderr)
        return 2

    summary_path = find_first(results_dir, ["summary.csv", "headline.csv"])
    rr_path = find_first(results_dir, ["rr_error.csv"])
    per_record_path = find_first(results_dir, ["per_record.csv"])
    paired_test_path = find_first(results_dir, ["paired_raw_150.csv", "paired_test.csv"])
    paired_ltaf_path = find_first(results_dir, ["paired_raw_150_ltafdb.csv", "paired_ltafdb.csv"])
    paired_oracle_path = find_first(results_dir, ["paired_oracle_75.csv"])

    print(f"Using summary: {summary_path}")
    print(f"Using rr: {rr_path}")
    print(f"Using per_record: {per_record_path}")
    print(f"Using paired test: {paired_test_path}")
    print(f"Using paired ltaf: {paired_ltaf_path}")

    if summary_path is None:
        print("ERROR: need summary.csv or headline.csv in --results", file=sys.stderr)
        return 2
    if rr_path is None:
        print("ERROR: need rr_error.csv in --results", file=sys.stderr)
        return 2

    out_dir.mkdir(parents=True, exist_ok=True)

    # Fig 1a: raw 150, the AAMI standard condition
    plot_three_panel_benchmark(summary_path, rr_path, per_record_path, str(out_dir), variant="raw", tol_ms=150.0, stem="fig_comprehensive_raw150", paired_path=paired_test_path)

    # Fig 1b: offset-corrected 75, only if that variant is in summary.csv
    try:
        summ_check = pd.read_csv(summary_path)
        has_offset75 = ((summ_check["variant"] == "offset_corrected") & (summ_check["tol_ms"] == 75.0)).any()
    except Exception as e:
        print(f"  Could not check for offset75 variant: {e}")
        has_offset75 = False
    if has_offset75:
        plot_three_panel_benchmark(summary_path, rr_path, per_record_path, str(out_dir), variant="offset_corrected", tol_ms=75.0, stem="fig_comprehensive_offset75", paired_path=paired_oracle_path)
    else:
        print("\n[three-panel] SKIP offset75: no offset_corrected/75 rows in summary.csv")

    # Fig 2: speed-accuracy + clean sweep, needs paired files
    if paired_test_path is None:
        print("\n[clean-sweep] SKIP: no paired_raw_150.csv in --results, cannot build Fig 2")
    else:
        plot_speed_accuracy_cleansweep(summary_path, paired_test_path, paired_ltaf_path, per_record_path, str(out_dir))

    print(f"\nDone. Check {out_dir}/ for new .png / .pdf timestamps.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())