"""Print the README runtime table (markdown) from a benchmark folder.

    python scripts/runtime_table.py results/benchmark_final

Detection time of each detector, summed over every record (raw +/-150 ms rows, one per record), with the
share of the total. Wall-clock time, worker count and CPU come from run_info.json when it is in the folder.
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import competitors as bt  # noqa: E402  (OURS, LABEL)


def fmt(sec):
    m = sec / 60
    return f"{m:,.1f} min ({m / 60:.1f} h)" if m >= 10 else f"{m:.1f} min ({sec:.0f} s)"


def main():
    folder = sys.argv[1] if len(sys.argv) > 1 else sys.exit(__doc__)
    R = pd.read_csv(os.path.join(folder, "per_record.csv"), dtype={"record": str}, keep_default_na=False, low_memory=False)
    R = R[(R.variant == "raw") & (R.tol_ms == R.tol_ms.max())]
    s = pd.to_numeric(R.seconds).groupby(R.method).sum().sort_values(ascending=False)
    tot = s.sum()
    ours = [m for m in ("fusion", "minel", "titen") if m in s.index]
    rest = [m for m in s.index if m not in ours]
    show, other = rest[:6], rest[6:]
    name = lambda m: f"`{bt.LABEL.get(m, m).replace(' (10 Hz, exploratory)', ' (10 Hz)')}`" if m not in bt.OURS else \
        f"**{bt.LABEL[m].replace(' (10 Hz, exploratory)', ' (10 Hz)')}**"
    n_rec = R.drop_duplicates(["db", "record"]).shape[0]
    lines = ["| Detector | Total time | Share |", "| :--- | :---: | :---: |"]
    for m in show:
        lines.append(f"| {name(m)} | {fmt(s[m])} | {100 * s[m] / tot:.1f}% |")
    lines.append(f"| {len(other)} other detectors combined | {fmt(s[other].sum())} | {100 * s[other].sum() / tot:.1f}% |")
    for m in ours:
        lines.append(f"| {name(m)} | **{fmt(s[m])}** | **{100 * s[m] / tot:.2f}%** |")
    lines.append(f"| **All {len(s)} detectors** | **{fmt(tot)}** | **100%** |")
    print("\n".join(lines))
    ri = os.path.join(folder, "run_info.json")
    print(f"\n({n_rec} records; times are detection only, summed over records.)")
    if os.path.exists(ri):
        i = json.load(open(ri))
        print(f"Run: {i['workers']} workers, {i['wall_clock_min']:.1f} min wall-clock, {i['cpu']}, "
              f"{i['logical_cpus']} logical CPUs. Ideal (total / workers): {tot / 60 / i['workers']:.1f} min.")


if __name__ == "__main__":
    main()
