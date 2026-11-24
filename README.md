# tenhz-ecg - QRS detection from 10 samples per second

**Accurate ECG QRS detection from a 10 Hz min-max representation of the signal.**

Three single-lead QRS detectors (`MinEl`, `TiTen` and `Fusion`) that operate on 10 values per second, compiled with Numba.

- **Accuracy:** Evaluated on 482 recordings (2,663 h, ~12 M beats) from 8 PhysioNet databases (`MITDB`, `SVDB`, `EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`, `LTAFDB`). `Fusion` scores significantly higher F1 than all 21 competing detectors, on both the test set and the held-out `LTAFDB`, which was not used for tuning. `MinEl` and `TiTen` rank directly below it and are not significantly worse than any competitor.
- **Speed:** `TiTen` is the fastest of the three; `MinEl` takes 1.5× and `Fusion` 1.9× its runtime. Against the next three most accurate detectors (by test F1), `TiTen` runs 36× faster than `hamilton2002` and 72× faster than both WFDB `xqrs` and `emrich2023`, with a mean RR interval error below 5 ms on the test set.
- **Determinism:** Identical detections on a Linux Xeon server and an Apple laptop.

**Method:** The ECG is band-pass filtered at its native sampling rate, and each 200 ms block is reduced to its minimum and maximum, giving 10 values per second. Detection runs on this 10 Hz stream; each detected beat is reported at the native-rate sample of its extremum.

## Installation

```bash
python3 -m venv .venv             # create a virtual environment
source .venv/bin/activate         # activate it
pip install -e .                  # core library (numpy, scipy, numba)
pip install -r requirements.txt   # competing detectors and plotting (pinned versions)
```

## Quickstart

```python
import tenhz

# Optional: compile the Numba kernels before timing
tenhz.warmup()

# Detect R-peaks (detector: "minel", "titen" or "fusion")
res = tenhz.detect(ecg, fs=250, detector="fusion")
print(res.locations)  # 0-based sample indices (int64)
print(res.peaks)      # ECG amplitude at each detected peak

# Sample indices only
locs = tenhz.qrs_locations(ecg, fs=250, detector="titen")
```

### Input requirements and notes
- Single-lead 1-D array in mV, at least **10 seconds** long, finite values only (no `NaN`/`Inf`).
- Filtering is zero-phase (`filtfilt`), so processing is offline.
- The first call compiles the Numba kernels and caches them locally; later runs load them from the cache.
- `detector="fusion"` combines `MinEl` and `TiTen` on a shared 7–32 Hz stream using an RR-plausibility Viterbi decoder and an adaptive per-window switch. Its parameters were tuned on `MITDB` (excluding records `118` and `119`) and `SVDB` only.

## Results

Top 7 of 24 detectors (all are shown in Figure 1). Beats are matched within ±150 ms.

- **Test set:** 274 recordings (`EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`).
- **Held out:** `LTAFDB`, 84 recordings, 1,961 h.
- **RR error:** test set only.
- **Runtime:** relative to `TiTen`, measured on an Intel Xeon Platinum 8358 @ 2.60 GHz.
- **Significant:** the 95% paired bootstrap interval of the F1 difference excludes 0 and the Holm-corrected Wilcoxon p-value over records is below 0.05.

| Detector | Test F1 (%) | LTAFDB F1 (%) | RR error MAE / RMSE (ms) | Runtime (TiTen = 1×) |
| :--- | :---: | :---: | :---: | :---: |
| **Fusion (10 Hz)** | **99.61** | **99.04** | **4.4 / 8.8** | **1.9×** |
| **MinEl (10 Hz)** | **99.54** | **98.97** | **4.3 / 8.9** | **1.5×** |
| **TiTen (10 Hz)** | **99.48** | **98.69** | **4.3 / 9.2** | **1.0×** |
| hamilton2002 | 99.39 | 98.61 | 10.3 / 18.0 | 36× |
| WFDB `xqrs` | 99.37 | 97.77 | 4.5 / 10.1 | 72× |
| emrich2023 | 99.35 | 98.27 | 4.5 / 9.6 | 72× |
| WFDB `gqrs` | 99.04 | 97.37 | 4.7 / 9.9 | 207× |

**Figure 1:** all 24 detectors.
- **A:** F1 on each of the 8 databases and the pooled test set.
- **B:** RR interval error, MAE / RMSE.
- **C:** runtime relative to `TiTen`.

<p align="center">
  <img src="figures/fig_comprehensive_raw150.png" alt="All detectors: F1 per database, RR error and runtime" width="100%">
</p>

**Figure 2:**
- **A:** test F1 against runtime.
- **B:** paired F1 difference between `Fusion` and each other detector, with 95% bootstrap confidence intervals, on the test set (274 recordings) and `LTAFDB` (84 recordings).

<p align="center">
  <img src="figures/fig_speed_accuracy_cleansweep.png" alt="Speed-accuracy trade-off and paired comparison" width="100%">
</p>

### Notes on the competing detectors

All competing detectors run on the same single lead at its native sampling rate. NeuroKit2 methods use the NeuroKit2 cleaning method of the same name where one exists, and `neurokit` cleaning otherwise.

- **`ssf`:** NeuroKit2's default threshold finds no beats in ECG given in mV. The threshold is set to 600 (mV/s)², chosen on the training databases with `scripts/tune_ssf.py`.
- **NeuroKit2 errors:** the installed NeuroKit2 0.3.0 raises on some records for two methods, `manikandan2012` (41 records) and `emrich2023` (1 record, `LTAFDB` 113, which has a flat stretch of signal). Those records are re-run with patched copies of the two functions in `scripts/competitors.py`; every other record uses the installed code unchanged.

## Reproducing the benchmark

**Runtime:** The full benchmark takes about 35 hours of total compute, almost all of it in the slower competing detectors. Run it with several workers unless you need single-process timings.

```bash
# Download the PhysioNet databases into data/
python scripts/download_physionet.py --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb

# Run all detectors on all 8 databases (10 parallel workers)
python scripts/benchmark_three.py --data data --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb --methods all --out results/benchmark_v1 --workers 10

# Build the figures
python scripts/plot_figures.py --results results/benchmark_v1 --out figures/

# Print the runtime table below
python scripts/runtime_table.py results/benchmark_v1
```

Use the pinned versions in `requirements.txt`. Other versions of the competing libraries, WFDB in particular, can change those detectors' results; `MinEl`, `TiTen` and `Fusion` are unaffected. `run_info.json` in the results folder records the library versions, CPU, worker count and wall-clock time of a run.

## Runtime

Intel Xeon Platinum 8358 @ 2.60 GHz, 30 vCPUs, 200 GB RAM. Times are detection only (reading files, saving detections and RR statistics are excluded), summed over all 482 recordings, run with 30 parallel workers.

| Detector | Total time | Share |
| :--- | :---: | :---: |
| `martinez2004` | 866.3 min (14.4 h) | 41.4% |
| `christov2004` | 266.5 min (4.4 h) | 12.7% |
| `zong2003` | 263.4 min (4.4 h) | 12.6% |
| `WFDB gqrs` | 141.9 min (2.4 h) | 6.8% |
| `engzeemod2012` | 115.4 min (1.9 h) | 5.5% |
| `khamis2016` | 106.4 min (1.8 h) | 5.1% |
| 15 other detectors combined | 330.8 min (5.5 h) | 15.8% |
| **Fusion (10 Hz)** | **1.6 min (98 s)** | **0.08%** |
| **MinEl (10 Hz)** | **1.3 min (81 s)** | **0.06%** |
| **TiTen (10 Hz)** | **0.8 min (50 s)** | **0.04%** |
| **All 24 detectors** | **2,094.5 min (34.9 h)** | **100%** |

Total wall-clock time for the full benchmark with 30 workers: **75 min**.

## Citation

If you use `tenhz` in your research, please cite:

```bibtex
@article{jayarathna2026tenhz,
  title   = {MinEl and TiTen: QRS Detection from a 10 Hz Min-Max Representation of the ECG},
  author  = {Jayarathna, Titus},
  year    = {2026}
}
```