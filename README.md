# tenhz-ecg

Fast, deterministic QRS detection from a **10 Hz min-max representation** of single-lead ECG (`MinEl`, `TiTen`, and `Fusion`), compiled with Numba.

* **State-of-the-art accuracy across 482 recordings (3,332 hours):** 8 PhysioNet databases (`MITDB`, `SVDB`, `EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`, held-out `LTAFDB`), outperforming 20 competing detectors. 124 tuning, 274 unseen test, 84 held-out.
* **10 Hz runtime efficiency:** `TiTen` is the 1.0x reference. `Fusion` is 1.8x `TiTen`, `MinEl` is 1.5x `TiTen`. `TiTen` is 24x faster than `hamilton2002` and 1775x faster than WFDB `xqrs`, while maintaining sub-5 ms RR interval MAE on the test set.
* **Deterministic:** Zero-phase bandpass filtering at the native sampling rate followed by block min-max decimation to 10 Hz.

## Installation

```bash
python3 -m venv .venv             # create virtual environment
source .venv/bin/activate         # activate 
pip install -e .                  # core library (numpy, scipy, numba)
pip install -e ".[data,test]"     # include PhysioNet export & pytest dependencies
pip install -r requirements.txt
```

## Quickstart

```python
import tenhz

# Optional: pre-compile Numba kernels before timing
tenhz.warmup()

# Detect R-peaks (options: "minel", "titen", "fusion")
res = tenhz.detect(ecg, fs=250, detector="fusion")
print(res.locations)  # 0-based sample indices (int64)
print(res.peaks)      # ECG amplitudes at detected peaks

# Sample indices only
locs = tenhz.qrs_locations(ecg, fs=250, detector="titen")
```

### Input Requirements & Notes
* Single-lead 1-D array (`mV`), minimum duration **10 seconds**, finite values only (no `NaN`/`Inf`).
* Offline zero-phase filtering (`filtfilt`). The first call compiles Numba kernels and caches them locally; subsequent runs load directly from cache.
* `detector="fusion"` combines `MinEl` and `TiTen` on a shared 7–32 Hz stream via an RR-plausibility Viterbi decoder and adaptive window switch, tuned strictly on `MITDB` (excluding records `118` and `119`) and `SVDB`.

## Benchmark Summary (Unseen Test Set, Raw ±150 ms AAMI Standard)

Top 7. All 23 in Figure 1. Test: 274 recordings (`EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`). Held-out: 84 recordings (`LTAFDB`, 1,961 hours). RR is test set only. Runtime: `TiTen` = 1.0x.

| Detector | Pooled Test F1 (%) | LTAFDB Held-Out F1 (%) | Test RR MAE / RMSE (ms) | Relative Runtime (TiTen = 1.0x) |
| :--- | :---: | :---: | :---: | :---: |
| **Fusion (10 Hz)** | **99.61** | **99.02** | **4.4 / 8.8** | **1.8x** |
| **MinEl (10 Hz)** | **99.54** | **98.97** | **4.3 / 8.9** | **1.5x** |
| **TiTen (10 Hz)** | **99.48** | **98.69** | **4.3 / 9.2** | **1.0x** |
| hamilton2002 | 99.39 | 98.61 | 10.3 / 18.0 | 23.5x |
| emrich2023 | 99.35 | 97.64 | 4.5 / 9.6 | 41.5x |
| WFDB `xqrs` | 99.33 | 97.60 | 4.5 / 10.2 | 1775x |
| WFDB `gqrs` | 99.04 | 97.37 | 4.7 / 9.9 | 134x |

Figure 1 — All 23 detectors, raw 150 ms. A: F1 across all 8 DBs plus pooled test. B: test RR timing error, MAE / RMSE in ms. C: relative runtime, 1.0x is TiTen.

<p align="center">
  <img src="figures/fig_comprehensive_raw150.png" alt="Comprehensive benchmark, raw 150 ms" width="100%">
</p>

Figure 2 — Speed-accuracy trade-off and paired clean sweep, raw 150 ms. A: pooled test F1 vs relative runtime. B: paired per-record F1 gain, Fusion minus each comparator, with 95% bootstrap CI, on test (274 recs) and held-out LTAFDB (84 recs).

<p align="center">
  <img src="figures/fig_speed_accuracy_cleansweep.png" alt="Speed-accuracy trade-off and paired clean sweep, raw 150 ms" width="100%">
</p>

## Reproducing the Benchmarks & Tests

```bash

# Download and export PhysioNet databases to data/
python scripts/download_physionet.py --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb

# Benchmark all 8 DBs
python scripts/benchmark_three.py --data data --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb --methods all --out results/benchmark_v1

# Build figures
python scripts/plot_figures.py --results results/benchmark_v1 --out figures/
```

## Citation

If you use `tenhz` in your research, please cite:

```bibtex
@article{jayarathna2026tenhz,
  title   = {MinEl and TiTen: QRS Detection from a 10 Hz Min-Max Representation of the ECG},
  author  = {Jayarathna, Titus},
  year    = {2026}
}
```
