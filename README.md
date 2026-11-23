# tenhz-ecg

Fast, deterministic QRS detection from a **10 Hz min-max representation** of single-lead ECG (`MinEl`, `TiTen`, and `Fusion`), compiled with Numba.

* **State-of-the-art accuracy across 482 recordings (3,332 hours):** Evaluated across 8 PhysioNet databases (`MITDB`, `SVDB`, `EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`, and held-out `LTAFDB`), outperforming 20 competing detectors.
* **10 Hz runtime efficiency:** **12×–871× faster** than standard Python and WFDB baselines while maintaining sub-5 ms RR interval MAE.
* **Deterministic:** Zero-phase bandpass filtering at the native sampling rate followed by block min-max decimation to 10 Hz, matching the reference implementation beat for beat.

## Installation

```bash
pip install -e .                  # core library (numpy, scipy, numba)
pip install -e ".[data,test]"     # include PhysioNet export & pytest dependencies
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

Evaluated across 274 unseen test recordings (`EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`) and 84 held-out long-term ambulatory recordings (`LTAFDB`, 1,961 hours):

| Detector | Pooled Test F1 (%) | LTAFDB Held-Out F1 (%) | Test RR MAE / RMSE (ms) | Relative Runtime |
| :--- | :---: | :---: | :---: | :---: |
| **Fusion (10 Hz)** | **99.61** | **99.04** | **4.4 / 8.8** | **1.0×** |
| **MinEl (10 Hz)** | **99.54** | **98.97** | **4.3 / 8.9** | **0.73×** |
| **TiTen (10 Hz)** | **99.48** | **98.69** | **4.3 / 9.2** | **0.49×** |
| hamilton2002 (NeuroKit2) | 99.39 | 98.61 | 10.3 / 18.0 | 12× |
| emrich2023 (NeuroKit2) | 99.35 | 97.64 | 4.5 / 9.6 | 20× |
| WFDB `xqrs` | 99.33 | 97.60 | 4.5 / 10.2 | 871× |
| WFDB `gqrs` | 99.04 | 97.37 | 4.7 / 9.9 | 66× |

<p align="center">
  <img src="figures/fig_comprehensive_raw150.png" alt="Comprehensive benchmark, raw 150 ms" width="100%">
</p>

## Reproducing the Benchmarks & Tests

```bash
# Run golden regression tests
pytest

# Download and export PhysioNet databases to data/
python scripts/download_physionet.py --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb

# Run benchmark evaluation
python scripts/benchmark_three.py --data data --dbs edb nsrdb qtdb twadb nstdb --methods all --out results/benchmark_v1
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