# tenhz-ecg - QRS detection from 10 samples per second

**Accurate ECG QRS detection from a 10 Hz min-max representation of the signal.**

Three single-lead QRS detectors (`MinEl`, `TiTen` and `Fusion`) that operate on 10 values per second, compiled with Numba.

- **Accuracy:** Evaluated on 482 recordings (2,663 h, ~12 M beats) from 8 PhysioNet databases (`MITDB`, `SVDB`, `EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`, `LTAFDB`). `Fusion` scores significantly higher F1 than all 21 competing detectors, on both the test set and the held-out `LTAFDB`, which was not used for tuning. `MinEl` and `TiTen` rank directly below it and are not significantly worse than any competitor.
- **Speed:** `TiTen` is the fastest of the three; `MinEl` takes 1.7× and `Fusion` 2.3× its runtime. Against the next three most accurate detectors (by test F1), `TiTen` runs 37-53× faster than `hamilton2002`, 71-109× faster than WFDB `xqrs` and 72-107× faster than `emrich2023`, with a mean RR interval error below 5 ms on the test set.
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
- R-peak snap: each detection is moved to the R extremum of its QRS, with one QRS polarity per record decided on the first minute with a clear polarity (at most 10 min scanned, otherwise upright). `fusion` snaps both detectors with `MinEl`'s polarity. Pass `snap=False` for detections without this option.
- Long recordings : run `titen` and `fusion` on multi-day Holter recordings as they are; they recover from large transients and signal dropouts by default. Pass `recover=False` only to reproduce the detections of 2.1.0 and earlier.
- Short recordings and artefact at the start: since 2.3.0 `titen` (and the `TiTen` branch of `fusion`) starts from a quiet 2 s block of the first 10 s, ignores candidates far above its signal level and re-initialises after 5 s without a beat, so a recording that starts with rail saturation is not lost. Pass `startup="2.2"` to reproduce the detections of 2.2.0. See [Changes in 2.3.0](#changes-in-230).
- `detector="fusion"` combines `MinEl` and `TiTen` on a shared 7–32 Hz stream using an RR-plausibility Viterbi decoder and an adaptive per-window switch. Its parameters were tuned on `MITDB` (excluding records `118` and `119`) and `SVDB` only.

## Results

Top 7 of 24 detectors (all are shown in Figure 1). Beats are matched within ±150 ms. Results use the default settings.

**Version:** the accuracy results and figures in this section are for version 2.3.0 (`TiTen` and `Fusion` rerun on 27 November 2026; the competitors and `MinEl` are unchanged from the 2.2.0 run at commit `c83dc21`). 2.3.0 changes only `TiTen`'s start-up; see [Changes in 2.3.0](#changes-in-230). Runtimes, here and in [Runtime](#runtime), are from the 2.2.0 run, measured under the same conditions as the competitors; the 2.3.0 start-up changes detection time by less than 1% (`TiTen` ×1.006, `Fusion` ×1.001, measured back to back).

- **Test set:** 274 recordings (`EDB`, `NSRDB`, `QTDB`, `TWADB`, `NSTDB`).
- **Held out:** `LTAFDB`, 84 recordings, 1,961 h; `TELE` (telehealth ECG, dry handheld electrodes, Khamis et al. 2016), 250 recordings, 2 h, scored outside its visual artefact mask.
- **RR error:** test set only.
- **Runtime:** relative to `TiTen`, measured on an Intel Xeon Platinum 8358 @ 2.60 GHz.
- **Significant:** the 95% paired bootstrap interval of the F1 difference excludes 0 and the Holm-corrected Wilcoxon p-value over records is below 0.05.

| Detector | Test F1 (%) | LTAFDB F1 (%) | TELE F1 (%) | RR error MAE / RMSE (ms) | Runtime (TiTen = 1×) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **Fusion (10 Hz)** | **99.60** | **99.05** | **94.85** | **4.2 / 8.7** | **2.3×** |
| **MinEl (10 Hz)** | **99.52** | **98.95** | **94.20** | **4.2 / 9.5** | **1.7×** |
| **TiTen (10 Hz)** | **99.49** | **98.70** | **90.82** | **4.1 / 8.9** | **1.0×** |
| hamilton2002 | 99.39 | 98.61 | 59.37 | 10.3 / 18.0 | 53× |
| WFDB `xqrs` | 99.37 | 97.77 | 59.16 | 4.5 / 10.1 | 109× |
| emrich2023 | 99.35 | 98.27 | 74.93 | 4.5 / 9.6 | 107× |
| WFDB `gqrs` | 99.04 | 97.37 | 95.30 | 4.7 / 9.9 | 317× |

On `TELE` three competitors outside this top 7 score higher than the 10 Hz detectors: neurokit 96.24, khamis2016 96.05 (the detector written for this database) and WFDB `gqrs` 95.30. The matched filter (py-ecg-detectors) has no `TELE` result.

**Figure 1:** all 24 detectors.
- **A:** F1 on each of the 9 databases and the pooled test set.
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

**Figure 3:** F1 against noise level on `NSTDB` (records 118 and 119 with added electrode-motion noise, 24 to −6 dB SNR), all 24 detectors. The dashed line is the best F1 any detector reaches at each noise level. `Fusion` is within 1.7 F1 points of it at every level and matches it at 24, 12 and −6 dB.

<p align="center">
  <img src="figures/fig_nstdb_snr.png" alt="F1 against noise level on NSTDB" width="100%">
</p>

### Notes on the competing detectors

All competing detectors run on the same single lead at its native sampling rate. NeuroKit2 methods use the NeuroKit2 cleaning method of the same name where one exists, and `neurokit` cleaning otherwise.

- **`ssf`:** NeuroKit2's default threshold finds no beats in ECG given in mV. The threshold is set to 600 (mV/s)², chosen on the training databases with `scripts/tune_ssf.py`.
- **NeuroKit2 errors:** the installed NeuroKit2 0.3.0 raises on some records for two methods, `manikandan2012` (41 records) and `emrich2023` (1 record, `LTAFDB` 113, which has a flat stretch of signal). Those records are re-run with patched copies of the two functions in `scripts/competitors.py`; every other record uses the installed code unchanged.

## Reproducing the benchmark

**Runtime:** The full benchmark takes about 54 hours of total compute, almost all of it in the slower competing detectors. Run it with several workers unless you need single-process timings.

```bash
# Download the PhysioNet databases into data/
python scripts/download_physionet.py --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb

# Run all detectors on all 8 databases (10 parallel workers)
python scripts/benchmark_three.py --data data --dbs mitdb svdb edb nsrdb qtdb twadb nstdb ltafdb --methods all --out results/benchmark_v1 --workers 10

# Build the figures
python scripts/plot_figures.py --results results/benchmark_v1 --out figures/

# Build the NSTDB noise-level figure
python scripts/plot_nstdb_snr.py --results results/benchmark_v1 --out figures/

# Print the runtime table below
python scripts/runtime_table.py results/benchmark_v1
```

Use the pinned versions in `requirements.txt`. Other versions of the competing libraries, WFDB in particular, can change those detectors' results; `MinEl`, `TiTen` and `Fusion` are unaffected. `run_info.json` in the results folder records the library versions, CPU, worker count and wall-clock time of a run.

## Runtime

Version 2.2.0 (commit `c83dc21`). Intel Xeon Platinum 8358 @ 2.60 GHz, 30 vCPUs, 222 GB RAM. Times are detection only (reading files, saving detections and RR statistics are excluded), summed over all 482 recordings, run with 30 parallel workers.

| Detector | Total time | Share |
| :--- | :---: | :---: |
| `martinez2004` | 892.1 min (14.9 h) | 27.5% |
| `christov2004` | 523.0 min (8.7 h) | 16.1% |
| `zong2003` | 501.9 min (8.4 h) | 15.5% |
| `WFDB gqrs` | 295.5 min (4.9 h) | 9.1% |
| `engzeemod2012` | 233.0 min (3.9 h) | 7.2% |
| `khamis2016` | 163.3 min (2.7 h) | 5.0% |
| 15 other detectors combined | 625.4 min (10.4 h) | 19.3% |
| **Fusion (10 Hz)** | **2.6 min (156 s)** | **0.08%** |
| **MinEl (10 Hz)** | **1.8 min (111 s)** | **0.06%** |
| **TiTen (10 Hz)** | **1.0 min (62 s)** | **0.03%** |
| **All 24 detectors** | **3,239.8 min (54.0 h)** | **100%** |

Total wall-clock time for the full benchmark with 30 workers: **114 min**.

## Citation

If you use `tenhz` in your research, please cite:

```bibtex
@article{jayarathna2026tenhz,
  title   = {MinEl and TiTen: QRS Detection from a 10 Hz Min-Max Representation of the ECG},
  author  = {Jayarathna, Titus},
  year    = {2026}
}
```

## Changes in 2.3.0

2.3.0 changes only `TiTen`'s start-up (and with it the `TiTen` branch of `Fusion`). `MinEl` and the `Fusion` decoder are unchanged.

**Why.** On `TELE`, `TiTen` 2.2.0 scored F1 62.49 (sensitivity 46.4%). 132 of the 250 recordings start with several seconds of rail-saturated signal (about 5.5 mV, inside the visual artefact mask). `TiTen` set its levels from the first 2 s, accepted the rail artefact as beats, and its 20 s watchdog, counted from the last accepted artefact, came too late for recordings of 20–40 s. Its median sensitivity was 7.2% on those recordings and 97.4% on the others. `MinEl` keeps no state from the start of a recording and was not affected.

**What changed** (`startup="robust"`, default):
1. The levels start from the 2 s block of the first 10 s with the second-lowest peak energy (order statistic (n−1)//4 of the block maxima), using the 2.2 rule on that block.
2. Candidates above 20× the current signal level are ignored: neither beat nor noise.
3. The watchdog fires after 5 s without a beat (was 20 s) and re-initialises the levels from the last 2 s with the start-up rule, instead of capping the signal level.

No new parameter was tuned. `startup="2.2"` reproduces 2.2.0 exactly.

**How it was chosen.** The design was fixed after the `TELE` failure was found, on training data only. It uses a synthetic set of 496 excerpts of 30 s from `MITDB` (without 118 and 119) and `SVDB`. 374 of them start with rail clipping: duration 1–15 s, rail at 3–10× the QRS amplitude, ranges chosen in advance from amplifier behaviour, not from `TELE`. Four variants were tried in a fixed order; this is the fourth. `TELE` was run once, after the choice. Rule B was relaxed after the result: the chosen variant reached synthetic F1 91.65, which passes the ≥ 90 bar but is 2.37 below `MinEl`, against a pre-set margin of 2.0. The full record is in `results/startup_230/PREREGISTRATION.md`.

F1 (%), beats matched within ±150 ms, default settings.

| Detector | MITDB | SVDB | Synthetic clipped starts | Test F1 | LTAFDB | TELE |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **TiTen 2.3.0** | **99.53** | **99.48** | **91.65** | **99.49** | **98.70** | **90.82** |
| TiTen 2.2.0 | 99.52 | 99.48 | 53.25 | 99.48 | 98.68 | 62.49 |
| **Fusion 2.3.0** | **99.71** | **99.74** | **93.64** | **99.60** | **99.05** | **94.85** |
| Fusion 2.2.0 | 99.71 | 99.74 | 93.19 | 99.60 | 99.05 | 93.66 |
| MinEl (unchanged) | 99.74 | 99.78 | 94.02 | 99.52 | 98.95 | 94.20 |

No test database or `LTAFDB` lost more than 0.1 F1. The largest change is `TiTen` on `TWADB`, −0.054.

## Changes in 2.2.0

2.2.0 changes only the `Fusion` decoder: its local median RR now comes from a first decoding pass instead of from agreed `MinEl` beats only (in stretches where the two detectors disagree that median was about twice the true RR, so true beats were skipped, mostly at high heart rates), and candidates much smaller than the local beat amplitude are penalised. The eight decoder values were tuned on the training databases (`MITDB` without 118 and 119, `SVDB`, clean and with synthetic noise) only. `MinEl` and `TiTen` are unchanged. The 2.1 decoder stays available: `tenhz.fusion.run_chain(ecg, fs, params=tenhz.fusion.FUSION_2_1)`.

F1 (%) with the default settings, beats matched within ±150 ms, same records as above. Fusion 2.1.1 values are from the 2.1.1 benchmark run.

| Detector | Test F1 | LTAFDB F1 | NSTDB F1 | LTAFDB sensitivity ≥ 180 bpm |
| :--- | :---: | :---: | :---: | :---: |
| **Fusion 2.2.0** | **99.60** | **99.05** | **92.19** | **81.5** |
| Fusion 2.1.1 | 99.60 | 99.01 | 91.95 | 76.9 |
| MinEl (unchanged) | 99.52 | 98.95 | 90.78 | 71.7 |
| TiTen (unchanged) | 99.48 | 98.68 | 91.69 | 75.1 |
