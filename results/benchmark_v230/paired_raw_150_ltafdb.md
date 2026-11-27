# Paired comparison: raw F1 at +/-150 ms [LTAFDB — held out, nothing tuned on it]

84 held-out LTAFDB records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 98.95)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.05 | -0.10 [-0.17, -0.03] | 2.5e-05 | 19 / 65 / 0 | reference worse |
| neurokit | 98.75 | +0.20 [+0.03, +0.41] | 0.197 | 44 / 40 / 0 |  |
| TiTen (10 Hz) | 98.70 | +0.25 [+0.14, +0.37] | 2.07e-05 | 64 / 20 / 0 | reference better |
| hamilton2002 | 98.61 | +0.34 [+0.06, +0.68] | 0.00174 | 60 / 24 / 0 | reference better |
| emrich2023 | 98.27 | +0.68 [+0.44, +0.99] | 4.87e-10 | 73 / 11 / 0 | reference better |
| manikandan2012 | 97.94 | +1.01 [+0.50, +1.63] | 0.00212 | 50 / 34 / 0 | reference better |
| pantompkins1985 | 97.88 | +1.07 [+0.19, +2.31] | 0.0206 | 55 / 29 / 0 | reference better |
| kalidas2017 | 97.81 | +1.14 [+0.16, +2.48] | 0.155 | 51 / 33 / 0 |  |
| WFDB xqrs | 97.77 | +1.18 [+0.45, +2.12] | 1.15e-05 | 64 / 20 / 0 | reference better |
| Matched filter (py-ecg-detectors) | 97.58 | +1.37 [+0.58, +2.29] | 6.8e-07 | 68 / 16 / 0 | reference better |
| WFDB gqrs | 97.37 | +1.58 [+0.80, +2.62] | 1.63e-05 | 63 / 21 / 0 | reference better |
| christov2004 | 96.58 | +2.36 [+1.15, +3.89] | 2.69e-10 | 75 / 9 / 0 | reference better |
| rodrigues2021 | 95.75 | +3.19 [+2.17, +4.40] | 1.85e-13 | 80 / 4 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.29 | +3.65 [+2.57, +4.85] | 6.08e-13 | 77 / 7 / 0 | reference better |
| elgendi2010 | 95.24 | +3.71 [+2.77, +4.72] | 2.38e-13 | 81 / 3 / 0 | reference better |
| ssf | 91.05 | +7.90 [+4.79, +11.57] | 2.58e-08 | 67 / 17 / 0 | reference better |
| engzeemod2012 | 90.97 | +7.98 [+5.42, +10.81] | 9.49e-10 | 69 / 15 / 0 | reference better |
| nabian2018 | 90.87 | +8.08 [+5.80, +10.99] | 8.08e-14 | 80 / 4 / 0 | reference better |
| elgendi2013 | 88.78 | +10.17 [+8.02, +12.57] | 1.24e-13 | 81 / 3 / 0 | reference better |
| zong2003 | 87.20 | +11.74 [+9.39, +14.18] | 3.93e-14 | 83 / 1 / 0 | reference better |
| khamis2016 | 85.05 | +13.89 [+10.90, +16.81] | 9.84e-14 | 79 / 5 / 0 | reference better |
| martinez2004 | 54.72 | +44.22 [+36.83, +52.12] | 3.93e-14 | 83 / 1 / 0 | reference better |
| gamboa2008 | 49.34 | +49.60 [+47.88, +51.37] | 3.93e-14 | 84 / 0 / 0 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 98.70)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.05 | -0.35 [-0.51, -0.22] | 1.86e-10 | 7 / 77 / 0 | reference worse |
| MinEl (10 Hz) | 98.95 | -0.25 [-0.37, -0.14] | 2.66e-05 | 20 / 64 / 0 | reference worse |
| neurokit | 98.75 | -0.05 [-0.24, +0.18] | 0.609 | 34 / 50 / 0 |  |
| hamilton2002 | 98.61 | +0.09 [-0.13, +0.35] | 1 | 43 / 41 / 0 |  |
| emrich2023 | 98.27 | +0.43 [+0.24, +0.68] | 1.11e-07 | 67 / 17 / 0 | reference better |
| manikandan2012 | 97.94 | +0.76 [+0.30, +1.35] | 1 | 38 / 46 / 0 |  |
| pantompkins1985 | 97.88 | +0.82 [-0.01, +2.05] | 1 | 41 / 43 / 0 |  |
| kalidas2017 | 97.81 | +0.89 [-0.02, +2.19] | 1 | 40 / 44 / 0 |  |
| WFDB xqrs | 97.77 | +0.93 [+0.26, +1.79] | 0.28 | 48 / 36 / 0 |  |
| Matched filter (py-ecg-detectors) | 97.58 | +1.12 [+0.37, +1.99] | 0.0362 | 53 / 31 / 0 | reference better |
| WFDB gqrs | 97.37 | +1.33 [+0.58, +2.35] | 0.0295 | 50 / 34 / 0 | reference better |
| christov2004 | 96.58 | +2.12 [+0.95, +3.59] | 1.67e-05 | 62 / 22 / 0 | reference better |
| rodrigues2021 | 95.75 | +2.95 [+1.94, +4.15] | 6.7e-11 | 75 / 9 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.29 | +3.41 [+2.34, +4.55] | 4.04e-11 | 72 / 12 / 0 | reference better |
| elgendi2010 | 95.24 | +3.46 [+2.52, +4.46] | 2.42e-12 | 77 / 7 / 0 | reference better |
| ssf | 91.05 | +7.65 [+4.57, +11.27] | 8.26e-07 | 63 / 21 / 0 | reference better |
| engzeemod2012 | 90.97 | +7.73 [+5.24, +10.54] | 7.45e-08 | 63 / 21 / 0 | reference better |
| nabian2018 | 90.87 | +7.83 [+5.57, +10.68] | 1.6e-13 | 81 / 3 / 0 | reference better |
| elgendi2013 | 88.78 | +9.92 [+7.73, +12.32] | 6.05e-13 | 80 / 4 / 0 | reference better |
| zong2003 | 87.20 | +11.50 [+9.12, +13.93] | 3.93e-13 | 76 / 8 / 0 | reference better |
| khamis2016 | 85.05 | +13.65 [+10.65, +16.54] | 9.92e-13 | 75 / 9 / 0 | reference better |
| martinez2004 | 54.72 | +43.97 [+36.68, +51.89] | 3.93e-14 | 84 / 0 / 0 | reference better |
| gamboa2008 | 49.34 | +49.35 [+47.56, +51.20] | 3.93e-14 | 84 / 0 / 0 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 99.05)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| MinEl (10 Hz) | 98.95 | +0.10 [+0.03, +0.17] | 1.67e-05 | 65 / 19 / 0 | reference better |
| neurokit | 98.75 | +0.30 [+0.12, +0.51] | 0.000823 | 57 / 27 / 0 | reference better |
| TiTen (10 Hz) | 98.70 | +0.35 [+0.22, +0.51] | 1.59e-10 | 77 / 7 / 0 | reference better |
| hamilton2002 | 98.61 | +0.44 [+0.15, +0.81] | 3.3e-07 | 69 / 15 / 0 | reference better |
| emrich2023 | 98.27 | +0.78 [+0.53, +1.10] | 2.95e-13 | 82 / 2 / 0 | reference better |
| manikandan2012 | 97.94 | +1.11 [+0.60, +1.74] | 6.85e-07 | 65 / 19 / 0 | reference better |
| pantompkins1985 | 97.88 | +1.17 [+0.27, +2.42] | 2.22e-05 | 64 / 20 / 0 | reference better |
| kalidas2017 | 97.81 | +1.24 [+0.24, +2.62] | 0.000251 | 61 / 23 / 0 | reference better |
| WFDB xqrs | 97.77 | +1.28 [+0.53, +2.25] | 3.15e-06 | 63 / 21 / 0 | reference better |
| Matched filter (py-ecg-detectors) | 97.58 | +1.47 [+0.64, +2.41] | 4.28e-08 | 68 / 16 / 0 | reference better |
| WFDB gqrs | 97.37 | +1.68 [+0.89, +2.72] | 3.3e-07 | 67 / 17 / 0 | reference better |
| christov2004 | 96.58 | +2.46 [+1.23, +4.04] | 4.96e-11 | 76 / 8 / 0 | reference better |
| rodrigues2021 | 95.75 | +3.29 [+2.28, +4.51] | 4.3e-14 | 83 / 1 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.29 | +3.76 [+2.66, +4.96] | 2.3e-13 | 79 / 5 / 0 | reference better |
| elgendi2010 | 95.24 | +3.81 [+2.86, +4.81] | 6.87e-14 | 82 / 2 / 0 | reference better |
| ssf | 91.05 | +8.00 [+4.87, +11.74] | 1.38e-09 | 68 / 16 / 0 | reference better |
| engzeemod2012 | 90.97 | +8.08 [+5.50, +10.92] | 2.39e-10 | 70 / 14 / 0 | reference better |
| nabian2018 | 90.87 | +8.18 [+5.88, +11.09] | 4.91e-14 | 81 / 3 / 0 | reference better |
| elgendi2013 | 88.78 | +10.27 [+8.10, +12.67] | 4.49e-14 | 81 / 3 / 0 | reference better |
| zong2003 | 87.20 | +11.85 [+9.49, +14.30] | 4.3e-14 | 82 / 2 / 0 | reference better |
| khamis2016 | 85.05 | +13.99 [+10.93, +16.96] | 2.56e-13 | 76 / 8 / 0 | reference better |
| martinez2004 | 54.72 | +44.32 [+36.96, +52.22] | 3.93e-14 | 84 / 0 / 0 | reference better |
| gamboa2008 | 49.34 | +49.70 [+47.97, +51.49] | 3.93e-14 | 84 / 0 / 0 | reference better |

