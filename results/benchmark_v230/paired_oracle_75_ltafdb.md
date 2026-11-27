# Paired comparison: offset-corrected (oracle) F1 at +/-75 ms [LTAFDB — held out, nothing tuned on it]

84 held-out LTAFDB records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 98.73)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 98.86 | -0.13 [-0.22, -0.06] | 9.36e-07 | 15 / 69 / 0 | reference worse |
| neurokit | 98.61 | +0.12 [-0.04, +0.32] | 0.459 | 43 / 41 / 0 |  |
| TiTen (10 Hz) | 98.50 | +0.23 [+0.12, +0.34] | 0.000185 | 61 / 23 / 0 | reference better |
| hamilton2002 | 98.27 | +0.46 [+0.05, +1.05] | 0.0105 | 58 / 26 / 0 | reference better |
| emrich2023 | 97.94 | +0.79 [+0.47, +1.25] | 2.43e-09 | 71 / 13 / 0 | reference better |
| manikandan2012 | 97.70 | +1.03 [+0.48, +1.71] | 0.00823 | 49 / 35 / 0 | reference better |
| kalidas2017 | 97.55 | +1.18 [+0.13, +2.58] | 0.383 | 51 / 33 / 0 |  |
| WFDB gqrs | 97.07 | +1.66 [+0.82, +2.74] | 2.25e-05 | 63 / 21 / 0 | reference better |
| WFDB xqrs | 96.99 | +1.74 [+0.58, +3.45] | 3.5e-06 | 62 / 22 / 0 | reference better |
| Matched filter (py-ecg-detectors) | 95.84 | +2.89 [+1.18, +5.03] | 1.14e-07 | 70 / 14 / 0 | reference better |
| christov2004 | 95.63 | +3.10 [+1.48, +5.17] | 1.81e-09 | 74 / 10 / 0 | reference better |
| rodrigues2021 | 95.50 | +3.23 [+2.17, +4.51] | 1.49e-11 | 74 / 10 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.11 | +3.62 [+2.52, +4.84] | 1.15e-11 | 75 / 9 / 0 | reference better |
| pantompkins1985 | 94.63 | +4.10 [+2.94, +5.56] | 8.86e-12 | 73 / 11 / 0 | reference better |
| elgendi2010 | 93.51 | +5.22 [+3.91, +6.67] | 2.06e-13 | 80 / 4 / 0 | reference better |
| ssf | 90.97 | +7.76 [+4.68, +11.37] | 7.48e-08 | 65 / 19 / 0 | reference better |
| nabian2018 | 90.29 | +8.44 [+6.00, +11.53] | 8.68e-14 | 80 / 4 / 0 | reference better |
| elgendi2013 | 88.52 | +10.21 [+8.02, +12.62] | 2.87e-13 | 80 / 4 / 0 | reference better |
| engzeemod2012 | 88.37 | +10.36 [+6.92, +14.29] | 6.19e-10 | 69 / 15 / 0 | reference better |
| khamis2016 | 84.80 | +13.92 [+10.89, +16.86] | 5.61e-13 | 77 / 7 / 0 | reference better |
| zong2003 | 83.68 | +15.04 [+12.33, +17.89] | 4.79e-14 | 83 / 1 / 0 | reference better |
| martinez2004 | 54.69 | +44.04 [+36.71, +51.98] | 3.93e-14 | 83 / 1 / 0 | reference better |
| gamboa2008 | 31.72 | +67.01 [+65.47, +68.54] | 3.93e-14 | 84 / 0 / 0 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 98.50)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 98.86 | -0.36 [-0.49, -0.25] | 2.58e-11 | 6 / 78 / 0 | reference worse |
| MinEl (10 Hz) | 98.73 | -0.23 [-0.34, -0.12] | 0.000259 | 23 / 61 / 0 | reference worse |
| neurokit | 98.61 | -0.10 [-0.30, +0.12] | 0.841 | 35 / 49 / 0 |  |
| hamilton2002 | 98.27 | +0.24 [-0.15, +0.82] | 0.844 | 47 / 37 / 0 |  |
| emrich2023 | 97.94 | +0.57 [+0.29, +0.95] | 5.66e-08 | 68 / 16 / 0 | reference better |
| manikandan2012 | 97.70 | +0.80 [+0.27, +1.47] | 0.841 | 39 / 45 / 0 |  |
| kalidas2017 | 97.55 | +0.96 [-0.09, +2.34] | 0.844 | 39 / 45 / 0 |  |
| WFDB gqrs | 97.07 | +1.43 [+0.60, +2.55] | 0.0398 | 50 / 34 / 0 | reference better |
| WFDB xqrs | 96.99 | +1.51 [+0.38, +3.16] | 0.0181 | 53 / 31 / 0 | reference better |
| Matched filter (py-ecg-detectors) | 95.84 | +2.67 [+0.99, +4.78] | 0.000134 | 62 / 22 / 0 | reference better |
| christov2004 | 95.63 | +2.88 [+1.26, +4.91] | 2.89e-06 | 66 / 18 / 0 | reference better |
| rodrigues2021 | 95.50 | +3.01 [+1.96, +4.28] | 3.71e-10 | 74 / 10 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.11 | +3.39 [+2.31, +4.58] | 3.13e-10 | 71 / 13 / 0 | reference better |
| pantompkins1985 | 94.63 | +3.88 [+2.74, +5.31] | 9.76e-11 | 73 / 11 / 0 | reference better |
| elgendi2010 | 93.51 | +4.99 [+3.70, +6.41] | 3.86e-13 | 78 / 6 / 0 | reference better |
| ssf | 90.97 | +7.53 [+4.49, +11.13] | 1.01e-06 | 63 / 21 / 0 | reference better |
| nabian2018 | 90.29 | +8.22 [+5.83, +11.27] | 1.53e-13 | 81 / 3 / 0 | reference better |
| elgendi2013 | 88.52 | +9.99 [+7.78, +12.39] | 8.36e-13 | 80 / 4 / 0 | reference better |
| engzeemod2012 | 88.37 | +10.14 [+6.75, +14.00] | 3.32e-08 | 64 / 20 / 0 | reference better |
| khamis2016 | 84.80 | +13.70 [+10.65, +16.63] | 3.49e-12 | 75 / 9 / 0 | reference better |
| zong2003 | 83.68 | +14.82 [+12.10, +17.67] | 1.21e-13 | 82 / 2 / 0 | reference better |
| martinez2004 | 54.69 | +43.81 [+36.53, +51.68] | 3.93e-14 | 84 / 0 / 0 | reference better |
| gamboa2008 | 31.72 | +66.79 [+65.23, +68.34] | 3.93e-14 | 84 / 0 / 0 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 98.86)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| MinEl (10 Hz) | 98.73 | +0.13 [+0.06, +0.22] | 4.68e-07 | 69 / 15 / 0 | reference better |
| neurokit | 98.61 | +0.26 [+0.08, +0.47] | 0.00172 | 57 / 27 / 0 | reference better |
| TiTen (10 Hz) | 98.50 | +0.36 [+0.25, +0.49] | 1.94e-11 | 78 / 6 / 0 | reference better |
| hamilton2002 | 98.27 | +0.60 [+0.16, +1.25] | 1.01e-07 | 70 / 14 / 0 | reference better |
| emrich2023 | 97.94 | +0.93 [+0.59, +1.40] | 5.27e-14 | 83 / 1 / 0 | reference better |
| manikandan2012 | 97.70 | +1.16 [+0.61, +1.85] | 4.68e-07 | 68 / 16 / 0 | reference better |
| kalidas2017 | 97.55 | +1.32 [+0.21, +2.76] | 0.000177 | 64 / 20 / 0 | reference better |
| WFDB gqrs | 97.07 | +1.79 [+0.93, +2.89] | 1.86e-07 | 67 / 17 / 0 | reference better |
| WFDB xqrs | 96.99 | +1.87 [+0.69, +3.66] | 1.31e-07 | 65 / 19 / 0 | reference better |
| Matched filter (py-ecg-detectors) | 95.84 | +3.02 [+1.27, +5.24] | 6.32e-10 | 74 / 10 / 0 | reference better |
| christov2004 | 95.63 | +3.24 [+1.57, +5.36] | 5.23e-11 | 76 / 8 / 0 | reference better |
| rodrigues2021 | 95.50 | +3.37 [+2.31, +4.69] | 1.86e-13 | 82 / 2 / 0 | reference better |
| Energy + find_peaks (baseline) | 95.11 | +3.75 [+2.64, +4.98] | 6.25e-13 | 79 / 5 / 0 | reference better |
| pantompkins1985 | 94.63 | +4.24 [+3.08, +5.72] | 4.94e-13 | 79 / 5 / 0 | reference better |
| elgendi2010 | 93.51 | +5.35 [+4.04, +6.80] | 5.77e-14 | 82 / 2 / 0 | reference better |
| ssf | 90.97 | +7.89 [+4.80, +11.54] | 2.64e-09 | 68 / 16 / 0 | reference better |
| nabian2018 | 90.29 | +8.58 [+6.09, +11.72] | 5.77e-14 | 82 / 2 / 0 | reference better |
| elgendi2013 | 88.52 | +10.35 [+8.17, +12.76] | 8.81e-14 | 81 / 3 / 0 | reference better |
| engzeemod2012 | 88.37 | +10.50 [+7.03, +14.45] | 1.29e-10 | 71 / 13 / 0 | reference better |
| khamis2016 | 84.80 | +14.06 [+10.98, +17.03] | 3.38e-13 | 77 / 7 / 0 | reference better |
| zong2003 | 83.68 | +15.18 [+12.44, +18.05] | 4.97e-14 | 82 / 2 / 0 | reference better |
| martinez2004 | 54.69 | +44.17 [+36.86, +52.09] | 3.93e-14 | 84 / 0 / 0 | reference better |
| gamboa2008 | 31.72 | +67.15 [+65.61, +68.67] | 3.93e-14 | 84 / 0 / 0 | reference better |

