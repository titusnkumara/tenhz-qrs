# Paired comparison: offset-corrected (oracle) F1 at +/-75 ms

274 test records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 99.38)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.53 | -0.15 [-0.27, -0.07] | 6.21e-07 | 40 / 103 / 131 | reference worse |
| TiTen (10 Hz) | 99.39 | -0.01 [-0.14, +0.12] | 0.524 | 99 / 69 / 106 |  |
| WFDB xqrs | 99.25 | +0.13 [-0.04, +0.35] | 0.819 | 109 / 92 / 73 |  |
| emrich2023 | 99.22 | +0.17 [-0.01, +0.38] | 0.819 | 98 / 88 / 88 |  |
| WFDB gqrs | 98.94 | +0.44 [+0.14, +0.83] | 1.19e-05 | 151 / 67 / 56 | reference better |
| hamilton2002 | 98.92 | +0.46 [+0.26, +0.70] | 4.04e-21 | 199 / 46 / 29 | reference better |
| kalidas2017 | 98.09 | +1.29 [-0.12, +4.09] | 1.44e-06 | 175 / 71 / 28 |  |
| Matched filter (py-ecg-detectors) | 97.08 | +2.30 [+1.26, +3.46] | 4.24e-29 | 224 / 25 / 25 | reference better |
| manikandan2012 | 96.61 | +2.77 [+1.98, +3.63] | 0.108 | 95 / 76 / 103 |  |
| pantompkins1985 | 96.04 | +3.35 [+1.64, +6.32] | 1.88e-27 | 215 / 36 / 23 | reference better |
| engzeemod2012 | 95.71 | +3.67 [+2.51, +4.90] | 1.57e-37 | 244 / 20 / 10 | reference better |
| neurokit | 94.69 | +4.69 [+3.75, +5.68] | 1.89e-08 | 154 / 73 / 47 | reference better |
| Energy + find_peaks (baseline) | 94.07 | +5.31 [+4.45, +6.24] | 4.25e-14 | 153 / 52 / 69 | reference better |
| khamis2016 | 93.16 | +6.22 [+3.04, +9.78] | 1.04e-08 | 129 / 65 / 80 | reference better |
| elgendi2010 | 92.89 | +6.49 [+5.61, +7.43] | 7.41e-22 | 178 / 38 / 58 | reference better |
| rodrigues2021 | 92.54 | +6.85 [+5.74, +7.99] | 1.39e-06 | 166 / 70 / 38 | reference better |
| ssf | 91.38 | +8.01 [+4.48, +11.77] | 5.14e-28 | 195 / 34 / 45 | reference better |
| nabian2018 | 91.02 | +8.36 [+6.45, +10.43] | 1.21e-37 | 244 / 14 / 16 | reference better |
| elgendi2013 | 90.36 | +9.02 [+5.75, +12.51] | 9.77e-39 | 234 / 11 / 29 | reference better |
| christov2004 | 90.17 | +9.21 [+6.52, +12.34] | 7.61e-28 | 204 / 37 / 33 | reference better |
| zong2003 | 82.91 | +16.47 [+13.70, +19.28] | 5.09e-39 | 244 / 11 / 19 | reference better |
| martinez2004 | 70.07 | +29.32 [+20.25, +40.89] | 1.26e-33 | 204 / 18 / 52 | reference better |
| gamboa2008 | 27.59 | +71.79 [+70.37, +73.15] | 2.62e-45 | 272 / 2 / 0 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 99.39)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.53 | -0.14 [-0.23, -0.07] | 3.47e-15 | 22 / 130 / 122 | reference worse |
| MinEl (10 Hz) | 99.38 | +0.01 [-0.12, +0.14] | 0.699 | 69 / 99 / 106 |  |
| WFDB xqrs | 99.25 | +0.14 [-0.03, +0.34] | 1 | 88 / 104 / 82 |  |
| emrich2023 | 99.22 | +0.17 [+0.02, +0.37] | 1 | 82 / 109 / 83 |  |
| WFDB gqrs | 98.94 | +0.45 [+0.17, +0.81] | 0.0437 | 124 / 88 / 62 | reference better |
| hamilton2002 | 98.92 | +0.47 [+0.24, +0.73] | 2.88e-16 | 189 / 53 / 32 | reference better |
| kalidas2017 | 98.09 | +1.30 [-0.09, +4.09] | 0.000557 | 153 / 79 / 42 |  |
| Matched filter (py-ecg-detectors) | 97.08 | +2.31 [+1.28, +3.46] | 6.25e-26 | 217 / 38 / 19 | reference better |
| manikandan2012 | 96.61 | +2.78 [+1.99, +3.64] | 1 | 85 / 99 / 90 |  |
| pantompkins1985 | 96.04 | +3.35 [+1.66, +6.36] | 6.67e-27 | 213 / 35 / 26 | reference better |
| engzeemod2012 | 95.71 | +3.68 [+2.53, +4.88] | 3.25e-37 | 241 / 22 / 11 | reference better |
| neurokit | 94.69 | +4.70 [+3.77, +5.68] | 7.72e-06 | 134 / 78 / 62 | reference better |
| Energy + find_peaks (baseline) | 94.07 | +5.32 [+4.47, +6.24] | 1.67e-11 | 147 / 58 / 69 | reference better |
| khamis2016 | 93.16 | +6.23 [+3.04, +9.77] | 6.63e-06 | 119 / 79 / 76 | reference better |
| elgendi2010 | 92.89 | +6.50 [+5.63, +7.44] | 1.51e-15 | 161 / 56 / 57 | reference better |
| rodrigues2021 | 92.54 | +6.85 [+5.77, +7.97] | 0.0026 | 147 / 98 / 29 | reference better |
| ssf | 91.38 | +8.01 [+4.50, +11.80] | 6.25e-26 | 186 / 38 / 50 | reference better |
| nabian2018 | 91.02 | +8.37 [+6.47, +10.47] | 1.03e-34 | 235 / 16 / 23 | reference better |
| elgendi2013 | 90.36 | +9.03 [+5.80, +12.54] | 6.51e-36 | 229 / 12 / 33 | reference better |
| christov2004 | 90.17 | +9.22 [+6.56, +12.34] | 3.27e-24 | 192 / 52 / 30 | reference better |
| zong2003 | 82.91 | +16.48 [+13.74, +19.34] | 3.26e-37 | 244 / 20 / 10 | reference better |
| martinez2004 | 70.07 | +29.32 [+20.22, +40.89] | 9.23e-33 | 204 / 18 / 52 | reference better |
| gamboa2008 | 27.59 | +71.80 [+70.38, +73.21] | 2.71e-45 | 271 / 3 / 0 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 99.53)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| TiTen (10 Hz) | 99.39 | +0.14 [+0.07, +0.23] | 2.21e-15 | 130 / 22 / 122 | reference better |
| MinEl (10 Hz) | 99.38 | +0.15 [+0.07, +0.27] | 1.55e-07 | 103 / 40 / 131 | reference better |
| WFDB xqrs | 99.25 | +0.28 [+0.12, +0.50] | 3.25e-07 | 129 / 59 / 86 | reference better |
| emrich2023 | 99.22 | +0.31 [+0.13, +0.54] | 7.06e-11 | 126 / 44 / 104 | reference better |
| WFDB gqrs | 98.94 | +0.59 [+0.31, +0.97] | 9.61e-14 | 162 / 47 / 65 | reference better |
| hamilton2002 | 98.92 | +0.61 [+0.39, +0.87] | 6.7e-30 | 214 / 27 / 33 | reference better |
| kalidas2017 | 98.09 | +1.44 [+0.05, +4.25] | 6.69e-17 | 193 / 42 / 39 | reference better |
| Matched filter (py-ecg-detectors) | 97.08 | +2.45 [+1.40, +3.63] | 1.61e-34 | 235 / 16 / 23 | reference better |
| manikandan2012 | 96.61 | +2.92 [+2.14, +3.78] | 7.06e-11 | 123 / 47 / 104 | reference better |
| pantompkins1985 | 96.04 | +3.49 [+1.78, +6.51] | 5.38e-35 | 230 / 20 / 24 | reference better |
| engzeemod2012 | 95.71 | +3.82 [+2.63, +5.06] | 6.12e-40 | 249 / 15 / 10 | reference better |
| neurokit | 94.69 | +4.84 [+3.91, +5.82] | 2.18e-16 | 165 / 49 / 60 | reference better |
| Energy + find_peaks (baseline) | 94.07 | +5.46 [+4.62, +6.38] | 6.52e-21 | 164 / 33 / 77 | reference better |
| khamis2016 | 93.16 | +6.37 [+3.19, +9.92] | 1.9e-14 | 136 / 50 / 88 | reference better |
| elgendi2010 | 92.89 | +6.64 [+5.77, +7.59] | 1.55e-25 | 180 / 29 / 65 | reference better |
| rodrigues2021 | 92.54 | +6.99 [+5.92, +8.14] | 3.31e-19 | 190 / 47 / 37 | reference better |
| ssf | 91.38 | +8.15 [+4.66, +11.92] | 1.41e-30 | 201 / 23 / 50 | reference better |
| nabian2018 | 91.02 | +8.51 [+6.58, +10.64] | 2.22e-39 | 246 / 8 / 20 | reference better |
| elgendi2013 | 90.36 | +9.17 [+5.90, +12.74] | 4.42e-39 | 235 / 11 / 28 | reference better |
| christov2004 | 90.17 | +9.36 [+6.66, +12.48] | 5.1e-30 | 210 / 32 / 32 | reference better |
| zong2003 | 82.91 | +16.62 [+13.85, +19.45] | 4.87e-40 | 246 / 11 / 17 | reference better |
| martinez2004 | 70.07 | +29.46 [+20.38, +41.04] | 1.19e-34 | 208 / 10 / 56 | reference better |
| gamboa2008 | 27.59 | +71.94 [+70.51, +73.36] | 2.65e-45 | 272 / 2 / 0 | reference better |

