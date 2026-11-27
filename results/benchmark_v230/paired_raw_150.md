# Paired comparison: raw F1 at +/-150 ms

274 test records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 99.52)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.60 | -0.08 [-0.12, -0.04] | 0.000575 | 47 / 91 / 136 | reference worse |
| TiTen (10 Hz) | 99.49 | +0.04 [-0.02, +0.11] | 0.00419 | 110 / 56 / 108 |  |
| hamilton2002 | 99.39 | +0.14 [-0.00, +0.31] | 5.34e-14 | 179 / 57 / 38 |  |
| WFDB xqrs | 99.37 | +0.15 [+0.01, +0.36] | 0.262 | 110 / 87 / 77 |  |
| emrich2023 | 99.35 | +0.17 [+0.03, +0.36] | 0.118 | 105 / 76 / 93 |  |
| WFDB gqrs | 99.04 | +0.49 [+0.22, +0.86] | 4.23e-09 | 161 / 58 / 55 | reference better |
| kalidas2017 | 98.18 | +1.35 [-0.04, +4.13] | 0.000378 | 147 / 76 / 51 |  |
| pantompkins1985 | 97.86 | +1.66 [+0.22, +4.48] | 5.8e-11 | 163 / 61 / 50 | reference better |
| manikandan2012 | 96.68 | +2.84 [+2.06, +3.71] | 0.0142 | 90 / 69 / 115 | reference better |
| engzeemod2012 | 96.15 | +3.38 [+2.28, +4.51] | 5e-38 | 243 / 20 / 11 | reference better |
| neurokit | 94.82 | +4.70 [+3.78, +5.66] | 2.11e-13 | 165 / 59 / 50 | reference better |
| Energy + find_peaks (baseline) | 94.13 | +5.39 [+4.54, +6.31] | 3.46e-16 | 154 / 46 / 74 | reference better |
| elgendi2010 | 93.64 | +5.88 [+4.97, +6.85] | 1.05e-22 | 180 / 29 / 65 | reference better |
| khamis2016 | 93.22 | +6.31 [+3.12, +9.83] | 6.4e-14 | 138 / 53 / 83 | reference better |
| Matched filter (py-ecg-detectors) | 93.17 | +6.35 [+4.64, +8.22] | 3.07e-32 | 207 / 27 / 40 | reference better |
| rodrigues2021 | 92.59 | +6.93 [+5.84, +8.09] | 2.62e-13 | 179 / 54 / 41 | reference better |
| christov2004 | 91.44 | +8.09 [+5.63, +10.99] | 1.03e-25 | 188 / 43 / 43 | reference better |
| ssf | 91.39 | +8.13 [+4.63, +11.91] | 1.11e-30 | 200 / 26 / 48 | reference better |
| nabian2018 | 91.14 | +8.38 [+6.45, +10.50] | 1.31e-37 | 242 / 16 / 16 | reference better |
| elgendi2013 | 90.44 | +9.08 [+5.84, +12.63] | 1.02e-39 | 238 / 9 / 27 | reference better |
| zong2003 | 83.75 | +15.77 [+13.06, +18.56] | 1.26e-36 | 230 / 18 / 26 | reference better |
| martinez2004 | 70.09 | +29.43 [+20.35, +41.02] | 1.15e-34 | 209 / 14 / 51 | reference better |
| gamboa2008 | 50.42 | +49.10 [+46.82, +51.55] | 2.9e-45 | 272 / 2 / 0 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 99.49)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| Fusion (10 Hz, exploratory) | 99.60 | -0.11 [-0.19, -0.06] | 1.73e-16 | 20 / 129 / 125 | reference worse |
| MinEl (10 Hz) | 99.52 | -0.04 [-0.11, +0.02] | 0.00628 | 56 / 110 / 108 |  |
| hamilton2002 | 99.39 | +0.10 [-0.04, +0.28] | 2.21e-07 | 161 / 74 / 39 |  |
| WFDB xqrs | 99.37 | +0.11 [-0.03, +0.32] | 0.516 | 81 / 108 / 85 |  |
| emrich2023 | 99.35 | +0.14 [+0.02, +0.30] | 0.433 | 77 / 112 / 85 |  |
| WFDB gqrs | 99.04 | +0.45 [+0.19, +0.81] | 0.0113 | 121 / 87 / 66 | reference better |
| kalidas2017 | 98.18 | +1.31 [-0.08, +4.09] | 0.656 | 111 / 92 / 71 |  |
| pantompkins1985 | 97.86 | +1.63 [+0.18, +4.43] | 6.55e-05 | 146 / 79 / 49 | reference better |
| manikandan2012 | 96.68 | +2.81 [+2.02, +3.67] | 0.802 | 79 / 99 / 96 |  |
| engzeemod2012 | 96.15 | +3.34 [+2.25, +4.47] | 7.51e-38 | 241 / 21 / 12 | reference better |
| neurokit | 94.82 | +4.66 [+3.75, +5.64] | 2.68e-05 | 132 / 78 / 64 | reference better |
| Energy + find_peaks (baseline) | 94.13 | +5.35 [+4.50, +6.26] | 2.89e-12 | 152 / 52 / 70 | reference better |
| elgendi2010 | 93.64 | +5.85 [+4.94, +6.81] | 1.99e-11 | 147 / 61 / 66 | reference better |
| khamis2016 | 93.22 | +6.27 [+3.07, +9.80] | 1.58e-07 | 124 / 77 / 73 | reference better |
| Matched filter (py-ecg-detectors) | 93.17 | +6.32 [+4.62, +8.20] | 1.44e-24 | 193 / 42 / 39 | reference better |
| rodrigues2021 | 92.59 | +6.89 [+5.81, +8.02] | 4.26e-05 | 154 / 91 / 29 | reference better |
| christov2004 | 91.44 | +8.05 [+5.59, +10.96] | 9.73e-19 | 168 / 64 / 42 | reference better |
| ssf | 91.39 | +8.09 [+4.58, +11.87] | 1.38e-27 | 189 / 33 / 52 | reference better |
| nabian2018 | 91.14 | +8.34 [+6.43, +10.44] | 3.01e-35 | 236 / 15 / 23 | reference better |
| elgendi2013 | 90.44 | +9.05 [+5.83, +12.57] | 6.58e-38 | 236 / 9 / 29 | reference better |
| zong2003 | 83.75 | +15.74 [+13.00, +18.52] | 5.23e-32 | 221 / 38 / 15 | reference better |
| martinez2004 | 70.09 | +29.39 [+20.31, +41.00] | 9.73e-34 | 209 / 16 / 49 | reference better |
| gamboa2008 | 50.42 | +49.07 [+46.80, +51.49] | 3.27e-45 | 270 / 4 / 0 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 99.60)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| MinEl (10 Hz) | 99.52 | +0.08 [+0.04, +0.12] | 0.000115 | 91 / 47 / 136 | reference better |
| TiTen (10 Hz) | 99.49 | +0.11 [+0.06, +0.19] | 9.87e-17 | 129 / 20 / 125 | reference better |
| hamilton2002 | 99.39 | +0.21 [+0.08, +0.38] | 3.8e-23 | 197 / 34 / 43 | reference better |
| WFDB xqrs | 99.37 | +0.23 [+0.08, +0.45] | 2.68e-06 | 129 / 58 / 87 | reference better |
| emrich2023 | 99.35 | +0.25 [+0.09, +0.43] | 1.26e-08 | 121 / 47 / 106 | reference better |
| WFDB gqrs | 99.04 | +0.56 [+0.29, +0.93] | 7.28e-16 | 165 / 42 / 67 | reference better |
| kalidas2017 | 98.18 | +1.42 [+0.03, +4.24] | 4.59e-10 | 157 / 54 / 63 | reference better |
| pantompkins1985 | 97.86 | +1.74 [+0.28, +4.57] | 4.83e-20 | 178 / 42 / 54 | reference better |
| manikandan2012 | 96.68 | +2.92 [+2.13, +3.77] | 2.12e-10 | 118 / 44 / 112 | reference better |
| engzeemod2012 | 96.15 | +3.45 [+2.35, +4.61] | 9.58e-41 | 249 / 14 / 11 | reference better |
| neurokit | 94.82 | +4.78 [+3.86, +5.74] | 3.94e-17 | 166 / 46 / 62 | reference better |
| Energy + find_peaks (baseline) | 94.13 | +5.47 [+4.62, +6.39] | 3e-20 | 163 / 31 / 80 | reference better |
| elgendi2010 | 93.64 | +5.96 [+5.04, +6.92] | 3.51e-23 | 175 / 27 / 72 | reference better |
| khamis2016 | 93.22 | +6.38 [+3.19, +9.93] | 1.08e-16 | 145 / 44 / 85 | reference better |
| Matched filter (py-ecg-detectors) | 93.17 | +6.43 [+4.71, +8.31] | 3.94e-30 | 204 / 27 / 43 | reference better |
| rodrigues2021 | 92.59 | +7.00 [+5.93, +8.15] | 8.01e-24 | 199 / 40 / 35 | reference better |
| christov2004 | 91.44 | +8.16 [+5.68, +11.05] | 2.2e-25 | 189 / 44 / 41 | reference better |
| ssf | 91.39 | +8.21 [+4.71, +11.97] | 6.56e-32 | 204 / 19 / 51 | reference better |
| nabian2018 | 91.14 | +8.46 [+6.53, +10.58] | 1.78e-39 | 246 / 8 / 20 | reference better |
| elgendi2013 | 90.44 | +9.16 [+5.89, +12.71] | 8.02e-40 | 239 / 8 / 27 | reference better |
| zong2003 | 83.75 | +15.85 [+13.13, +18.65] | 1.88e-36 | 229 / 22 / 23 | reference better |
| martinez2004 | 70.09 | +29.51 [+20.40, +41.10] | 1.87e-35 | 213 / 9 / 52 | reference better |
| gamboa2008 | 50.42 | +49.18 [+46.88, +51.61] | 2.93e-45 | 271 / 3 / 0 | reference better |

