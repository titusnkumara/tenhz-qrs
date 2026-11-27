# Paired comparison: offset-corrected (oracle) F1 at +/-75 ms [TELE — held out, nothing tuned on it]

250 held-out TELE records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 93.83)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -1.87 [-3.08, -0.68] | 4.6e-07 | 39 / 112 / 99 | reference worse |
| neurokit | 95.69 | -1.86 [-2.62, -1.05] | 0.00231 | 60 / 101 / 89 | reference worse |
| WFDB gqrs | 94.91 | -1.08 [-1.91, -0.26] | 0.0792 | 70 / 94 / 86 |  |
| Fusion (10 Hz, exploratory) | 94.47 | -0.65 [-1.26, -0.08] | 0.0792 | 54 / 76 / 120 |  |
| TiTen (10 Hz) | 90.22 | +3.60 [+2.32, +5.04] | 1.56e-06 | 111 / 64 / 75 | reference better |
| manikandan2012 | 89.23 | +4.59 [+3.57, +5.75] | 1.04e-15 | 142 / 37 / 71 | reference better |
| elgendi2010 | 88.85 | +4.97 [+3.60, +6.49] | 1.23e-09 | 146 / 56 / 48 | reference better |
| Energy + find_peaks (baseline) | 88.38 | +5.44 [+4.19, +6.97] | 1.81e-16 | 159 / 44 / 47 | reference better |
| rodrigues2021 | 84.51 | +9.31 [+7.68, +10.98] | 9.63e-18 | 153 / 36 / 61 | reference better |
| zong2003 | 78.72 | +15.11 [+12.87, +17.41] | 1.81e-27 | 204 / 25 / 21 | reference better |
| emrich2023 | 73.92 | +19.91 [+16.11, +24.25] | 2.3e-20 | 147 / 35 / 68 | reference better |
| nabian2018 | 73.73 | +20.10 [+14.88, +26.70] | 1.34e-23 | 193 / 33 / 24 | reference better |
| kalidas2017 | 59.01 | +34.82 [+29.40, +40.68] | 3.58e-34 | 207 / 17 / 26 | reference better |
| WFDB xqrs | 58.39 | +35.43 [+29.14, +41.52] | 2.82e-28 | 168 / 22 / 60 | reference better |
| hamilton2002 | 53.80 | +40.02 [+35.02, +45.04] | 8.35e-35 | 213 / 12 / 25 | reference better |
| pantompkins1985 | 53.33 | +40.50 [+35.37, +45.64] | 3.89e-36 | 222 / 8 / 20 | reference better |
| elgendi2013 | 52.64 | +41.18 [+36.37, +46.11] | 7.51e-39 | 232 / 12 / 6 | reference better |
| ssf | 51.93 | +41.89 [+35.27, +48.28] | 2.26e-33 | 198 / 10 / 42 | reference better |
| martinez2004 | 48.24 | +45.58 [+39.34, +51.87] | 2.27e-35 | 211 / 9 / 30 | reference better |
| christov2004 | 44.67 | +49.16 [+46.04, +52.15] | 1.42e-37 | 227 / 14 / 9 | reference better |
| gamboa2008 | 44.29 | +49.54 [+47.08, +52.05] | 1.66e-38 | 230 / 9 / 11 | reference better |
| engzeemod2012 | 41.82 | +52.01 [+47.01, +57.31] | 9.04e-39 | 232 / 4 / 14 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 90.22)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -5.48 [-7.55, -3.49] | 6.74e-12 | 34 / 133 / 83 | reference worse |
| neurokit | 95.69 | -5.47 [-7.06, -3.94] | 3.01e-10 | 51 / 124 / 75 | reference worse |
| WFDB gqrs | 94.91 | -4.69 [-6.26, -3.24] | 2.26e-09 | 56 / 121 / 73 | reference worse |
| Fusion (10 Hz, exploratory) | 94.47 | -4.25 [-5.59, -3.07] | 3.64e-11 | 26 / 112 / 112 | reference worse |
| MinEl (10 Hz) | 93.83 | -3.60 [-5.04, -2.32] | 1.95e-06 | 64 / 111 / 75 | reference worse |
| manikandan2012 | 89.23 | +0.99 [-0.40, +2.31] | 0.113 | 111 / 82 / 57 |  |
| elgendi2010 | 88.85 | +1.37 [-0.55, +3.26] | 0.194 | 117 / 91 / 42 |  |
| Energy + find_peaks (baseline) | 88.38 | +1.84 [+0.03, +3.60] | 0.00318 | 126 / 75 / 49 | reference better |
| rodrigues2021 | 84.51 | +5.71 [+4.05, +7.31] | 0.000225 | 138 / 73 / 39 | reference better |
| zong2003 | 78.72 | +11.50 [+8.68, +14.25] | 1.74e-13 | 181 / 53 / 16 | reference better |
| emrich2023 | 73.92 | +16.31 [+12.55, +20.38] | 2.09e-12 | 135 / 61 / 54 | reference better |
| nabian2018 | 73.73 | +16.50 [+10.99, +23.57] | 1.49e-11 | 181 / 48 / 21 | reference better |
| kalidas2017 | 59.01 | +31.21 [+26.04, +36.85] | 6.31e-34 | 206 / 17 / 27 | reference better |
| WFDB xqrs | 58.39 | +31.83 [+25.67, +37.98] | 1.09e-25 | 160 / 31 / 59 | reference better |
| hamilton2002 | 53.80 | +36.42 [+31.53, +41.16] | 3.9e-35 | 212 / 16 / 22 | reference better |
| pantompkins1985 | 53.33 | +36.89 [+31.76, +41.99] | 9.33e-36 | 217 / 14 / 19 | reference better |
| elgendi2013 | 52.64 | +37.58 [+32.64, +42.66] | 2.24e-34 | 221 / 21 / 8 | reference better |
| ssf | 51.93 | +38.29 [+31.92, +44.73] | 2.28e-31 | 193 / 21 / 36 | reference better |
| martinez2004 | 48.24 | +41.98 [+35.65, +48.26] | 8.1e-34 | 205 / 17 / 28 | reference better |
| christov2004 | 44.67 | +45.56 [+42.44, +48.65] | 1.9e-34 | 219 / 19 / 12 | reference better |
| gamboa2008 | 44.29 | +45.93 [+42.76, +48.96] | 1.26e-33 | 215 / 28 / 7 | reference better |
| engzeemod2012 | 41.82 | +48.41 [+43.53, +53.47] | 3.04e-38 | 231 / 5 / 14 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 94.47)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -1.23 [-2.48, +0.03] | 0.000969 | 44 / 88 / 118 |  |
| neurokit | 95.69 | -1.21 [-2.12, -0.31] | 0.492 | 70 / 82 / 98 |  |
| WFDB gqrs | 94.91 | -0.44 [-1.31, +0.47] | 0.878 | 83 / 77 / 90 |  |
| MinEl (10 Hz) | 93.83 | +0.65 [+0.08, +1.26] | 0.167 | 76 / 54 / 120 |  |
| TiTen (10 Hz) | 90.22 | +4.25 [+3.07, +5.59] | 2.73e-11 | 112 / 26 / 112 | reference better |
| manikandan2012 | 89.23 | +5.24 [+4.09, +6.51] | 8.55e-17 | 147 / 30 / 73 | reference better |
| elgendi2010 | 88.85 | +5.62 [+4.15, +7.22] | 2.14e-10 | 148 / 45 / 57 | reference better |
| Energy + find_peaks (baseline) | 88.38 | +6.09 [+4.72, +7.71] | 8.55e-17 | 161 / 35 / 54 | reference better |
| rodrigues2021 | 84.51 | +9.96 [+8.26, +11.73] | 4.09e-17 | 161 / 28 / 61 | reference better |
| zong2003 | 78.72 | +15.76 [+13.55, +17.94] | 3.85e-26 | 203 / 25 / 22 | reference better |
| emrich2023 | 73.92 | +20.56 [+16.71, +24.88] | 9.85e-23 | 151 / 21 / 78 | reference better |
| nabian2018 | 73.73 | +20.75 [+15.59, +27.19] | 1.22e-22 | 203 / 25 / 22 | reference better |
| kalidas2017 | 59.01 | +35.47 [+29.96, +41.49] | 1.69e-34 | 207 / 10 / 33 | reference better |
| WFDB xqrs | 58.39 | +36.08 [+29.74, +42.28] | 3.67e-29 | 170 / 15 / 65 | reference better |
| hamilton2002 | 53.80 | +40.67 [+35.71, +45.65] | 3.61e-35 | 212 / 8 / 30 | reference better |
| pantompkins1985 | 53.33 | +41.15 [+36.02, +46.25] | 2.75e-36 | 220 / 4 / 26 | reference better |
| elgendi2013 | 52.64 | +41.83 [+37.02, +46.67] | 2.46e-38 | 229 / 8 / 13 | reference better |
| ssf | 51.93 | +42.54 [+36.11, +48.84] | 3.56e-33 | 196 / 9 / 45 | reference better |
| martinez2004 | 48.24 | +46.23 [+40.14, +52.52] | 6.07e-35 | 207 / 7 / 36 | reference better |
| christov2004 | 44.67 | +49.81 [+46.59, +52.93] | 4.02e-36 | 226 / 13 / 11 | reference better |
| gamboa2008 | 44.29 | +50.19 [+47.73, +52.72] | 2.45e-37 | 227 / 12 / 11 | reference better |
| engzeemod2012 | 41.82 | +52.66 [+47.58, +58.04] | 5.58e-38 | 228 / 3 / 19 | reference better |

