# Paired comparison: raw F1 at +/-150 ms [TELE — held out, nothing tuned on it]

250 held-out TELE records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 94.20)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -2.04 [-2.77, -1.26] | 0.000291 | 54 / 102 / 94 | reference worse |
| khamis2016 | 96.05 | -1.85 [-3.04, -0.65] | 2.09e-07 | 37 / 110 / 103 | reference worse |
| WFDB gqrs | 95.30 | -1.10 [-1.93, -0.28] | 0.0303 | 71 / 93 / 86 | reference worse |
| Fusion (10 Hz, exploratory) | 94.85 | -0.65 [-1.24, -0.09] | 0.0465 | 52 / 75 / 123 | reference worse |
| elgendi2010 | 91.83 | +2.38 [+1.32, +3.55] | 5.11e-05 | 129 / 62 / 59 | reference better |
| TiTen (10 Hz) | 90.82 | +3.38 [+2.12, +4.79] | 9.37e-06 | 108 / 65 / 77 | reference better |
| manikandan2012 | 89.96 | +4.24 [+3.24, +5.36] | 8.62e-15 | 139 / 37 / 74 | reference better |
| Energy + find_peaks (baseline) | 88.62 | +5.58 [+4.30, +7.08] | 8.79e-17 | 158 / 45 / 47 | reference better |
| rodrigues2021 | 85.07 | +9.13 [+7.55, +10.77] | 6.22e-17 | 151 / 38 / 61 | reference better |
| zong2003 | 82.82 | +11.38 [+9.25, +13.71] | 7.4e-24 | 192 / 32 / 26 | reference better |
| nabian2018 | 75.34 | +18.86 [+13.68, +25.48] | 6.15e-23 | 191 / 35 / 24 | reference better |
| emrich2023 | 74.93 | +19.27 [+15.42, +23.53] | 5.93e-20 | 145 / 38 / 67 | reference better |
| kalidas2017 | 61.99 | +32.21 [+26.72, +38.17] | 5.33e-33 | 203 / 20 / 27 | reference better |
| pantompkins1985 | 59.98 | +34.23 [+29.49, +39.28] | 7.92e-34 | 213 / 11 / 26 | reference better |
| gamboa2008 | 59.75 | +34.45 [+32.54, +36.46] | 1.29e-36 | 226 / 14 / 10 | reference better |
| hamilton2002 | 59.37 | +34.84 [+30.06, +39.94] | 2.67e-33 | 208 / 15 / 27 | reference better |
| WFDB xqrs | 59.16 | +35.04 [+28.67, +41.23] | 4.87e-28 | 169 / 21 / 60 | reference better |
| elgendi2013 | 53.92 | +40.28 [+35.48, +45.24] | 8.69e-39 | 231 / 13 / 6 | reference better |
| christov2004 | 52.23 | +41.98 [+38.98, +44.82] | 3.02e-35 | 215 / 21 / 14 | reference better |
| ssf | 51.95 | +42.25 [+35.55, +48.61] | 2.48e-33 | 198 / 11 / 41 | reference better |
| martinez2004 | 48.48 | +45.72 [+39.59, +52.04] | 2.61e-35 | 210 / 10 / 30 | reference better |
| engzeemod2012 | 42.93 | +51.28 [+46.27, +56.70] | 1.96e-38 | 229 / 5 / 16 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 90.82)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -5.42 [-6.96, -3.95] | 3.53e-11 | 46 / 123 / 81 | reference worse |
| khamis2016 | 96.05 | -5.23 [-7.23, -3.35] | 1.82e-12 | 32 / 131 / 87 | reference worse |
| WFDB gqrs | 95.30 | -4.48 [-6.00, -3.08] | 2.31e-09 | 54 / 121 / 75 | reference worse |
| Fusion (10 Hz, exploratory) | 94.85 | -4.03 [-5.27, -2.91] | 9.8e-11 | 23 / 113 / 114 | reference worse |
| MinEl (10 Hz) | 94.20 | -3.38 [-4.79, -2.12] | 9.37e-06 | 65 / 108 / 77 | reference worse |
| elgendi2010 | 91.83 | -1.01 [-2.73, +0.69] | 0.394 | 95 / 102 / 53 |  |
| manikandan2012 | 89.96 | +0.86 [-0.54, +2.14] | 0.128 | 109 / 80 / 61 |  |
| Energy + find_peaks (baseline) | 88.62 | +2.19 [+0.49, +3.95] | 0.000566 | 127 / 72 / 51 | reference better |
| rodrigues2021 | 85.07 | +5.75 [+4.12, +7.36] | 0.000233 | 137 / 72 / 41 | reference better |
| zong2003 | 82.82 | +8.00 [+5.27, +10.82] | 1.27e-09 | 171 / 60 / 19 | reference better |
| nabian2018 | 75.34 | +15.48 [+10.16, +22.60] | 5.56e-11 | 179 / 49 / 22 | reference better |
| emrich2023 | 74.93 | +15.89 [+12.08, +20.04] | 8.5e-12 | 133 / 62 / 55 | reference better |
| kalidas2017 | 61.99 | +28.83 [+23.75, +34.38] | 5.74e-32 | 197 / 21 / 32 | reference better |
| pantompkins1985 | 59.98 | +30.84 [+26.23, +35.43] | 8.93e-33 | 204 / 23 / 23 | reference better |
| gamboa2008 | 59.75 | +31.07 [+28.51, +33.51] | 1.38e-27 | 208 / 35 / 7 | reference better |
| hamilton2002 | 59.37 | +31.45 [+26.83, +36.28] | 8.85e-33 | 202 / 23 / 25 | reference better |
| WFDB xqrs | 59.16 | +31.66 [+25.42, +37.84] | 6.5e-25 | 158 / 35 / 57 | reference better |
| elgendi2013 | 53.92 | +36.89 [+31.87, +41.99] | 1.7e-34 | 223 / 18 / 9 | reference better |
| christov2004 | 52.23 | +38.59 [+35.51, +41.64] | 1.83e-31 | 206 / 29 / 15 | reference better |
| ssf | 51.95 | +38.87 [+32.45, +45.31] | 2.26e-31 | 192 / 22 / 36 | reference better |
| martinez2004 | 48.48 | +42.34 [+36.11, +48.60] | 7.56e-34 | 205 / 17 / 28 | reference better |
| engzeemod2012 | 42.93 | +47.89 [+42.89, +53.09] | 1.63e-37 | 226 / 7 / 17 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 94.85)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -1.39 [-2.27, -0.52] | 0.224 | 65 / 79 / 106 |  |
| khamis2016 | 96.05 | -1.20 [-2.39, +0.01] | 0.000521 | 42 / 85 / 123 |  |
| WFDB gqrs | 95.30 | -0.45 [-1.31, +0.46] | 0.753 | 83 / 75 / 92 |  |
| MinEl (10 Hz) | 94.20 | +0.65 [+0.09, +1.24] | 0.139 | 75 / 52 / 123 |  |
| elgendi2010 | 91.83 | +3.03 [+1.81, +4.37] | 4.91e-06 | 129 / 53 / 68 | reference better |
| TiTen (10 Hz) | 90.82 | +4.03 [+2.91, +5.27] | 7.35e-11 | 113 / 23 / 114 | reference better |
| manikandan2012 | 89.96 | +4.89 [+3.86, +6.12] | 3.81e-16 | 145 / 29 / 76 | reference better |
| Energy + find_peaks (baseline) | 88.62 | +6.23 [+4.87, +7.80] | 1.14e-17 | 164 / 32 / 54 | reference better |
| rodrigues2021 | 85.07 | +9.78 [+8.12, +11.56] | 6.43e-17 | 158 / 31 / 61 | reference better |
| zong2003 | 82.82 | +12.03 [+9.91, +14.27] | 6.49e-23 | 194 / 31 / 25 | reference better |
| nabian2018 | 75.34 | +19.51 [+14.47, +26.13] | 7.87e-22 | 200 / 26 / 24 | reference better |
| emrich2023 | 74.93 | +19.92 [+16.08, +24.24] | 1.45e-22 | 150 / 21 / 79 | reference better |
| kalidas2017 | 61.99 | +32.86 [+27.36, +38.92] | 2.65e-33 | 201 / 11 / 38 | reference better |
| pantompkins1985 | 59.98 | +34.88 [+30.08, +39.75] | 1.42e-34 | 209 / 9 / 32 | reference better |
| gamboa2008 | 59.75 | +35.10 [+33.16, +37.21] | 4.9e-34 | 224 / 14 / 12 | reference better |
| hamilton2002 | 59.37 | +35.49 [+30.69, +40.64] | 2.95e-34 | 208 / 10 / 32 | reference better |
| WFDB xqrs | 59.16 | +35.69 [+29.32, +41.87] | 4.76e-29 | 169 / 15 / 66 | reference better |
| elgendi2013 | 53.92 | +40.93 [+36.15, +45.83] | 3.03e-38 | 227 / 9 / 14 | reference better |
| christov2004 | 52.23 | +42.63 [+39.52, +45.67] | 3.55e-34 | 217 / 18 / 15 | reference better |
| ssf | 51.95 | +42.90 [+36.44, +49.30] | 3.74e-33 | 195 / 10 / 45 | reference better |
| martinez2004 | 48.48 | +46.37 [+40.29, +52.72] | 6.59e-35 | 207 / 7 / 36 | reference better |
| engzeemod2012 | 42.93 | +51.93 [+46.82, +57.40] | 3.01e-37 | 225 / 3 / 22 | reference better |

