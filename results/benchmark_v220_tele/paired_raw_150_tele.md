# Paired comparison: raw F1 at +/-150 ms [TELE — held out, nothing tuned on it]

250 held-out TELE records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 94.20)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -2.04 [-2.77, -1.26] | 0.000291 | 54 / 102 / 94 | reference worse |
| khamis2016 | 96.05 | -1.85 [-3.04, -0.65] | 1.74e-07 | 37 / 110 / 103 | reference worse |
| WFDB gqrs | 95.30 | -1.10 [-1.93, -0.28] | 0.0175 | 71 / 93 / 86 | reference worse |
| Fusion (10 Hz, exploratory) | 93.66 | +0.54 [+0.08, +1.04] | 0.0175 | 63 / 43 / 144 | reference better |
| elgendi2010 | 91.83 | +2.38 [+1.32, +3.55] | 5.11e-05 | 129 / 62 / 59 | reference better |
| manikandan2012 | 89.96 | +4.24 [+3.24, +5.36] | 7.39e-15 | 139 / 37 / 74 | reference better |
| Energy + find_peaks (baseline) | 88.62 | +5.58 [+4.30, +7.08] | 7.69e-17 | 158 / 45 / 47 | reference better |
| rodrigues2021 | 85.07 | +9.13 [+7.55, +10.77] | 5.53e-17 | 151 / 38 / 61 | reference better |
| zong2003 | 82.82 | +11.38 [+9.25, +13.71] | 6.79e-24 | 192 / 32 / 26 | reference better |
| nabian2018 | 75.34 | +18.86 [+13.68, +25.48] | 5.6e-23 | 191 / 35 / 24 | reference better |
| emrich2023 | 74.93 | +19.27 [+15.42, +23.53] | 5.34e-20 | 145 / 38 / 67 | reference better |
| TiTen (10 Hz) | 62.49 | +31.71 [+26.17, +37.73] | 8.63e-30 | 177 / 24 / 49 | reference better |
| kalidas2017 | 61.99 | +32.21 [+26.72, +38.17] | 5.33e-33 | 203 / 20 / 27 | reference better |
| pantompkins1985 | 59.98 | +34.23 [+29.49, +39.28] | 7.92e-34 | 213 / 11 / 26 | reference better |
| gamboa2008 | 59.75 | +34.45 [+32.54, +36.46] | 1.29e-36 | 226 / 14 / 10 | reference better |
| hamilton2002 | 59.37 | +34.84 [+30.06, +39.94] | 2.67e-33 | 208 / 15 / 27 | reference better |
| WFDB xqrs | 59.16 | +35.04 [+28.67, +41.23] | 4.5e-28 | 169 / 21 / 60 | reference better |
| elgendi2013 | 53.92 | +40.28 [+35.48, +45.24] | 8.69e-39 | 231 / 13 / 6 | reference better |
| christov2004 | 52.23 | +41.98 [+38.98, +44.82] | 3.02e-35 | 215 / 21 / 14 | reference better |
| ssf | 51.95 | +42.25 [+35.55, +48.61] | 2.48e-33 | 198 / 11 / 41 | reference better |
| martinez2004 | 48.48 | +45.72 [+39.59, +52.04] | 2.61e-35 | 210 / 10 / 30 | reference better |
| engzeemod2012 | 42.93 | +51.28 [+46.27, +56.70] | 1.96e-38 | 229 / 5 / 16 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 62.49)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -33.75 [-39.99, -27.80] | 2.39e-29 | 28 / 182 / 40 | reference worse |
| khamis2016 | 96.05 | -33.56 [-40.12, -27.16] | 6.89e-29 | 25 / 184 / 41 | reference worse |
| WFDB gqrs | 95.30 | -32.81 [-39.10, -26.82] | 4.6e-28 | 43 / 176 / 31 | reference worse |
| MinEl (10 Hz) | 94.20 | -31.71 [-37.73, -26.17] | 1.39e-29 | 24 / 177 / 49 | reference worse |
| Fusion (10 Hz, exploratory) | 93.66 | -31.17 [-37.01, -25.91] | 1.21e-29 | 7 / 179 / 64 | reference worse |
| elgendi2010 | 91.83 | -29.34 [-35.59, -23.55] | 6.77e-25 | 53 / 171 / 26 | reference worse |
| manikandan2012 | 89.96 | -27.47 [-33.13, -22.29] | 6.53e-26 | 42 / 169 / 39 | reference worse |
| Energy + find_peaks (baseline) | 88.62 | -26.14 [-31.98, -20.52] | 1.12e-22 | 51 / 164 / 35 | reference worse |
| rodrigues2021 | 85.07 | -22.58 [-28.44, -17.48] | 1.11e-21 | 71 / 164 / 15 | reference worse |
| zong2003 | 82.82 | -20.33 [-27.43, -13.23] | 6.26e-19 | 82 / 157 / 11 | reference worse |
| nabian2018 | 75.34 | -12.85 [-22.73, -1.48] | 4.44e-12 | 97 / 145 / 8 | reference worse |
| emrich2023 | 74.93 | -12.44 [-17.56, -7.95] | 1.71e-05 | 74 / 107 / 69 | reference worse |
| kalidas2017 | 61.99 | +0.50 [-2.80, +3.46] | 1.71e-05 | 131 / 44 / 75 |  |
| pantompkins1985 | 59.98 | +2.51 [-2.47, +7.48] | 1.16e-06 | 143 / 45 / 62 |  |
| gamboa2008 | 59.75 | +2.74 [-4.27, +9.30] | 0.00027 | 102 / 143 / 5 |  |
| hamilton2002 | 59.37 | +3.12 [-1.79, +8.21] | 3.51e-06 | 140 / 42 / 68 |  |
| WFDB xqrs | 59.16 | +3.33 [-3.14, +10.89] | 0.00327 | 105 / 53 / 92 |  |
| elgendi2013 | 53.92 | +8.56 [+1.88, +15.91] | 0.00508 | 145 / 73 / 32 | reference better |
| christov2004 | 52.23 | +10.26 [+6.11, +13.87] | 0.0365 | 107 / 119 / 24 | reference better |
| ssf | 51.95 | +10.54 [+1.09, +19.82] | 0.0365 | 120 / 83 / 47 | reference better |
| martinez2004 | 48.48 | +14.01 [+6.65, +21.99] | 1.29e-07 | 135 / 42 / 73 | reference better |
| engzeemod2012 | 42.93 | +19.56 [+15.46, +23.56] | 9.53e-19 | 165 / 20 / 65 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 93.66)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| neurokit | 96.24 | -2.58 [-3.57, -1.59] | 2.83e-05 | 58 / 107 / 85 | reference worse |
| khamis2016 | 96.05 | -2.39 [-3.77, -0.99] | 1.17e-08 | 40 / 113 / 97 | reference worse |
| WFDB gqrs | 95.30 | -1.64 [-2.63, -0.59] | 0.00268 | 70 / 100 / 80 | reference worse |
| MinEl (10 Hz) | 94.20 | -0.54 [-1.04, -0.08] | 0.00877 | 43 / 63 / 144 | reference worse |
| elgendi2010 | 91.83 | +1.84 [+0.57, +3.22] | 0.00719 | 118 / 74 / 58 | reference better |
| manikandan2012 | 89.96 | +3.70 [+2.72, +4.79] | 4.01e-10 | 137 / 43 / 70 | reference better |
| Energy + find_peaks (baseline) | 88.62 | +5.04 [+3.74, +6.51] | 2.41e-12 | 145 / 50 / 55 | reference better |
| rodrigues2021 | 85.07 | +8.59 [+7.07, +10.21] | 4.4e-13 | 151 / 44 / 55 | reference better |
| zong2003 | 82.82 | +10.84 [+8.58, +13.30] | 4.83e-19 | 186 / 40 / 24 | reference better |
| nabian2018 | 75.34 | +18.32 [+13.04, +25.21] | 1.62e-17 | 186 / 43 / 21 | reference better |
| emrich2023 | 74.93 | +18.73 [+15.01, +22.93] | 9.44e-19 | 144 / 37 / 69 | reference better |
| TiTen (10 Hz) | 62.49 | +31.17 [+25.91, +37.01] | 7.16e-30 | 179 / 7 / 64 | reference better |
| kalidas2017 | 61.99 | +31.67 [+26.33, +37.59] | 9.57e-33 | 198 / 18 / 34 | reference better |
| pantompkins1985 | 59.98 | +33.69 [+28.87, +38.64] | 7.23e-34 | 212 / 11 / 27 | reference better |
| gamboa2008 | 59.75 | +33.91 [+31.88, +36.00] | 2.2e-32 | 219 / 20 / 11 | reference better |
| hamilton2002 | 59.37 | +34.30 [+29.45, +39.34] | 4.17e-33 | 206 / 15 / 29 | reference better |
| WFDB xqrs | 59.16 | +34.50 [+28.15, +40.76] | 1.04e-28 | 169 / 18 / 63 | reference better |
| elgendi2013 | 53.92 | +39.74 [+34.86, +44.80] | 1.04e-37 | 224 / 16 / 10 | reference better |
| christov2004 | 52.23 | +41.44 [+38.55, +44.21] | 1.51e-33 | 214 / 19 / 17 | reference better |
| ssf | 51.95 | +41.71 [+35.08, +48.12] | 6.93e-33 | 194 / 13 / 43 | reference better |
| martinez2004 | 48.48 | +45.18 [+39.00, +51.44] | 1.15e-34 | 205 / 12 / 33 | reference better |
| engzeemod2012 | 42.93 | +50.74 [+45.76, +56.08] | 2.79e-37 | 225 / 5 / 20 | reference better |

