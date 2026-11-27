# Paired comparison: offset-corrected (oracle) F1 at +/-75 ms [TELE — held out, nothing tuned on it]

250 held-out TELE records. Difference = reference minus comparator, in F1 points, pooled over records; CI from 2000 paired bootstrap resamples of records (the same resample for every detector, stratified by database); p from a Wilcoxon signed-rank test on per-record F1, Holm-corrected over the comparators of that reference. The two tests answer different questions: the CI is on pooled F1 (every beat counts equally, so long records weigh more), the Wilcoxon test on the typical record (every record counts equally). A comparison is marked only when both agree: CI excludes 0 and p(Holm) < 0.05.

## MinEl (10 Hz) against each other detector (reference F1 93.83)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -1.87 [-3.08, -0.68] | 3.68e-07 | 39 / 112 / 99 | reference worse |
| neurokit | 95.69 | -1.86 [-2.62, -1.05] | 0.00231 | 60 / 101 / 89 | reference worse |
| WFDB gqrs | 94.91 | -1.08 [-1.91, -0.26] | 0.0396 | 70 / 94 / 86 | reference worse |
| Fusion (10 Hz, exploratory) | 93.28 | +0.54 [+0.07, +1.04] | 0.0141 | 65 / 44 / 141 | reference better |
| manikandan2012 | 89.23 | +4.59 [+3.57, +5.75] | 8.88e-16 | 142 / 37 / 71 | reference better |
| elgendi2010 | 88.85 | +4.97 [+3.60, +6.49] | 1.02e-09 | 146 / 56 / 48 | reference better |
| Energy + find_peaks (baseline) | 88.38 | +5.44 [+4.19, +6.97] | 1.58e-16 | 159 / 44 / 47 | reference better |
| rodrigues2021 | 84.51 | +9.31 [+7.68, +10.98] | 8.56e-18 | 153 / 36 / 61 | reference better |
| zong2003 | 78.72 | +15.11 [+12.87, +17.41] | 1.66e-27 | 204 / 25 / 21 | reference better |
| emrich2023 | 73.92 | +19.91 [+16.11, +24.25] | 2.07e-20 | 147 / 35 / 68 | reference better |
| nabian2018 | 73.73 | +20.10 [+14.88, +26.70] | 1.22e-23 | 193 / 33 / 24 | reference better |
| TiTen (10 Hz) | 62.11 | +31.72 [+26.13, +37.76] | 6.76e-30 | 178 / 25 / 47 | reference better |
| kalidas2017 | 59.01 | +34.82 [+29.40, +40.68] | 3.58e-34 | 207 / 17 / 26 | reference better |
| WFDB xqrs | 58.39 | +35.43 [+29.14, +41.52] | 2.6e-28 | 168 / 22 / 60 | reference better |
| hamilton2002 | 53.80 | +40.02 [+35.02, +45.04] | 8.35e-35 | 213 / 12 / 25 | reference better |
| pantompkins1985 | 53.33 | +40.50 [+35.37, +45.64] | 3.89e-36 | 222 / 8 / 20 | reference better |
| elgendi2013 | 52.64 | +41.18 [+36.37, +46.11] | 7.51e-39 | 232 / 12 / 6 | reference better |
| ssf | 51.93 | +41.89 [+35.27, +48.28] | 2.26e-33 | 198 / 10 / 42 | reference better |
| martinez2004 | 48.24 | +45.58 [+39.34, +51.87] | 2.27e-35 | 211 / 9 / 30 | reference better |
| christov2004 | 44.67 | +49.16 [+46.04, +52.15] | 1.42e-37 | 227 / 14 / 9 | reference better |
| gamboa2008 | 44.29 | +49.54 [+47.08, +52.05] | 1.66e-38 | 230 / 9 / 11 | reference better |
| engzeemod2012 | 41.82 | +52.01 [+47.01, +57.31] | 9.04e-39 | 232 / 4 / 14 | reference better |

## TiTen (10 Hz) against each other detector (reference F1 62.11)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -33.59 [-40.12, -27.18] | 4.83e-29 | 25 / 187 / 38 | reference worse |
| neurokit | 95.69 | -33.58 [-39.81, -27.61] | 1.68e-29 | 30 / 184 / 36 | reference worse |
| WFDB gqrs | 94.91 | -32.80 [-39.09, -26.82] | 3.18e-28 | 43 / 177 / 30 | reference worse |
| MinEl (10 Hz) | 93.83 | -31.72 [-37.76, -26.13] | 1.09e-29 | 25 / 178 / 47 | reference worse |
| Fusion (10 Hz, exploratory) | 93.28 | -31.17 [-37.06, -25.86] | 9.4e-30 | 8 / 180 / 62 | reference worse |
| manikandan2012 | 89.23 | -27.12 [-32.75, -21.87] | 7.98e-26 | 42 / 170 / 38 | reference worse |
| elgendi2010 | 88.85 | -26.74 [-33.20, -21.00] | 3.08e-23 | 60 / 168 / 22 | reference worse |
| Energy + find_peaks (baseline) | 88.38 | -26.27 [-32.11, -20.66] | 4.49e-23 | 50 / 166 / 34 | reference worse |
| rodrigues2021 | 84.51 | -22.41 [-28.22, -17.31] | 1.52e-21 | 71 / 165 / 14 | reference worse |
| zong2003 | 78.72 | -16.61 [-23.80, -9.56] | 4.5e-17 | 86 / 155 / 9 | reference worse |
| emrich2023 | 73.92 | -11.81 [-16.79, -7.26] | 5.43e-05 | 77 / 104 / 69 | reference worse |
| nabian2018 | 73.73 | -11.62 [-21.32, -0.37] | 1.61e-10 | 100 / 141 / 9 | reference worse |
| kalidas2017 | 59.01 | +3.10 [-0.17, +6.14] | 3.47e-07 | 141 / 37 / 72 |  |
| WFDB xqrs | 58.39 | +3.71 [-2.71, +11.20] | 0.00213 | 109 / 47 / 94 |  |
| hamilton2002 | 53.80 | +8.31 [+2.71, +14.22] | 5.27e-08 | 146 / 38 / 66 | reference better |
| pantompkins1985 | 53.33 | +8.78 [+2.81, +15.38] | 1.38e-09 | 153 / 38 / 59 | reference better |
| elgendi2013 | 52.64 | +9.46 [+2.81, +16.85] | 0.00223 | 145 / 74 / 31 | reference better |
| ssf | 51.93 | +10.18 [+0.78, +19.43] | 0.0653 | 121 / 82 / 47 |  |
| martinez2004 | 48.24 | +13.86 [+6.43, +21.89] | 1.09e-07 | 136 / 41 / 73 | reference better |
| christov2004 | 44.67 | +17.44 [+13.42, +20.90] | 1 | 133 / 95 / 22 |  |
| gamboa2008 | 44.29 | +17.82 [+10.28, +24.93] | 1 | 112 / 133 / 5 |  |
| engzeemod2012 | 41.82 | +20.29 [+16.07, +24.38] | 1.23e-19 | 167 / 20 / 63 | reference better |

## Fusion (10 Hz, exploratory) against each other detector (reference F1 93.28)

| Comparator | Its F1 | Difference [95% paired CI] | p (Holm) | records: reference better / worse / tied | |
|---|---|---|---|---|---|
| khamis2016 | 95.70 | -2.42 [-3.85, -0.98] | 2.54e-08 | 40 / 116 / 94 | reference worse |
| neurokit | 95.69 | -2.41 [-3.41, -1.39] | 0.000153 | 61 / 111 / 78 | reference worse |
| WFDB gqrs | 94.91 | -1.63 [-2.63, -0.63] | 0.00136 | 70 / 102 / 78 | reference worse |
| MinEl (10 Hz) | 93.83 | -0.54 [-1.04, -0.07] | 0.00707 | 44 / 65 / 141 | reference worse |
| manikandan2012 | 89.23 | +4.05 [+3.02, +5.17] | 6.84e-11 | 140 / 43 / 67 | reference better |
| elgendi2010 | 88.85 | +4.43 [+2.98, +6.04] | 8.12e-06 | 133 / 70 / 47 | reference better |
| Energy + find_peaks (baseline) | 88.38 | +4.90 [+3.61, +6.38] | 1.21e-11 | 142 / 53 / 55 | reference better |
| rodrigues2021 | 84.51 | +8.77 [+7.22, +10.37] | 1.55e-13 | 154 / 41 / 55 | reference better |
| zong2003 | 78.72 | +14.56 [+12.21, +17.05] | 5.88e-23 | 197 / 32 / 21 | reference better |
| emrich2023 | 73.92 | +19.37 [+15.71, +23.51] | 4.29e-19 | 145 / 37 / 68 | reference better |
| nabian2018 | 73.73 | +19.55 [+14.19, +26.35] | 7.33e-19 | 189 / 40 / 21 | reference better |
| TiTen (10 Hz) | 62.11 | +31.17 [+25.86, +37.06] | 5.56e-30 | 180 / 8 / 62 | reference better |
| kalidas2017 | 59.01 | +34.27 [+29.06, +40.13] | 3.94e-34 | 204 / 16 / 30 | reference better |
| WFDB xqrs | 58.39 | +34.89 [+28.61, +41.08] | 9.81e-29 | 169 / 18 / 63 | reference better |
| hamilton2002 | 53.80 | +39.48 [+34.46, +44.44] | 1.69e-34 | 210 / 12 / 28 | reference better |
| pantompkins1985 | 53.33 | +39.95 [+34.71, +44.99] | 7.79e-36 | 219 / 8 / 23 | reference better |
| elgendi2013 | 52.64 | +40.64 [+35.74, +45.62] | 7.82e-38 | 224 / 16 / 10 | reference better |
| ssf | 51.93 | +41.35 [+34.73, +47.76] | 5.64e-33 | 195 / 12 / 43 | reference better |
| martinez2004 | 48.24 | +45.04 [+38.83, +51.29] | 1.28e-34 | 205 / 12 / 33 | reference better |
| christov2004 | 44.67 | +48.62 [+45.63, +51.52] | 4.53e-36 | 224 / 13 / 13 | reference better |
| gamboa2008 | 44.29 | +48.99 [+46.37, +51.58] | 1.19e-36 | 225 / 15 / 10 | reference better |
| engzeemod2012 | 41.82 | +51.46 [+46.51, +56.77] | 7.29e-38 | 229 / 3 / 18 | reference better |

