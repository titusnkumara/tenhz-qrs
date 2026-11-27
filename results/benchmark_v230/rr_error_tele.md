# RR-interval error over consecutive matched beats (held out: TELE)

dRR = (detected interval) - (reference interval) for two consecutive reference beats that were both matched (+/-150 ms, one-to-one). A constant offset between a detector's marker and the annotation cancels exactly, so this measures interval timing, which is what rhythm and HRV analysis use; it is not an oracle. It does not cancel a change of offset from one beat to the next (morphology-dependent marker position), and it is conditional on detection: read it next to Se and PPV. Median is read from 1 ms bins.

| Detector | median abs dRR (ms) | mean abs (ms) | RMSE (ms) | >=20 ms (%) | >=40 ms (%) | >=75 ms (%) | median TELE | pairs | Se | PPV |
|---|---|---|---|---|---|---|---|---|---|---|
| khamis2016 | 2.5 | 4.8 | 12.3 | 3.9 | 1.4 | 0.6 | 2.5 | 6165 | 97.06 | 95.06 |
| ssf | 2.5 | 4.4 | 7.7 | 3.6 | 0.5 | 0.1 | 2.5 | 1986 | 35.30 | 98.34 |
| martinez2004 | 2.5 | 6.1 | 14.1 | 5.3 | 2.3 | 0.9 | 2.5 | 1923 | 33.77 | 85.93 |
| MinEl (10 Hz) | 4.5 | 7.5 | 14.6 | 8.5 | 3.1 | 0.7 | 4.5 | 5817 | 92.90 | 95.54 |
| emrich2023 | 4.5 | 9.8 | 20.6 | 11.3 | 4.9 | 2.1 | 4.5 | 3880 | 63.60 | 91.19 |
| WFDB gqrs | 4.5 | 9.9 | 19.0 | 15.8 | 5.0 | 1.0 | 4.5 | 5888 | 94.34 | 96.29 |
| WFDB xqrs | 4.5 | 9.1 | 20.5 | 11.8 | 4.8 | 2.2 | 4.5 | 2764 | 43.78 | 91.18 |
| Energy + find_peaks (baseline) | 4.5 | 7.8 | 15.3 | 10.6 | 3.3 | 0.6 | 4.5 | 5025 | 81.77 | 96.74 |
| nabian2018 | 4.5 | 8.5 | 22.3 | 7.8 | 3.8 | 2.3 | 4.5 | 4052 | 72.33 | 78.61 |
| rodrigues2021 | 4.5 | 8.6 | 17.9 | 11.0 | 4.1 | 1.1 | 4.5 | 4844 | 86.52 | 83.67 |
| TiTen (10 Hz) | 4.5 | 7.3 | 15.4 | 7.4 | 2.9 | 1.0 | 4.5 | 5239 | 86.20 | 95.97 |
| engzeemod2012 | 4.5 | 9.9 | 24.6 | 9.7 | 5.6 | 3.5 | 4.5 | 1495 | 28.82 | 84.12 |
| elgendi2013 | 4.5 | 10.7 | 23.5 | 15.4 | 5.5 | 2.8 | 4.5 | 2394 | 43.93 | 69.80 |
| neurokit | 4.5 | 9.1 | 18.1 | 10.9 | 5.4 | 1.3 | 4.5 | 6053 | 96.17 | 96.31 |
| Fusion (10 Hz, exploratory) | 4.5 | 7.3 | 14.0 | 8.3 | 2.8 | 0.6 | 4.5 | 5840 | 93.28 | 96.48 |
| kalidas2017 | 6.5 | 17.0 | 38.4 | 16.6 | 10.3 | 7.4 | 6.5 | 3029 | 48.91 | 84.63 |
| elgendi2010 | 6.5 | 21.2 | 42.4 | 25.0 | 19.5 | 6.4 | 6.5 | 5875 | 94.72 | 89.10 |
| zong2003 | 6.5 | 19.4 | 39.1 | 23.7 | 14.7 | 8.3 | 6.5 | 5580 | 89.06 | 77.40 |
| manikandan2012 | 6.5 | 12.5 | 21.8 | 18.4 | 6.2 | 1.5 | 6.5 | 5373 | 86.69 | 93.49 |
| hamilton2002 | 24.5 | 35.0 | 55.8 | 55.0 | 26.6 | 13.0 | 24.5 | 3020 | 49.14 | 74.98 |
| gamboa2008 | 26.5 | 59.8 | 89.4 | 53.1 | 46.0 | 36.5 | 26.5 | 4717 | 84.44 | 46.24 |
| pantompkins1985 | 30.5 | 41.4 | 61.5 | 60.5 | 38.5 | 18.4 | 30.5 | 3185 | 50.57 | 73.69 |
| christov2004 | 32.5 | 48.9 | 68.4 | 59.7 | 46.2 | 28.5 | 32.5 | 5438 | 85.67 | 37.56 |
