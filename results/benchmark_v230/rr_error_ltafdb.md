# RR-interval error over consecutive matched beats (held out: LTAFDB)

dRR = (detected interval) - (reference interval) for two consecutive reference beats that were both matched (+/-150 ms, one-to-one). A constant offset between a detector's marker and the annotation cancels exactly, so this measures interval timing, which is what rhythm and HRV analysis use; it is not an oracle. It does not cancel a change of offset from one beat to the next (morphology-dependent marker position), and it is conditional on detection: read it next to Se and PPV. Median is read from 1 ms bins.

| Detector | median abs dRR (ms) | mean abs (ms) | RMSE (ms) | >=20 ms (%) | >=40 ms (%) | >=75 ms (%) | median LTAFDB | pairs | Se | PPV |
|---|---|---|---|---|---|---|---|---|---|---|
| MinEl (10 Hz) | 0.5 | 4.7 | 11.4 | 2.4 | 1.6 | 0.5 | 0.5 | 8845665 | 98.98 | 98.92 |
| martinez2004 | 0.5 | 3.3 | 5.3 | 0.1 | 0.0 | 0.0 | 0.5 | 2887883 | 39.21 | 90.56 |
| WFDB gqrs | 0.5 | 4.6 | 12.7 | 4.8 | 1.1 | 0.4 | 0.5 | 8906893 | 99.45 | 95.37 |
| khamis2016 | 0.5 | 5.6 | 12.7 | 4.8 | 1.1 | 0.6 | 0.5 | 8966065 | 99.82 | 74.09 |
| engzeemod2012 | 0.5 | 7.0 | 21.5 | 3.2 | 2.5 | 2.2 | 0.5 | 8479498 | 96.11 | 86.35 |
| WFDB xqrs | 0.5 | 5.9 | 16.2 | 4.7 | 3.4 | 1.0 | 0.5 | 8813425 | 98.36 | 97.19 |
| emrich2023 | 0.5 | 4.2 | 11.0 | 1.5 | 0.7 | 0.5 | 0.5 | 8722116 | 98.39 | 98.15 |
| nabian2018 | 0.5 | 3.8 | 10.0 | 0.8 | 0.5 | 0.4 | 0.5 | 7773792 | 91.09 | 90.64 |
| TiTen (10 Hz) | 0.5 | 3.9 | 8.6 | 1.1 | 0.5 | 0.3 | 0.5 | 8756194 | 98.43 | 98.97 |
| rodrigues2021 | 0.5 | 4.5 | 10.7 | 2.0 | 0.7 | 0.3 | 0.5 | 8716501 | 98.36 | 93.28 |
| Energy + find_peaks (baseline) | 0.5 | 4.2 | 9.9 | 2.7 | 0.6 | 0.2 | 0.5 | 8572317 | 96.82 | 93.81 |
| kalidas2017 | 0.5 | 4.7 | 10.3 | 2.8 | 0.7 | 0.3 | 0.5 | 8650868 | 96.60 | 99.04 |
| ssf | 0.5 | 3.6 | 7.3 | 0.8 | 0.3 | 0.1 | 0.5 | 7600345 | 86.57 | 96.01 |
| elgendi2013 | 0.5 | 5.0 | 12.2 | 4.4 | 1.1 | 0.4 | 0.5 | 8694197 | 98.05 | 81.11 |
| neurokit | 0.5 | 3.8 | 8.4 | 1.2 | 0.5 | 0.2 | 0.5 | 8800626 | 98.56 | 98.94 |
| elgendi2010 | 0.5 | 12.2 | 31.9 | 11.5 | 10.7 | 3.4 | 0.5 | 8687040 | 98.14 | 92.50 |
| Fusion (10 Hz, exploratory) | 0.5 | 4.3 | 10.1 | 1.9 | 1.1 | 0.3 | 0.5 | 8851370 | 99.03 | 99.06 |
| manikandan2012 | 7.5 | 8.6 | 14.1 | 9.8 | 1.2 | 0.3 | 7.5 | 8859983 | 99.06 | 96.83 |
| Matched filter (py-ecg-detectors) | 7.5 | 9.1 | 21.3 | 10.4 | 4.5 | 1.9 | 7.5 | 8785453 | 98.02 | 97.14 |
| hamilton2002 | 7.5 | 10.9 | 19.0 | 23.2 | 3.0 | 0.4 | 7.5 | 8796455 | 98.37 | 98.85 |
| pantompkins1985 | 7.5 | 23.0 | 36.1 | 37.0 | 30.6 | 5.1 | 7.5 | 8690578 | 96.94 | 98.83 |
| christov2004 | 7.5 | 7.1 | 17.5 | 6.1 | 2.5 | 1.4 | 7.5 | 8933464 | 99.61 | 93.73 |
| zong2003 | 7.5 | 13.0 | 31.9 | 12.7 | 8.4 | 5.3 | 7.5 | 8748781 | 98.38 | 78.31 |
| gamboa2008 | 70.5 | 76.3 | 104.5 | 63.3 | 56.2 | 48.1 | 70.5 | 5462413 | 77.25 | 36.25 |
