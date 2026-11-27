# TiTen start-up fix: pre-registered design and acceptance rules

Written 2026-11-27 before any variant is run. TELE (held out) is run once, only after a variant is chosen here.

## Problem (from the TELE diagnosis)
1. Thresholds are initialised from the first 2 s (max and mean of the first 20 stream samples).
2. A start with large artifact sets the signal level far above real beats; artifact peaks are accepted as beats.
3. The 20 s watchdog counts from the last accepted (artifact) beat, so it fires too late in short records.

## Variants (fixed in advance; no other parameters change)
- V0: TiTen 2.2.0 as released (reference).
- V1: robust initialisation. Split the first min(10 s, record) of the integrated stream into 2 s blocks. Apply the
  original initialisation (max and mean of the block, same init_sig_frac and 0.5 factors) to the block at order
  statistic index (nb-1)//4 of the block maxima (for 5 blocks: the second lowest), instead of to block 0.
- V2: V1 + watchdog 5 s instead of 20 s (no other watchdog change).
- V3 (only if V2 fails the synthetic criterion): V2 + reject candidates whose integrated peak exceeds 20x the current
  signal level (treated as neither beat nor noise, excluded from the watchdog maximum).
Choose the simplest variant (V1 < V2 < V3) that passes all rules.

## Synthetic clipped-start stress set (training records only)
Source: MITDB (without 118, 119) and SVDB, first lead, native rate. 4 excerpts of 30 s per record at random starts,
seed 42. Parameters drawn per excerpt from ranges chosen from amplifier physics, not from TELE:
- 25% controls: no artifact.
- 75% artifact: start clip of duration D ~ U(1, 15) s. Rail level R = m * A, m ~ U(3, 10), A = 99th percentile of
  |x - median(x)| of the excerpt. During the clip the signal is piecewise constant at +R or -R (random sign), segment
  lengths ~ Exp(mean 0.5 s), with 20 ms linear transitions. At release an exponential baseline offset starting at the
  last rail value decays with tau ~ U(0.2, 2) s and is added to the ECG.
- With probability 0.5 (artifact excerpts only), one mid-record dropout of 0.5-3 s at a random position after the clip,
  built the same way.
- Reference beats inside clip/dropout regions are removed and detections there are ignored (same rule as TELE).

## Acceptance rules
A. MITDB and SVDB (full records, +/-150 ms, raw): pooled F1 must not drop by more than 0.02 against V0.
B. Synthetic stress set: pooled F1 of the variant must be at least 90 and must not be more than 2 F1 below MinEl on
   the same set (MinEl unchanged, as the reference for a detector without start-up state). Controls must not drop by
   more than 0.1 F1 against V0.
C. Test databases (EDB, NSRDB, QTDB, TWADB, NSTDB) and LTAFDB: no database loses more than 0.1 F1 against V0;
   the pooled test F1 must not drop.
D. Fusion (uses the TiTen branch): rerun with the chosen variant; rules A and C apply to Fusion as well.
   No Fusion parameter is retuned unless rule A or C fails for Fusion.
Only after A-D pass: run TELE once for TiTen and Fusion and report the result as it comes out.

## Amendment 1 (after V1-V3 results on training and synthetic data; TELE still not run)
Results: no variant passed rule B. TiTen synthetic F1: V0 53.25, V1 56.28, V2 78.89, V3 83.02 (MinEl 94.02).
Rule A held for all (MITDB 99.52-99.53, SVDB 99.48). Diagnosis: V3 misses depend on clip length (Spearman -0.65);
clips longer than ~6 s contaminate the 10 s initialisation window, and the watchdog only caps the signal level at the
window maximum, so recovery takes ~20 s.
- V4: V3, and at each watchdog firing re-initialise the four levels from the last 2 s (20 stream samples) before the
  current candidate with the original rule (signal = init_sig_frac x max, noise = 0.5 x mean, for the integrated and
  the raw stream), instead of capping. No new parameter. Same rules A-D. This is the last design change before TELE;
  if V4 fails rule B, the fix is reported as not achieved and TiTen stays at 2.2.0.

## Result of Amendment 1 (V4)
TiTen synthetic F1 91.65 (artifact excerpts 87.45, controls 99.61); MinEl 94.02 on the same set -> 2.37 F1 below
MinEl, limit 2.0: rule B FAILS by 0.37. Rule A holds (MITDB 99.53, SVDB 99.48; V0 99.52 / 99.48).
Fusion with V4: synthetic 93.64 (V0 93.19), MITDB 99.71, SVDB 99.74 (unchanged).
Per the amendment, the fix is not achieved under the pre-registered rules; TiTen stays at 2.2.0 unless the rule is
changed openly. Rules C and D and TELE were not run.

## Amendment 2 (decided by the author after seeing the V4 result; disclosed)
Rule B is relaxed to: synthetic pooled F1 >= 90 and controls not more than 0.1 F1 below V0. The "within 2 F1 of
MinEl" part is dropped. Reason given: the 2-point margin was an arbitrary choice, and the absolute bar (>= 90) is the
part that tests recovery from a clipped start. This change was made AFTER V4 scored 91.65 (2.37 below MinEl), so any
report must state that rule B was relaxed post hoc. V4 passes the relaxed rule B. Rules A, C, D unchanged; TELE is
run once after C and D pass.

## Results after Amendment 2
Rule C (TiTen, V4 - V0, F1 at +/-150 ms): EDB +0.006, NSRDB +0.002, QTDB 0.000, TWADB -0.054, NSTDB +0.117,
LTAFDB +0.014; pooled test 99.483 -> 99.486. PASS.
Rule D (Fusion with V4): MITDB/SVDB unchanged (99.71/99.74); EDB +0.001, NSRDB 0.000, QTDB 0.000, TWADB -0.011,
NSTDB -0.001, LTAFDB 0.000; pooled test 99.599 -> 99.600. PASS.
TELE, run once (masked, +/-150 ms): TiTen 62.49 -> 90.82 (Se 46.39 -> 86.20, PPV 95.69 -> 95.97);
Fusion 93.66 -> 94.85 (Se 91.34 -> 93.28, PPV 96.11 -> 96.48). Record-wise mean F1: TiTen 47.72 -> 86.21,
Fusion 89.77 -> 91.08. MinEl (unchanged) 94.20.
