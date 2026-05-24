# EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01 - Report

**Date:** 2026-05-05
**Problem:** Erdos #20 sunflower core closure
**Scope:** Exact per-core closure instrumentation on cheap calibration regimes
**Classification:** PER_CORE_SIGNAL_PRESENT
**Claim ceiling:** shadow signature, not universal law

## Meaning

The aggregate Hessian run could see jamming, but it could not say where the closure pressure lived. This run adds that missing channel for the cheap exact regimes: a core is an actual subset C, and a petal extension is a candidate w-set of the form C union P.

The per-core signal is present in this limited sense: near the aggregate jamming point, fixed-core channels become locally expensive, and the expense depends on core size. That makes the earlier aggregate trace more specific. It is still only a precursor, because the expensive regimes and the floor-normalized comparison are not executed here.

## Operational Definitions

- `core_id`: the sorted 1-based elements of a subset C of [n].
- `petal channel`: fixed C with candidate petals P in the complement of C, so the candidate set is C union P.
- `I_core_local(C,s,m)`: `-log2(local valid petal extensions / unused petal extensions through C)` at family size m.
- `I_core_global(C,s,m)`: the same denominator, but the numerator requires the candidate extension to be globally sunflower-free.
- `aggregate I_close`: the prior continuity check using `D(m+1)/D(m)` divided by unused ambient sites.

## Executed Targets

| w | n | N | families | M | m* | aggregate I | all-unused I | strongest s | strongest I_core | seconds |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 6 | 15 | 1198 | 6 | 4 | 5.53138 | 3.20945 | 1 | 0.732247 | 0.01 |
| 2 | 7 | 21 | 4264 | 6 | 4 | 5.80171 | 3.47978 | 0 | 0.779341 | 0.04 |
| 2 | 8 | 28 | 12265 | 6 | 4 | 6.08443 | 3.7625 | 0 | 1.02296 | 0.128 |
| 2 | 9 | 36 | 30193 | 6 | 4 | 6.35484 | 4.03292 | 0 | 1.22081 | 0.37 |
| 2 | 10 | 45 | 66130 | 6 | 4 | 6.60793 | 4.286 | 0 | 1.38845 | 0.91 |
| 3 | 4 | 4 | 16 | 4 | 2 | 1.58496 | 0 | 0 | 0 | 0 |
| 3 | 5 | 10 | 388 | 6 | 4 | 3.54432 | 1.22239 | 2 | 0.326228 | 0.006 |
| 3 | 6 | 20 | 33652 | 10 | 6 | 4.46447 | 1.65711 | 2 | 0.459083 | 0.798 |
| 3 | 7 | 35 | 2485795 | 12 | 7 | 5.13675 | 2.13675 | 2 | 0.354808 | 83.538 |
| 4 | 5 | 5 | 32 | 5 | 3 | 2 | 0 | 0 | 0 | 0.001 |
| 4 | 6 | 15 | 5789 | 9 | 5 | 3.50658 | 0.921612 | 3 | 0.199987 | 0.257 |

## Core-Size Slice At Target m

| w | n | m | s | petal | channels | local I | global I | median finite I | blocked frac |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 10 | 4 | 0 | 2 | 37170 | 1.38845 | 4.286 | 1.45066 | 0 |
| 2 | 10 | 4 | 1 | 1 | 371700 | 0.347652 | 4.286 | 0 | 0.250847 |
| 3 | 7 | 7 | 0 | 3 | 801440 | 0 | 2.13675 | 0 | 0 |
| 3 | 7 | 7 | 1 | 2 | 5610080 | 0.215674 | 2.13675 | 0.137504 | 0.00025 |
| 3 | 7 | 7 | 2 | 1 | 16830240 | 0.354808 | 2.13675 | 0 | 0.290702 |
| 4 | 6 | 5 | 0 | 4 | 1773 | 0 | 0.921612 | 0 | 0 |
| 4 | 6 | 5 | 1 | 3 | 10638 | 0 | 0.921612 | 0 | 0 |
| 4 | 6 | 5 | 2 | 2 | 26595 | 0 | 0.921612 | 0 | 0 |
| 4 | 6 | 5 | 3 | 1 | 35460 | 0.199987 | 0.921612 | 0 | 0.258883 |

## Classification

The classifier returns **PER_CORE_SIGNAL_PRESENT**.

This means the per-core channel is no longer blocked by missing instrumentation for the cheap regimes. It does not mean a Leg-4 pass. A pass would require the floor-normalized numerator, geometry controls, and a larger symmetry-reduced sweep.

The comparison to aggregate `I_close` should be read narrowly. Aggregate `I_close` asks how the density of states changes with family size. `I_core` asks how a fixed core channel closes as petals accumulate. Agreement in the late regime is a useful precursor; mismatch is expected because the denominators are different.

## Blockers

- `w=3,n=8` and `w=4,n=7` are too large for this exact Python per-core pass.
- `w=4,n=8` was already beyond the prior exhaustive aggregate run.
- No Mendoza-floor or construction-normalized numerator is predeclared here.
- No Abbott-Hansen-Sauer baseline family is instrumented as a control.
- The runner is exact enumeration, not the transfer-matrix rewrite needed for scale.

## Claim Limits

A-axis remains A0. The artifact measures an internal diagnostic channel. It does not change the status of Erdos #20, does not give new lower-bound progress, and does not justify public power-morphism language.

## Artifacts

- `EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json`
- `EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_REPORT.md`
- `EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.sha256`
- `sunflower_per_core_closure_analysis.py`
