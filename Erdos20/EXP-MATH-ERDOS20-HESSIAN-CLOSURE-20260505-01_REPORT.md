# EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01 - Report

**Date:** 2026-05-05
**Problem:** Erdos #20 sunflower core closure
**Scope:** Aggregate discrete Hessian/free-energy diagnostic from saved April 2026 lattice-gas data
**Classification:** HESSIAN_INCONCLUSIVE
**Claim ceiling:** shadow signature, not universal law

## Meaning

The saved sunflower data contains a real aggregate edge effect: log2 density-of-states bends in the high-family-size tail, and the closure-pressure quantity I_close rises as the mean extension rate approaches jamming. That is the right kind of shadow to inspect for a core-closure observable.

The result is still INCONCLUSIVE because the existing artifacts are aggregate by family size. They do not say which core carried the pressure, which petal channel closed, or whether a per-core Hessian would remain after ordinary geometry is accounted for.

## Method

For each available density-of-states profile D(m), the script computes:

- logD(m) = log2 D(m)
- delta logD(m) = logD(m+1) - logD(m)
- second delta logD(m) = logD(m+1) - 2 logD(m) + logD(m-1), centered at m
- p_safe(m) = g(m) / (N - m), when aggregate growth rates g(m) exist
- I_close(m) = -log2 p_safe(m)

This is a discrete free-energy/closure-pressure proxy. It is not a per-core measurement and not a Hessian of the actual sunflower core state.

## Growth-Rate Profiles

| w | n | N | M | jamming m* | tail min second logD | I_close at m* | max finite I_close |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 4 | 4 | 4 | 2 | -1.41504 | 1.58496 | 2 |
| 3 | 5 | 10 | 6 | 4 | -1.88864 | 3.54432 | 5.16992 |
| 3 | 6 | 20 | 10 | 6 | -1.71847 | 4.46447 | 7.37938 |
| 3 | 7 | 35 | 12 | 7 | -1.06676 | 5.13676 | 8.16993 |
| 4 | 7 | 35 | 15 | 9 | -1.47236 | 4.77074 | 10.043 |

## Strongest Aggregate Curvature Tails

| w | n | N | M | tail min second logD | growth rates |
| --- | --- | --- | --- | --- | --- |
| 2 | 10 | 45 | 6 | -3.11594 | no |
| 2 | 9 | 36 | 6 | -2.99509 | no |
| 2 | 8 | 28 | 6 | -2.86769 | no |
| 2 | 7 | 21 | 6 | -2.73733 | no |
| 2 | 6 | 15 | 6 | -2.61844 | no |
| 2 | 5 | 10 | 5 | -2.58996 | no |
| 2 | 4 | 6 | 4 | -2.50815 | no |
| 3 | 5 | 10 | 6 | -1.88864 | yes |

## Classification

The classifier returns **HESSIAN_INCONCLUSIVE**.

Why not HESSIAN_PRESENT: the PASS-style signal would need per-core closure instrumentation, a defined floor-normalized ratio, and a geometry-exhausting sweep. None of those are present in the saved April aggregate data.

Why not HESSIAN_FAIL: the aggregate profiles do show curvature and closure pressure near jamming, so the fantasy is not empty. It has a measurable bulk trace worth instrumenting properly.

Blocking facts:

- available data are aggregate by family size, not per-core closure states
- growth rates are missing for w=2 and for w=3,n=8
- w=4 stops at n=7 under exhaustive enumeration
- no Mendoza-floor numerator or per-core R profile is defined in these inputs
- Abbott-Hansen-Sauer dominates the small-n lower-bound framing; A-axis remains A0

## Leg-4 Relation

Core closure asks how much bookkeeping the shared intersection must carry so petals are not treated as independent fragments. The aggregate I_close profile measures the cost of finding a safe unused addition at family size m. That makes it a useful precursor to the Leg-4 observable.

It does not execute Leg 4. A real core-closure Leg-4 run needs I_core(C,s,m) distributions by core size, plus a predeclared Mendoza-floor or construction-normalized comparison.

## Claim Limits

A-axis remains A0. Abbott-Hansen-Sauer dominates the small-n lower-bound framing, so these measurements are not a new lower-bound story. The script computes a diagnostic observable; it does not close the conjecture, change the public status of #20, or justify publication language.

## Artifacts

- `EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_RESULTS.json`
- `EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_REPORT.md`
- `EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_RESULTS.sha256`
- `sunflower_core_hessian_analysis.py`
