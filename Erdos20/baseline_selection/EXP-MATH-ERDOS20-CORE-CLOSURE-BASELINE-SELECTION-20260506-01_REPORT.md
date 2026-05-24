# EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01 Report

**Scope:** saved-artifact diagnostic only. No new enumeration, no D1, no scorecard, no git, no public docs.

**Recommendation:** `floor_normalized_I_core_local` as the next-run baseline, with observed `I_core_local(s,m)` as the numerator and AHS as a secondary external control.

## Ranking

| Candidate | Quantity | Score | Meaning |
|---|---:|---:|---|
| floor_normalized_I_core_local | `I_core_local(s,m) / I_floor(s,m)` | 9 | Best Leg-4 target if I_floor is predeclared before the run; preserves the measured core channel. |
| observed_I_core_local | `I_core_local(s,m)` | 7 | Best immediate numerator, but not a Leg-4 baseline because it lacks a denominator. |
| ahs_normalized_closure_cost | `I_core_local(s,m) / I_AHS(s,m) or closure excess over an AHS-style construction baseline` | 5 | Necessary external control, but not ready as the primary numerator because no AHS construction artifact exists here. |

## Current Signal Used

The best saved precursor remains `w=3, fixed core size s=2` with classification `STABLE_BY_SLOPE_ONLY`, slope `-0.076613`, and late-window CV `0.182257`.

## Claim Ceiling

A-axis A0; shadow signature, not universal law; no theorem, no lower-bound progress, no Leg-4 pass.
