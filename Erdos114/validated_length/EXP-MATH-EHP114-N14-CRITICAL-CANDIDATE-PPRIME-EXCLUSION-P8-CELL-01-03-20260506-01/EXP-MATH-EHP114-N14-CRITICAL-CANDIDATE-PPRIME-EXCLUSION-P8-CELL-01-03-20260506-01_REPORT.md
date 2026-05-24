# EHP114 n=14 Critical-Candidate Affine Diagnostic

Experiment: `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-01-03-20260506-01`

## Verdict

- Status: `"CRITICAL_CANDIDATE_PPRIME_PARTIAL"`
- Processed critical candidates: `52`
- Affine excluded regions: `0`
- Gradient affine regular regions: `11`
- Affine partition closed regions: `11`
- Still critical candidates: `41`
- Param subdivision: `8`
- Test mode: `"pprime"`
- Total validated length upper: `16.813915919651034`
- Exact length cap: `20.672796062619668`
- Margin to cap: `3.8588801429686335`
- First failed condition: `"2484:4:root/y0/y0 remains critical after affine-gradient parameter tiles"`

## Interpretation

This diagnostic tracks the two root-affine parameters linearly and places nonlinear products into interval remainders. It does not promote branch length. It asks whether any of the 56 critical candidates from L21 are actually regular or excluded once root-parameter dependency is partially preserved.

## Claim Ceiling

Local n=14 hard-cell critical-candidate affine diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
