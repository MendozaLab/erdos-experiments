# EHP114 n=14 Direct Validated-Length Hard-Cell Diagnostic

Experiment: `EXP-MATH-EHP114-N14-VALIDATED-LENGTH-HARD-CELL-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

## Verdict

- Status: `VALIDATED_LENGTH_FAIL_BUDGET`
- Degree: `14`
- eps: `0.1`
- Root-affine subcell: `(6, 4)`
- Exact length cap: `20.672796062619668`
- Outside extent excluded: `true` (margin `1.0071176513872322`)
- Patch count: `49858`
- Excluded boxes: `12340542`
- Unresolved boxes: `0`
- Ownership duplicate count: `0`
- Sum patch length upper: `119.57917312364832`
- Endpoint/overlap tax: `0.0`
- Total validated length upper: `119.57917312364832`
- Margin to cap: `-98.90637706102865`

## Interpretation

This run demotes marching squares to a diagnostic and attempts a direct interval implicit-function graph enclosure. Normal drift is reported only as an explanatory diagnostic; it is not the proof anchor. The prior `relative error / budget = 2.730675196761203` is therefore no longer load-bearing for the pass/fail decision.

## Claim Ceiling

This is a local hard-cell diagnostic. It is not a proof of Erdős #114 and not a global n=14 proof.

## Next Blocker

The direct graph-patch upper bound is certified but too coarse for the exact-length cap. Need smaller boxes, tighter derivative bounds, or a sharper patch integration rule.
