# EHP114 n=14 Regular Residual Decomposition

Experiment: `EXP-MATH-EHP114-N14-X-ADAPTIVE-REGULAR-RESIDUAL-CLOSURE-CELL-02-07-20260506-03`

Source: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-02-07-20260506-01`

## Verdict

- Status: `"REGULAR_RESIDUAL_FAIL_WALL_SEPARATION"`
- Processed regular regions: `12`
- Certified regular regions: `0`
- Excluded regular regions: `11`
- Wall separation failures: `1`
- Ownership duplicates: `0`
- Adaptive max depth: `5`
- X-adaptive max depth: `12`
- Adaptive leaf boxes: `7082`
- Adaptive unresolved leaves: `6902`
- Regular slice length upper: `-0.0`
- Total validated length upper: `18.789609307756976`
- Exact length cap: `20.672796062619668`
- Margin to cap: `1.8831867548626917`
- First failed condition: `"2716:0:root/y0/y1 failed because x-wall F intervals do not have strict opposite signs; adaptive y-subdivision left 6902 of 7057 leaf boxes unresolved"`

## Interpretation

This diagnostic attempts the x-dominant regular regions found by the L21 global critical-point diagnostic. It does not process or promote the critical-candidate regions. When `max_depth` is positive, failed regular walls may be subdivided in the y direction. When `x_adaptive_depth` is positive, failed stable-`Fx` leaves may also be subdivided in the x direction; closure is still reported against the original owned source region, with child split paths preserved for audit. A pass means the regular slice has a theorem-shaped x-chart certificate with wall separation or monotone exclusion and ownership; it is still not a global n=14 certificate.

## Claim Ceiling

Local n=14 hard-cell regular residual decomposition diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
