# EHP114 n=14 Regular Residual Decomposition

Experiment: `EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-04-04-20260506-01`

Source: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-04-04-20260506-01`

## Verdict

- Status: `"REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF"`
- Processed regular regions: `12`
- Certified regular regions: `1`
- Excluded regular regions: `11`
- Wall separation failures: `0`
- Ownership duplicates: `0`
- Adaptive max depth: `5`
- X-adaptive max depth: `8`
- Adaptive leaf boxes: `219`
- Adaptive unresolved leaves: `0`
- Regular slice length upper: `0.05573343241906941`
- Total validated length upper: `19.25493537193725`
- Exact length cap: `20.672796062619668`
- Margin to cap: `1.4178606906824172`
- First failed condition: `"none"`

## Interpretation

This diagnostic attempts the x-dominant regular regions found by the L21 global critical-point diagnostic. It does not process or promote the critical-candidate regions. When `max_depth` is positive, failed regular walls may be subdivided in the y direction. When `x_adaptive_depth` is positive, failed stable-`Fx` leaves may also be subdivided in the x direction; if strict wall separation still does not close, the stable-`Fx` graph-envelope theorem is allowed as a conservative length upper bound over the owned y-interval. Closure is still reported against the original owned source region, with child split paths preserved where subdivision is load-bearing. A pass means the regular slice has a theorem-shaped x-chart certificate with wall separation, monotone exclusion, or stable-`Fx` graph-envelope ownership; it is still not a global n=14 certificate.

## Claim Ceiling

Local n=14 hard-cell regular residual decomposition diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
