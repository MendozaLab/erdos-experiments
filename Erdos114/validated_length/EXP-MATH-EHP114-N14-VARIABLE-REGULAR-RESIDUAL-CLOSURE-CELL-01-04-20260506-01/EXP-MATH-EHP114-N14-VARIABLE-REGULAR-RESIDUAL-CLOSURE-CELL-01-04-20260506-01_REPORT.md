# EHP114 n=14 Regular Residual Decomposition

Experiment: `EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-01-04-20260506-01`

Source: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-01-04-20260506-01`

## Verdict

- Status: `"REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF"`
- Processed regular regions: `12`
- Certified regular regions: `0`
- Excluded regular regions: `12`
- Wall separation failures: `0`
- Ownership duplicates: `0`
- Adaptive max depth: `5`
- Adaptive leaf boxes: `36`
- Adaptive unresolved leaves: `0`
- Regular slice length upper: `-0.0`
- Total validated length upper: `15.53773753688207`
- Exact length cap: `20.672796062619668`
- Margin to cap: `5.135058525737598`
- First failed condition: `"none"`

## Interpretation

This diagnostic attempts the x-dominant regular regions found by the L21 global critical-point diagnostic. It does not process or promote the critical-candidate regions. When `max_depth` is positive, failed regular walls may be subdivided in the y direction, but closure is still reported against the original owned source region. A pass means the regular slice has a theorem-shaped x-chart certificate with wall separation or monotone exclusion and ownership; it is still not a global n=14 certificate.

## Claim Ceiling

Local n=14 hard-cell regular residual decomposition diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
