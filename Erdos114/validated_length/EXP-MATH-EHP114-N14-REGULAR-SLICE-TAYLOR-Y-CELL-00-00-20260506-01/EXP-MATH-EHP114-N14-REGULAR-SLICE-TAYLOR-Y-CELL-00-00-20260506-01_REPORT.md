# EHP114 n=14 Regular Residual Decomposition

Experiment: `EXP-MATH-EHP114-N14-REGULAR-SLICE-TAYLOR-Y-CELL-00-00-20260506-01`

Source: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01`

## Verdict

- Status: `"REGULAR_RESIDUAL_FAIL_WALL_SEPARATION"`
- Processed regular regions: `8`
- Certified regular regions: `0`
- Excluded regular regions: `0`
- Wall separation failures: `8`
- Ownership duplicates: `0`
- Regular slice length upper: `-0.0`
- Total validated length upper: `15.410486279252329`
- Exact length cap: `20.672796062619668`
- Margin to cap: `5.262309783367339`
- First failed condition: `"2286:0:root/y0/y0 failed because x-wall F intervals do not have strict opposite signs"`

## Interpretation

This diagnostic only attempts the eight x-dominant regular regions found by the L21 global critical-point diagnostic. It does not process or promote the 56 critical-candidate regions. A pass means the regular slice has a theorem-shaped x-chart certificate with wall separation and ownership; it is still not a global n=14 certificate.

## Claim Ceiling

Local n=14 hard-cell regular residual decomposition diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
