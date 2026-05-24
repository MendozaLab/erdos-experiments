# EHP114 n=14 Normal-Collar Critical-Point Exclusion Pilot

Experiment: `EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-CELL-01-01-20260506-01`

Source: `EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01`

## Verdict

- Status: `"NORMAL_COLLAR_FAIL_CRITICAL_POINT_EXCLUSION"`
- Processed regions: `64`
- Normal-collar certified regions: `0`
- Excluded regions: `0`
- Critical-exclusion failures: `56`
- Wall-separation failures: `8`
- Remaining unresolved regions: `64`
- Resolved length upper: `-0.0`
- Total validated length upper: `18.308490790487795`
- Exact length cap: `20.672796062619668`
- Margin to cap: `2.3643052721318725`
- First failed condition: `"normal_wall_sign_separation_failed"`

## Interpretation

This run rotates each unresolved region into midpoint-gradient normal/tangent coordinates. It tests whether the level curve can be certified as a graph over the tangent direction by bounding the normal derivative away from zero and showing opposite signed normal walls.

## Next Blocker

Normal derivative cannot be bounded away from zero on some rotated collars. Next route is analytic critical-point exclusion or higher-order Taylor remainder control.

## Claim Ceiling

Local n=14 hard-cell normal-collar critical-point exclusion pilot only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
