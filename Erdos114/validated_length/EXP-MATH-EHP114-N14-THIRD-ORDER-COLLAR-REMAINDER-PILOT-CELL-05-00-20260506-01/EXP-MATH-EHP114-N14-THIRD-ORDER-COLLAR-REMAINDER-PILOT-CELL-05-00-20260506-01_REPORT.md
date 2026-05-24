# EHP114 n=14 Third-Order Root-Collar Remainder Pilot

Experiment: `EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-05-00-20260506-01`

Source: `EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01`

## Verdict

- Status: `"THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"`
- Processed regions: `64`
- Third-order certified regions: `0`
- Critical-exclusion passes: `15`
- Wall-separation passes: `0`
- Remaining unresolved regions: `64`
- Resolved length upper: `-0.0`
- Total validated length upper: `16.022905622722618`
- Exact length cap: `20.672796062619668`
- Margin to cap: `4.6498904398970495`
- First failed condition: `"third_order_wall_separation_failed"`

## Failing Inequality

sup_r |F(0,r)| + normal_remainder = 0.9644581177979158 is not below S*lower(|F_n|) = 0.03499913103312247; margin -0.9294589867647933

## Interpretation

This run keeps the midpoint-gradient normal/tangent coordinates from the normal-collar artifact and replaces the first-order collar test with an explicit Taylor remainder check. The recurrence computes `p`, `p'`, `p''`, and `p'''`; the logged third-directional upper bound is used to audit whether critical-point exclusion and wall separation can be certified on each residual region.

## Next Blocker

The Taylor critical-exclusion inequality still does not close. Retire sampled local-collar escalation and write the analytic critical-point theorem target.

## Claim Ceiling

Local n=14 hard-cell third-order collar remainder pilot only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
