# EHP114 n=14 Third-Order Root-Collar Remainder Pilot

Experiment: `EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-04-01-20260506-01`

Source: `EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01`

## Verdict

- Status: `"THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"`
- Processed regions: `64`
- Third-order certified regions: `0`
- Critical-exclusion passes: `22`
- Wall-separation passes: `0`
- Remaining unresolved regions: `64`
- Resolved length upper: `-0.0`
- Total validated length upper: `16.5271690923772`
- Exact length cap: `20.672796062619668`
- Margin to cap: `4.145626970242468`
- First failed condition: `"third_order_wall_separation_failed"`

## Failing Inequality

sup_r |F(0,r)| + normal_remainder = 0.9644188626456115 is not below S*lower(|F_n|) = 0.03498348366033147; margin -0.9294353789852801

## Interpretation

This run keeps the midpoint-gradient normal/tangent coordinates from the normal-collar artifact and replaces the first-order collar test with an explicit Taylor remainder check. The recurrence computes `p`, `p'`, `p''`, and `p'''`; the logged third-directional upper bound is used to audit whether critical-point exclusion and wall separation can be certified on each residual region.

## Next Blocker

The Taylor critical-exclusion inequality still does not close. Retire sampled local-collar escalation and write the analytic critical-point theorem target.

## Claim Ceiling

Local n=14 hard-cell third-order collar remainder pilot only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
