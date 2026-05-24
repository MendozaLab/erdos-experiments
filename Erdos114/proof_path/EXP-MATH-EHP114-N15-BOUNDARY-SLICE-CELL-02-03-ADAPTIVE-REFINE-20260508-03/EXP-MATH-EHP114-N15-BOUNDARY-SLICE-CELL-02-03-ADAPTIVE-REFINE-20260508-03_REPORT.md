# EHP114 n=15 CELL-02-03 Boundary-Slice Collar Pilot

Experiment: `EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03`

Source: `EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-02`

## Verdict

- Final verdict: `"BOUNDARY_SLICE_PARTIAL_NEEDS_REDUCTION"`
- Collar status: `"THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"`
- Processed regions: `285`
- Third-order certified regions: `0`
- Critical-exclusion passes: `108`
- Wall-separation passes: `0`
- Remaining unresolved regions: `177`
- Resolved length upper: `-0.0`
- Total validated length upper: `17.64690787614151`
- Exact length cap: `20.672796062619668`
- Margin to cap: `3.0258881864781593`
- First failed condition: `"third_order_wall_separation_failed"`

## Failing Inequality

sup_r |F(0,r)| + normal_remainder = 0.37211795118157914 is not below S*lower(|F_n|) = 0.010856799748727618; margin -0.36126115143285153

## Interpretation

This run reuses the n=14 CELL-02-03 residual boxes as a spatial testbed, but evaluates them with a degree-15 two-mode Fourier/root-slice chart. The chart axes are `m7_sin_tangent` and `m0_radial_boundary`. This is a tractability probe for the hard collar/root-isolation analogue, not full n=15 certification.

## Next Blocker

The Taylor critical-exclusion inequality still does not close. Build the analytic critical-point blocker target before spending full compute.

## Claim Ceiling

Review-only n=15 CELL-02-03 boundary-slice/collar pilot. Not full n=15 certification, not public state, and not a formal-status promotion.
