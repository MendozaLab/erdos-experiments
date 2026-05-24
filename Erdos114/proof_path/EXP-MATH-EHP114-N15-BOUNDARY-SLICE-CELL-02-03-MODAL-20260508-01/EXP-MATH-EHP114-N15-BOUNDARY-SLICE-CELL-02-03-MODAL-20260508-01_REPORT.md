# EHP114 n=15 CELL-02-03 Boundary-Slice Collar Pilot

Experiment: `EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-MODAL-20260508-01`

Source: `EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-02-03-20260506-01`

## Verdict

- Final verdict: `"BOUNDARY_SLICE_PARTIAL_NEEDS_REDUCTION"`
- Collar status: `"THIRD_ORDER_FAIL_CRITICAL_EXCLUSION"`
- Processed regions: `64`
- Third-order certified regions: `0`
- Critical-exclusion passes: `4`
- Wall-separation passes: `0`
- Remaining unresolved regions: `54`
- Resolved length upper: `-0.0`
- Total validated length upper: `17.64690787614151`
- Exact length cap: `20.672796062619668`
- Margin to cap: `3.0258881864781593`
- First failed condition: `"third_order_critical_exclusion_failed"`

## Failing Inequality

|F_n(0,0)|_lower - variation_bound = 4.140781566853416 - 47.055308974763044 = -42.91452740790963 <= 0

## Interpretation

This run reuses the n=14 CELL-02-03 residual boxes as a spatial testbed, but evaluates them with a degree-15 two-mode Fourier/root-slice chart. The chart axes are `m7_sin_tangent` and `m0_radial_boundary`. This is a tractability probe for the hard collar/root-isolation analogue, not full n=15 certification.

## Next Blocker

The Taylor critical-exclusion inequality still does not close. Build the analytic critical-point blocker target before spending full compute.

## Claim Ceiling

Review-only n=15 CELL-02-03 boundary-slice/collar pilot. Not full n=15 certification, not public state, and not a formal-status promotion.
