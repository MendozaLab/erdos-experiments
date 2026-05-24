# EHP114 n=14 Global Critical-Point Exclusion Target

Experiment: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-CELL-01-07-20260506-01`

Source: `EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01`

## Verdict

- Status: `"GLOBAL_CRITICAL_POINT_TARGET_REGULAR_REGIONS_FOUND"`
- Processed regions: `64`
- Regular regions: `8`
- Critical-candidate regions: `56`
- Excluded regions: `0`
- Gradient lower-bound candidate: `34.42575161732093`
- Total validated length upper: `16.033405318919826`
- Exact length cap: `20.672796062619668`
- Margin to cap: `4.639390743699842`
- First failed condition: `"2484:2:root/y0/y0 remains critical candidate because both Fx and Fy intervals contain zero under root-affine uncertainty"`

## Interpretation

This diagnostic asks whether any L18 residual region already admits a global root-affine lower bound for `|grad F|` in ordinary coordinates. It does not promote branch length. Regular regions are candidates for a future analytic domain decomposition; critical-candidate regions identify where a sharper critical-point theorem is still needed.

## Claim Ceiling

Local n=14 hard-cell global critical-point exclusion target diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
