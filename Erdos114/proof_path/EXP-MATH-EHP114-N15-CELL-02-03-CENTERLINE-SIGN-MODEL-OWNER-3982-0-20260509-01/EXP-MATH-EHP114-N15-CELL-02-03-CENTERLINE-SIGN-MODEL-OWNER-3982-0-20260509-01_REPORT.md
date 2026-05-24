# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `3982:0`'s
16 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_PARTIAL_PASS_3982_0`

## Headline numbers

- Failures in family: `16`
- Closed by WS-01 (conservative 2R bound): `4`
- Closed by WS-01 (tight R only): `11`
- Still failing: `12`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.005520414975142834`
- Median required factor (RHS/LHS) after WS-01 at 2R: `0.7574843985487314`

## What WS-01 does

The original wall test bounds `sup_r |F(0,r)|` by the raw absolute interval
hull of `|F|` on the center strip. That hull pays for full root-affine
dependency at once and is dominated by the linear-in-R coefficient
`|F_t(box-midpoint)|` which can be O(1).

WS-01 picks a validated zero point `z*(u)` on the F=0 curve in each box and
uses the tangent direction `t(u)` to that curve as the new center-strip
coordinate. By construction `F(z*) = 0` and `F_t(z*) = 0`, so the linear-in-R
term vanishes:

```
new_LHS = 0.5 |F_tt| R^2 + (T3/6) R^3 + R_n
```

Branch-point validation per box uses two interval-arithmetic facts already
recorded in the source boundary-slice run:

1. `f_interval` straddles zero (IVT => F has a zero in the box).
2. `center_fn_abs_lower > 0` (grad F nonvanishing => zero curve is a smooth
   1-manifold => t(u) is well-defined).

`F_t(z*) = 0` then holds by construction, not by any numerical estimate.

`R` is taken as `tangent_radius_upper` from the source data; we report the
conservative `2R` figure as the canonical pass condition (worst case z*
sitting near a box corner). The tight-R figure is reported as a sensitivity
diagnostic.

## Per-failure outcomes (owner family 3982:0)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 2476 | `root/y0/y0/y0/y0` | 0.001767 | 0.4012 | 0.00229 | 0.002803 | `CLOSED_BY_WS01` |
| 1 | 2476 | `root/y0/y0/y0/y1` | 0.001763 | 0.3894 | 0.0023 | 0.00265 | `CLOSED_BY_WS01` |
| 2 | 2476 | `root/y0/y0/y1/y0` | 0.001757 | 0.3775 | 0.002307 | 0.002499 | `CLOSED_BY_WS01` |
| 3 | 2476 | `root/y0/y0/y1/y1` | 0.001751 | 0.3653 | 0.00231 | 0.00235 | `CLOSED_BY_WS01` |
| 4 | 2476 | `root/y0/y1/y0/y0` | 0.001744 | 0.3529 | 0.00231 | 0.002204 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 5 | 2476 | `root/y0/y1/y0/y1` | 0.001736 | 0.3404 | 0.002305 | 0.002061 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 6 | 2476 | `root/y0/y1/y1/y0` | 0.001727 | 0.3276 | 0.002296 | 0.001922 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 7 | 2476 | `root/y0/y1/y1/y1` | 0.001717 | 0.3146 | 0.002281 | 0.001786 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 8 | 2476 | `root/y1/y0/y0/y0` | 0.001705 | 0.3086 | 0.002261 | 0.001655 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 9 | 2476 | `root/y1/y0/y0/y1` | 0.001693 | 0.3085 | 0.002236 | 0.001529 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 10 | 2476 | `root/y1/y0/y1/y0` | 0.00168 | 0.3081 | 0.002205 | 0.001408 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 11 | 2476 | `root/y1/y0/y1/y1` | 0.001665 | 0.3084 | 0.002177 | 0.001286 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 12 | 2476 | `root/y1/y1/y0/y0` | 0.00165 | 0.3116 | 0.002158 | 0.00115 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 13 | 2476 | `root/y1/y1/y0/y1` | 0.001633 | 0.3155 | 0.002152 | 0.001006 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 14 | 2476 | `root/y1/y1/y1/y0` | 0.001614 | 0.319 | 0.002139 | 0.0008671 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 15 | 2476 | `root/y1/y1/y1/y1` | 0.001595 | 0.3221 | 0.002129 | 0.0007354 | `STILL_FAILING_AFTER_WS01` |

## Interpretation

On owner family 3982:0 of CELL-02-03's n=15 wall-separation failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~181.15x miss to a ~0.76x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 4/16 boxes close by the conservative 2R bound and 11 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

Partial pass on owner family 3982:0. Track B1 viable but not uniform; the failing boxes need either smaller cells, a refined T3 bound, or a Rust interval-certified z* via interval Newton. Recommended next: characterize the still-failing rows (which inequality term dominates? f_tt_upper * R^2, T3 * R^3, or Rn?) and decide whether to refine cells or upgrade to interval-Newton.

## Source provenance

- Boundary-slice source: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03/EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.json`
- SHA-256 status: `PASS`
- Wall-separation-target packet: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01/EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01_RESULTS.json`
- SHA-256 status: `PASS`
- Math reference: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md` (lines 122-146)

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets, the existing
4488:3 centerline-sign-model artifact) were not attempted. The build script
`build_n15_cell0203_centerline_sign_model_replicate.py` writes only this
experiment folder. No prior script or artifact was modified.
