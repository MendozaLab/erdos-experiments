# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3104-0-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `3104:0`'s
16 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_BLOCKED_3104_0`

## Headline numbers

- Failures in family: `16`
- Closed by WS-01 (conservative 2R bound): `0`
- Closed by WS-01 (tight R only): `16`
- Still failing: `16`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.003451896959881573`
- Median required factor (RHS/LHS) after WS-01 at 2R: `0.678683167940443`

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

## Per-failure outcomes (owner family 3104:0)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 1375 | `root/y0/y0/y0/y0` | 0.001838 | 0.4219 | 0.003397 | 0.003267 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 1 | 1375 | `root/y0/y0/y0/y1` | 0.00186 | 0.428 | 0.003282 | 0.002943 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 2 | 1375 | `root/y0/y0/y1/y0` | 0.00188 | 0.433 | 0.003152 | 0.002654 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 3 | 1375 | `root/y0/y0/y1/y1` | 0.001898 | 0.437 | 0.003009 | 0.002396 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 4 | 1375 | `root/y0/y1/y0/y0` | 0.001914 | 0.4399 | 0.002853 | 0.002165 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 5 | 1375 | `root/y0/y1/y0/y1` | 0.001927 | 0.4417 | 0.002687 | 0.001958 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 6 | 1375 | `root/y0/y1/y1/y0` | 0.001939 | 0.4425 | 0.002514 | 0.001772 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 7 | 1375 | `root/y0/y1/y1/y1` | 0.001948 | 0.4424 | 0.002337 | 0.001603 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 8 | 1375 | `root/y1/y0/y0/y0` | 0.001956 | 0.4414 | 0.002157 | 0.001448 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 9 | 1375 | `root/y1/y0/y0/y1` | 0.001961 | 0.4396 | 0.001978 | 0.001306 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 10 | 1375 | `root/y1/y0/y1/y0` | 0.001964 | 0.4369 | 0.001801 | 0.001175 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 11 | 1375 | `root/y1/y0/y1/y1` | 0.001964 | 0.4335 | 0.001629 | 0.001053 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 12 | 1375 | `root/y1/y1/y0/y0` | 0.001963 | 0.4294 | 0.001464 | 0.0009391 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 13 | 1375 | `root/y1/y1/y0/y1` | 0.001959 | 0.4246 | 0.001306 | 0.0008317 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 14 | 1375 | `root/y1/y1/y1/y0` | 0.001954 | 0.4193 | 0.001158 | 0.0007307 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 15 | 1375 | `root/y1/y1/y1/y1` | 0.001946 | 0.4133 | 0.00102 | 0.0006349 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |

## Interpretation

On owner family 3104:0 of CELL-02-03's n=15 wall-separation failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~289.70x miss to a ~0.68x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 0/16 boxes close by the conservative 2R bound and 16 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

BLOCKED on owner family 3104:0. Moving-frame collar cannot establish a validated branch point and/or all wall LHS values still exceed RHS even with the WS-01 rewrite. The cell decomposition for n=15 needs to be reconsidered (smaller boxes, or a different analytic reduction).

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
