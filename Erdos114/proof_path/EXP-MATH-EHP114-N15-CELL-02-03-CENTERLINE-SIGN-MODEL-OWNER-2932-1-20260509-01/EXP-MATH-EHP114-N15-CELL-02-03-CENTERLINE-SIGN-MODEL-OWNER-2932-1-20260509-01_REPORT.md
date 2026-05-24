# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `2932:1`'s
1 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_FULL_PASS_2932_1`

## Headline numbers

- Failures in family: `1`
- Closed by WS-01 (conservative 2R bound): `1`
- Closed by WS-01 (tight R only): `0`
- Still failing: `0`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.008529360272181926`
- Median required factor (RHS/LHS) after WS-01 at 2R: `2.0443904669620077`

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

## Per-failure outcomes (owner family 2932:1)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 1100 | `root/y0/y0/y0/y0` | 0.001079 | 0.7115 | 0.002968 | 0.006069 | `CLOSED_BY_WS01` |

## Interpretation

On owner family 2932:1 of CELL-02-03's n=15 wall-separation failures (1 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~117.24x miss to a ~2.04x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 1/1 boxes close by the conservative 2R bound and 0 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

Track B1-replicate confirmed full pass on owner family 2932:1. Combined with the 4488:3 result (16/16 closed), Track B1 viable across the worst-three owner families. Recommended next: replicate on the remaining matched-count families (3104:0, 3546:0, 3982:0, 4404:0 with 16 each), then move toward Rust interval-certification of the moving-frame collar that emits z* via interval Newton on F(z)=0.

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
