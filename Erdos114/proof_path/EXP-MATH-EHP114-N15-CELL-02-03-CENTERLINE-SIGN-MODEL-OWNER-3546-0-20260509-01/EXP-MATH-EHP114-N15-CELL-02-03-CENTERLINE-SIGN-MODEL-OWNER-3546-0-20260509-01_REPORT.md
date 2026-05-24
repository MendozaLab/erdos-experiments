# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3546-0-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `3546:0`'s
16 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_PARTIAL_PASS_3546_0`

## Headline numbers

- Failures in family: `16`
- Closed by WS-01 (conservative 2R bound): `6`
- Closed by WS-01 (tight R only): `10`
- Still failing: `10`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.005025676692895417`
- Median required factor (RHS/LHS) after WS-01 at 2R: `0.8439888666232092`

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

## Per-failure outcomes (owner family 3546:0)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 1926 | `root/y0/y0/y0/y0` | 0.001769 | 0.3365 | 0.002603 | 0.003518 | `CLOSED_BY_WS01` |
| 1 | 1926 | `root/y0/y0/y0/y1` | 0.001772 | 0.3396 | 0.002519 | 0.003242 | `CLOSED_BY_WS01` |
| 2 | 1926 | `root/y0/y0/y1/y0` | 0.001775 | 0.3424 | 0.002437 | 0.002985 | `CLOSED_BY_WS01` |
| 3 | 1926 | `root/y0/y0/y1/y1` | 0.001777 | 0.3449 | 0.00236 | 0.002748 | `CLOSED_BY_WS01` |
| 4 | 1926 | `root/y0/y1/y0/y0` | 0.001779 | 0.3477 | 0.002294 | 0.002523 | `CLOSED_BY_WS01` |
| 5 | 1926 | `root/y0/y1/y0/y1` | 0.001781 | 0.3515 | 0.002244 | 0.002304 | `CLOSED_BY_WS01` |
| 6 | 1926 | `root/y0/y1/y1/y0` | 0.001783 | 0.3554 | 0.002204 | 0.002096 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 7 | 1926 | `root/y0/y1/y1/y1` | 0.001784 | 0.359 | 0.002165 | 0.001903 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 8 | 1926 | `root/y1/y0/y0/y0` | 0.001785 | 0.3624 | 0.002128 | 0.001722 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 9 | 1926 | `root/y1/y0/y0/y1` | 0.001787 | 0.3655 | 0.002092 | 0.001554 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 10 | 1926 | `root/y1/y0/y1/y0` | 0.001788 | 0.3684 | 0.002058 | 0.001398 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 11 | 1926 | `root/y1/y0/y1/y1` | 0.001789 | 0.3711 | 0.002026 | 0.001251 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 12 | 1926 | `root/y1/y1/y0/y0` | 0.00179 | 0.3735 | 0.001995 | 0.001115 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 13 | 1926 | `root/y1/y1/y0/y1` | 0.001791 | 0.3758 | 0.001965 | 0.0009872 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 14 | 1926 | `root/y1/y1/y1/y0` | 0.001791 | 0.3778 | 0.001937 | 0.000868 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 15 | 1926 | `root/y1/y1/y1/y1` | 0.001792 | 0.3797 | 0.001909 | 0.0007566 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |

## Interpretation

On owner family 3546:0 of CELL-02-03's n=15 wall-separation failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~198.98x miss to a ~0.84x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 6/16 boxes close by the conservative 2R bound and 10 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

Partial pass on owner family 3546:0. Track B1 viable but not uniform; the failing boxes need either smaller cells, a refined T3 bound, or a Rust interval-certified z* via interval Newton. Recommended next: characterize the still-failing rows (which inequality term dominates? f_tt_upper * R^2, T3 * R^3, or Rn?) and decide whether to refine cells or upgrade to interval-Newton.

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
