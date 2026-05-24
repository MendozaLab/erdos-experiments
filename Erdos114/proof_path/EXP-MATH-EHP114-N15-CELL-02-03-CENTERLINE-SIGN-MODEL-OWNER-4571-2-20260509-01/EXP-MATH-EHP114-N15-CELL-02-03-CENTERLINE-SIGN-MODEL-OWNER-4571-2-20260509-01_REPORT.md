# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `4571:2`'s
11 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_PARTIAL_PASS_4571_2`

## Headline numbers

- Failures in family: `11`
- Closed by WS-01 (conservative 2R bound): `4`
- Closed by WS-01 (tight R only): `1`
- Still failing: `7`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.0039039881908065933`
- Median required factor (RHS/LHS) after WS-01 at 2R: `0.36307964164003076`

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

## Per-failure outcomes (owner family 4571:2)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 3852 | `root/y0/y0/y0/y0` | 0.001932 | 1.649 | 0.0186 | 0.000861 | `STILL_FAILING_AFTER_WS01` |
| 1 | 3852 | `root/y0/y0/y0/y1` | 0.00188 | 1.596 | 0.01856 | 0.0003842 | `STILL_FAILING_AFTER_WS01` |
| 2 | 3852 | `root/y0/y0/y1/y0` | 0.001822 | 1.529 | 0.01818 | 8.063e-05 | `STILL_FAILING_AFTER_WS01` |
| 3 | 3852 | `root/y1/y0/y0/y0` | 0.001353 | 0.9188 | 0.01047 | 0.0005001 | `STILL_FAILING_AFTER_WS01` |
| 4 | 3852 | `root/y1/y0/y0/y1` | 0.001254 | 0.7958 | 0.00881 | 0.001475 | `STILL_FAILING_AFTER_WS01` |
| 5 | 3852 | `root/y1/y0/y1/y0` | 0.001149 | 0.6717 | 0.007222 | 0.002622 | `STILL_FAILING_AFTER_WS01` |
| 6 | 3852 | `root/y1/y0/y1/y1` | 0.001038 | 0.5484 | 0.005761 | 0.003868 | `CLOSED_BY_WS01_TIGHT_R_ONLY` |
| 7 | 3852 | `root/y1/y1/y0/y0` | 0.0009226 | 0.4364 | 0.004464 | 0.005126 | `CLOSED_BY_WS01` |
| 8 | 3852 | `root/y1/y1/y0/y1` | 0.0008014 | 0.4148 | 0.003394 | 0.006281 | `CLOSED_BY_WS01` |
| 9 | 3852 | `root/y1/y1/y1/y0` | 0.0006752 | 0.392 | 0.002513 | 0.007307 | `CLOSED_BY_WS01` |
| 10 | 3852 | `root/y1/y1/y1/y1` | 0.0005441 | 0.3687 | 0.001807 | 0.008152 | `CLOSED_BY_WS01` |

## Interpretation

On owner family 4571:2 of CELL-02-03's n=15 wall-separation failures (11 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~256.15x miss to a ~0.36x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 4/11 boxes close by the conservative 2R bound and 1 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

Partial pass on owner family 4571:2. Track B1 viable but not uniform; the failing boxes need either smaller cells, a refined T3 bound, or a Rust interval-certified z* via interval Newton. Recommended next: characterize the still-failing rows (which inequality term dominates? f_tt_upper * R^2, T3 * R^3, or Rn?) and decide whether to refine cells or upgrade to interval-Newton.

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
