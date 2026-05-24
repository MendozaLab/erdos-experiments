# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01

## Scope

Internal experiment artifact only. Track B1 of the EHP114 bridge program. Not
a proof of Erdos #114, not an n=15 certificate, and not a CELL-02-03 closure.
Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the Branch-Centered
Moving-Frame Collar rewrite) on owner family `4488:3`'s 16 wall-
separation failures.

## Status

`CENTERLINE_SIGN_MODEL_FULL_PASS_4488_3`

## Headline numbers

- Failures in family: `16`
- Closed by WS-01 (conservative 2R bound): `16`
- Closed by WS-01 (tight R only): `0`
- Still failing: `0`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.012902527624671575`
- Median required factor (RHS/LHS) after WS-01 at 2R: `4.17417952861563`

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

## Per-failure outcomes (owner family 4488:3)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 3577 | `root/y0/y0/y0/y0` | 0.0004866 | 0.6605 | 0.001759 | 0.01928 | `CLOSED_BY_WS01` |
| 1 | 3577 | `root/y0/y0/y0/y1` | 0.0005397 | 0.7914 | 0.002117 | 0.01985 | `CLOSED_BY_WS01` |
| 2 | 3577 | `root/y0/y0/y1/y0` | 0.0005908 | 0.93 | 0.002518 | 0.02041 | `CLOSED_BY_WS01` |
| 3 | 3577 | `root/y0/y0/y1/y1` | 0.00064 | 1.076 | 0.002969 | 0.02095 | `CLOSED_BY_WS01` |
| 4 | 3577 | `root/y0/y1/y0/y0` | 0.0006871 | 1.229 | 0.00347 | 0.02149 | `CLOSED_BY_WS01` |
| 5 | 3577 | `root/y0/y1/y0/y1` | 0.0007323 | 1.39 | 0.004031 | 0.02202 | `CLOSED_BY_WS01` |
| 6 | 3577 | `root/y0/y1/y1/y0` | 0.0007756 | 1.557 | 0.004638 | 0.02256 | `CLOSED_BY_WS01` |
| 7 | 3577 | `root/y0/y1/y1/y1` | 0.0008171 | 1.73 | 0.005288 | 0.02314 | `CLOSED_BY_WS01` |
| 8 | 3577 | `root/y1/y0/y0/y0` | 0.0008568 | 1.909 | 0.005975 | 0.02374 | `CLOSED_BY_WS01` |
| 9 | 3577 | `root/y1/y0/y0/y1` | 0.0008947 | 2.094 | 0.006696 | 0.02438 | `CLOSED_BY_WS01` |
| 10 | 3577 | `root/y1/y0/y1/y0` | 0.0009309 | 2.285 | 0.007445 | 0.02508 | `CLOSED_BY_WS01` |
| 11 | 3577 | `root/y1/y0/y1/y1` | 0.0009655 | 2.48 | 0.008225 | 0.02582 | `CLOSED_BY_WS01` |
| 12 | 3577 | `root/y1/y1/y0/y0` | 0.0009984 | 2.681 | 0.009109 | 0.02655 | `CLOSED_BY_WS01` |
| 13 | 3577 | `root/y1/y1/y0/y1` | 0.00103 | 2.885 | 0.01006 | 0.02728 | `CLOSED_BY_WS01` |
| 14 | 3577 | `root/y1/y1/y1/y0` | 0.00106 | 3.094 | 0.01105 | 0.02803 | `CLOSED_BY_WS01` |
| 15 | 3577 | `root/y1/y1/y1/y1` | 0.001088 | 3.306 | 0.01205 | 0.02883 | `CLOSED_BY_WS01` |

## Interpretation

On owner family 4488:3 of CELL-02-03's n=15 wall-separation failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~77.50x miss to a ~4.17x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 16/16 boxes close by the conservative 2R bound and 0 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics.

## Next dependency

B2 (signed wall endpoints) and B3 (owner-family local model on next-worst groups, e.g., 4571:2 with 11 failures and 2484:4 with 16 failures) are unblocked. Track B1 viable. Recommended: replicate this experiment on owner family 4571:2 next, then run a Rust-level interval certification of the moving-frame collar that emits z* via interval Newton on F(z)=0.

## Source provenance

- Boundary-slice source: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03/EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03_RESULTS.json`
- SHA-256 status: `PASS`
- Wall-separation-target packet: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01/EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01_RESULTS.json`
- SHA-256 status: `PASS`
- Math reference: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md` (lines 122-146)

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets) were not
attempted. The new build script
`build_n15_cell0203_centerline_sign_model.py` writes only this experiment
folder. No prior script or artifact was modified.
