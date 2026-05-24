# EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01

## Scope

Internal experiment artifact only. Track B1-replicate of the EHP114 bridge
program. Not a proof of Erdos #114, not an n=15 certificate, and not a
CELL-02-03 closure. Only a diagnostic of WS-01-CENTER-STRIP-CANCELLATION (the
Branch-Centered Moving-Frame Collar rewrite) on owner family `2484:4`'s
16 wall-separation failures. Python numerical
demonstration at the same claim level as the 4488:3 run; not Rust interval-
certified.

## Status

`CENTERLINE_SIGN_MODEL_FULL_PASS_2484_4`

## Headline numbers

- Failures in family: `16`
- Closed by WS-01 (conservative 2R bound): `16`
- Closed by WS-01 (tight R only): `0`
- Still failing: `0`
- Branch point not found: `0`
- Median required factor (RHS/LHS) before WS-01: `0.028722739683563148`
- Median required factor (RHS/LHS) after WS-01 at 2R: `8.866047942748601`

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

## Per-failure outcomes (owner family 2484:4)

| # | src_idx | split_path | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|------------|---|------------------|--------------------|----------|---------|
| 0 | 275 | `root/y0/y0/y0/y0` | 0.0006029 | 0.3721 | 0.0008745 | 0.01086 | `CLOSED_BY_WS01` |
| 1 | 275 | `root/y0/y0/y0/y1` | 0.0005358 | 0.3297 | 0.0007417 | 0.01137 | `CLOSED_BY_WS01` |
| 2 | 275 | `root/y0/y0/y1/y0` | 0.0004688 | 0.2858 | 0.0006304 | 0.01187 | `CLOSED_BY_WS01` |
| 3 | 275 | `root/y0/y0/y1/y1` | 0.0004501 | 0.2554 | 0.0006174 | 0.01227 | `CLOSED_BY_WS01` |
| 4 | 275 | `root/y0/y1/y0/y0` | 0.0005149 | 0.2511 | 0.0007667 | 0.01253 | `CLOSED_BY_WS01` |
| 5 | 275 | `root/y0/y1/y0/y1` | 0.0005777 | 0.2472 | 0.0009415 | 0.01275 | `CLOSED_BY_WS01` |
| 6 | 275 | `root/y0/y1/y1/y0` | 0.0006383 | 0.3183 | 0.001142 | 0.01295 | `CLOSED_BY_WS01` |
| 7 | 275 | `root/y0/y1/y1/y1` | 0.0006967 | 0.3925 | 0.00137 | 0.01312 | `CLOSED_BY_WS01` |
| 8 | 275 | `root/y1/y0/y0/y0` | 0.0007529 | 0.4694 | 0.001628 | 0.01327 | `CLOSED_BY_WS01` |
| 9 | 275 | `root/y1/y0/y0/y1` | 0.0008069 | 0.5489 | 0.001913 | 0.0134 | `CLOSED_BY_WS01` |
| 10 | 275 | `root/y1/y0/y1/y0` | 0.0008588 | 0.6306 | 0.00222 | 0.01349 | `CLOSED_BY_WS01` |
| 11 | 275 | `root/y1/y0/y1/y1` | 0.0009084 | 0.7142 | 0.002568 | 0.01356 | `CLOSED_BY_WS01` |
| 12 | 275 | `root/y1/y1/y0/y0` | 0.000956 | 0.7994 | 0.002952 | 0.01361 | `CLOSED_BY_WS01` |
| 13 | 275 | `root/y1/y1/y0/y1` | 0.001001 | 0.8858 | 0.003357 | 0.01365 | `CLOSED_BY_WS01` |
| 14 | 275 | `root/y1/y1/y1/y0` | 0.001045 | 0.9731 | 0.003778 | 0.0137 | `CLOSED_BY_WS01` |
| 15 | 275 | `root/y1/y1/y1/y1` | 0.001086 | 1.061 | 0.004211 | 0.01377 | `CLOSED_BY_WS01` |

## Interpretation

On owner family 2484:4 of CELL-02-03's n=15 wall-separation failures (16 boxes), the WS-01 Branch-Centered Moving-Frame Collar rewrite reduces the median wall-LHS-over-RHS deficit from ~34.82x miss to a ~8.87x clearance (factor of safety) when the worst-case z*-off-center radius 2R is used. 16/16 boxes close by the conservative 2R bound and 0 more close only with the tight-R bound. 0 boxes have no interval-validated branch point. Per-box validation: F_iv contains 0 (branch exists by IVT) and center |F_n| lower bound is large positive (grad F nonvanishing, so the zero curve is smooth and t(u) is well-defined). F_t(z*)=0 holds by construction, not by numerics. This is a Python numerical demonstration at the same claim level as the 4488:3 run; NOT a Rust interval-certified result.

## Next dependency

Track B1-replicate confirmed full pass on owner family 2484:4. Combined with the 4488:3 result (16/16 closed), Track B1 viable across the worst-three owner families. Recommended next: replicate on the remaining matched-count families (3104:0, 3546:0, 3982:0, 4404:0 with 16 each), then move toward Rust interval-certification of the moving-frame collar that emits z* via interval Newton on F(z)=0.

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
