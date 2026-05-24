# EXP-MATH-ERDOS-1038-WS01-APPLICABILITY-PROBE-20260514-01

## Scope

Diagnostic, not a proof. Cross-problem applicability probe for the WS-01
(Branch-Centered Moving-Frame Collar) wall-separation rewrite, generalized
from Erdos #114's CELL-02-03 owner-family-4488:3 reference run.

## Verdict

`WS01_APPLICABLE`

## Headline numbers

- Failure boxes processed: `12`
- Branch-point validation rate: `1.0`
- Closure rate at conservative 2R: `1.0`
- Closure rate at tight R (cumulative): `1.0`
- Structural-difference rate: `0.0`
- Closed at 2R: `12`
- Closed at tight R only: `0`
- Still failing: `0`
- Branch point not found: `0`
- Branch point not well-defined: `0`
- Median required factor (RHS/LHS) before WS-01: `0.9999002056184547`
- Median required factor (RHS/LHS) after WS-01 at 2R: `2503.7242366608934`

## Schema validation

All rows passed schema validation.

## What WS-01 does

The original wall test bounds `sup_r |F(0,r)|` by the raw absolute interval
hull of `|F|` on the center strip. That hull carries the linear-in-R
coefficient `|F_t(box-midpoint)|`, which can be O(1).

WS-01 picks a validated zero point `z*` on the F=0 curve in each box and
uses the tangent direction `t` to that curve as the new coordinate. By
construction `F(z*) = 0` and `F_t(z*) = 0`, so the linear-in-R term
vanishes:

```
new_LHS = 0.5 |F_tt| R^2 + (T3/6) R^3 + R_n
```

Branch-point validation per box uses two interval-arithmetic facts:

1. `f_interval` straddles zero (IVT => F has a zero in the box).
2. `center_fn_abs_lower > 0` (grad F nonvanishing => zero curve is a smooth
   1-manifold => `t` is well-defined).

`F_t(z*) = 0` then holds by construction, not by any numerical estimate.

`R` is taken as `tangent_radius_upper` from the source data. The probe
reports the conservative `2R` figure as the canonical pass condition
(worst case z* near a box corner) and the tight `R` figure as a sensitivity
diagnostic.

## Per-failure outcomes

| # | src_idx | owner | R | wall_LHS_before | wall_LHS_after_2R | wall_rhs | outcome |
|---|---------|-------|---|-----------------|-------------------|----------|---------|
| 0 | 0 | `w=0.100/d=0.100` | 0.0001 | 0.001028 | 4.216e-07 | 0.001027 | `CLOSED_BY_WS01` |
| 1 | 1 | `w=0.100/d=0.150` | 0.0001 | 0.001024 | 4.164e-07 | 0.001024 | `CLOSED_BY_WS01` |
| 2 | 2 | `w=0.100/d=0.200` | 0.0001 | 0.00102 | 4.113e-07 | 0.00102 | `CLOSED_BY_WS01` |
| 3 | 3 | `w=0.130/d=0.100` | 0.0001 | 0.00102 | 4.109e-07 | 0.00102 | `CLOSED_BY_WS01` |
| 4 | 4 | `w=0.130/d=0.150` | 0.0001 | 0.001012 | 4.007e-07 | 0.001012 | `CLOSED_BY_WS01` |
| 5 | 5 | `w=0.130/d=0.200` | 0.0001 | 0.001005 | 3.907e-07 | 0.001005 | `CLOSED_BY_WS01` |
| 6 | 6 | `w=0.160/d=0.100` | 0.0001 | 0.00102 | 4.109e-07 | 0.00102 | `CLOSED_BY_WS01` |
| 7 | 7 | `w=0.160/d=0.150` | 0.0001 | 0.001012 | 4.007e-07 | 0.001012 | `CLOSED_BY_WS01` |
| 8 | 8 | `w=0.160/d=0.200` | 0.0001 | 0.001005 | 3.907e-07 | 0.001005 | `CLOSED_BY_WS01` |
| 9 | 9 | `w=0.190/d=0.100` | 0.0001 | 0.00102 | 4.109e-07 | 0.00102 | `CLOSED_BY_WS01` |
| 10 | 10 | `w=0.190/d=0.150` | 0.0001 | 0.001012 | 4.007e-07 | 0.001012 | `CLOSED_BY_WS01` |
| 11 | 11 | `w=0.190/d=0.200` | 0.0001 | 0.001005 | 3.907e-07 | 0.001005 | `CLOSED_BY_WS01` |

## Interpretation

On 12 wall-separation failure boxes, WS-01 closes 12 at the conservative 2R bound and an additional 0 at the tight-R bound. Branch-point validation rate 100.00%, closure rate at tight R 100.00%. The architecture appears applicable to this problem's boundary-slice failures.

## Next action

Proceed with WS-01 for this problem. Run owner-family-by-owner-family probes to confirm uniform closure across the full wall-failure set, then attempt Rust-level interval certification of the moving-frame collar.

## Source provenance

- Input JSON: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos1038/proof_path/_staging_1038_boundary_slice_2026-05-14.json`
- SHA-256: `13b3168dc1689d7d5718341e6e0b36cffb0cbf5f39890d9fb9ccf1be30055c2e`

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets,
reference implementation script) were not attempted. The probe writes
only this experiment folder.
