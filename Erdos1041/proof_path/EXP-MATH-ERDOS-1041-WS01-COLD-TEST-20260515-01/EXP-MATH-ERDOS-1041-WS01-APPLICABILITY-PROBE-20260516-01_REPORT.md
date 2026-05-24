# EXP-MATH-ERDOS-1041-WS01-APPLICABILITY-PROBE-20260516-01

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
- Median required factor (RHS/LHS) before WS-01: `0.13243962915477572`
- Median required factor (RHS/LHS) after WS-01 at 2R: `343.7622214399209`

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
| 0 | 0 | `theta=0.100/krot=0.000` | 0.0001 | 0.01189 | 3.617e-06 | 0.001303 | `CLOSED_BY_WS01` |
| 1 | 1 | `theta=0.100/krot=0.050` | 0.0001 | 0.01124 | 3.354e-06 | 0.001234 | `CLOSED_BY_WS01` |
| 2 | 2 | `theta=0.100/krot=0.100` | 0.0001 | 0.0106 | 3.267e-06 | 0.001211 | `CLOSED_BY_WS01` |
| 3 | 3 | `theta=0.500/krot=0.000` | 0.0001 | 0.01197 | 5.853e-06 | 0.001817 | `CLOSED_BY_WS01` |
| 4 | 4 | `theta=0.500/krot=0.050` | 0.0001 | 0.01234 | 6.399e-06 | 0.001927 | `CLOSED_BY_WS01` |
| 5 | 5 | `theta=0.500/krot=0.100` | 0.0001 | 0.01232 | 6.641e-06 | 0.001975 | `CLOSED_BY_WS01` |
| 6 | 6 | `theta=0.950/krot=0.000` | 0.0001 | 0.0099 | 4.213e-06 | 0.001454 | `CLOSED_BY_WS01` |
| 7 | 7 | `theta=0.950/krot=0.050` | 0.0001 | 0.00935 | 3.719e-06 | 0.001331 | `CLOSED_BY_WS01` |
| 8 | 8 | `theta=0.950/krot=0.100` | 0.0001 | 0.00892 | 3.404e-06 | 0.00125 | `CLOSED_BY_WS01` |
| 9 | 9 | `theta=1.400/krot=0.000` | 0.0001 | 0.01262 | 4.295e-06 | 0.001471 | `CLOSED_BY_WS01` |
| 10 | 10 | `theta=1.400/krot=0.050` | 0.0001 | 0.01362 | 4.957e-06 | 0.001623 | `CLOSED_BY_WS01` |
| 11 | 11 | `theta=1.400/krot=0.100` | 0.0001 | 0.01424 | 5.668e-06 | 0.001777 | `CLOSED_BY_WS01` |

## Interpretation

On 12 wall-separation failure boxes, WS-01 closes 12 at the conservative 2R bound and an additional 0 at the tight-R bound. Branch-point validation rate 100.00%, closure rate at tight R 100.00%. The architecture appears applicable to this problem's boundary-slice failures.

## Next action

Proceed with WS-01 for this problem. Run owner-family-by-owner-family probes to confirm uniform closure across the full wall-failure set, then attempt Rust-level interval certification of the moving-frame collar.

## Source provenance

- Input JSON: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos1041/proof_path/INPUT_boundary_slice_1041_2026-05-15.json`
- SHA-256: `fc2fc8e04d0937b7c167fe906caa17910a643c4c9b8754c7699ed59f0218c709`

## Safety

This artifact is review-only. Forbidden writes (D1, curated morphisms,
public pages, proof registries, existing EHP114 finite packets,
reference implementation script) were not attempted. The probe writes
only this experiment folder.
