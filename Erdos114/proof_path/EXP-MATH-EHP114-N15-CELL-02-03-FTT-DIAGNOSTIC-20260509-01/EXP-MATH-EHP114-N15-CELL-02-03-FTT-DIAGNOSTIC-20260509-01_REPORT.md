# EXP-MATH-EHP114-N15-CELL-02-03-FTT-DIAGNOSTIC-20260509-01

## Status

- FTT_DIAGNOSTIC_CRITICAL_POINT_HYPOTHESIS_REFUTED

## Claim Ceiling

Internal diagnostic only. Not a proof of Erdos #114, not a CELL-02-03 closure. Just the F_tt upper-bound distribution on 8 residual failing boxes vs closing-family baseline. Same Python-numerical claim level as the existing WS-01 work; not Rust interval-certified. Lower bound unavailable upstream.

## Hypothesis Under Test

8 residual WS-01 failures contain or are adjacent to critical points where |F_tt| vanishes, making moving-frame collar inherently insufficient.

## Data Availability

- `f_tt_abs_upper`: AVAILABLE — recorded per box by the upstream Rust binary `ehp114_n15_boundary_slice_cell_02_03` and consumed verbatim by the Python build scripts.
- `f_tt_abs_lower`: NOT AVAILABLE — the upstream Rust binary does not emit this field. Recovering it requires modifying the binary and re-running the boundary-slice. Out of scope for this Python diagnostic.

Consequence: this diagnostic can only test the WEAK form of the critical-point hypothesis (does the upper bound collapse on failing boxes?). The STRONG form (does |F_tt| straddle zero?) needs the lower bound and is reported as `null` in the output.

## Comparison: Failing vs Closing |F_tt| Upper Bounds

| Group | n | min | median | mean | max |
|---|---|---|---|---|---|
| Failing (residual at FULL_PASS-2R) | 8 | 329.0 | 1640.9 | 1557.9 | 2120.3 |
| Closing (FULL_PASS-2R) | 49 | 296.1 | 700.1 | 966.1 | 2902.6 |
| Ratio (failing / closing) | — | 1.11 | 2.34 | 1.61 | 0.73 |

## Per-Failing-Box |F_tt| Upper Bound

| owner_key | src_idx | split_path | F_tt upper | wall_lhs_after_2R | wall_RHS | req_factor (RHS/LHS) |
|---|---|---|---|---|---|---|
| 4571:2 | 3852 | root/y0/y0/y0/y0 | 2020.4 | 1.859799e-02 | 8.609650e-04 | 0.0463 |
| 4571:2 | 3852 | root/y0/y0/y0/y1 | 2084.2 | 1.856277e-02 | 3.841991e-04 | 0.0207 |
| 4571:2 | 3852 | root/y0/y0/y1/y0 | 2120.3 | 1.818041e-02 | 8.062605e-05 | 0.0044 |
| 4571:2 | 3852 | root/y1/y0/y0/y0 | 1717.2 | 1.047152e-02 | 5.000813e-04 | 0.0478 |
| 4571:2 | 3852 | root/y1/y0/y0/y1 | 1564.6 | 8.809893e-03 | 1.475039e-03 | 0.1674 |
| 4571:2 | 3852 | root/y1/y0/y1/y0 | 1399.5 | 7.222397e-03 | 2.622305e-03 | 0.3631 |
| 4571:2 | 3852 | root/y1/y0/y1/y1 | 1227.6 | 5.761015e-03 | 3.868040e-03 | 0.6714 |
| 3982:0 | 2476 | root/y1/y1/y1/y1 | 329.0 | 2.129087e-03 | 7.354214e-04 | 0.3454 |

## Verdict

REFUTED (weak form). The 8 residual failing boxes have a median |F_tt| upper bound of 1640.9, vs 700.1 for closing-family boxes (ratio 2.34x). Failing boxes do NOT have collapsing |F_tt| upper bounds — if anything they are larger. The wall test is failing for some other reason (small wall RHS / S * lower(|F_n|), large T3 cubic remainder, or large normal_remainder R_n). Modal Tier-2 Rust port + interval-Newton recovery story remains plausible. STRONG-FORM REFUTATION (proving |F_tt| bounded uniformly BELOW on the failing boxes) requires adding f_tt_abs_lower to the upstream Rust binary — recommend folding that into Tier-2 as a small free-rider task.

## Implications for Modal Tier-2 Spend

- **Spend justified:** conditional
- **Next dependency:** Tier-2 Rust port with f_tt_abs_lower emitter folded in. Lower bound proves strong-form refutation; without it, the weak-form refutation here only shows the upper bound is not collapsing, which is necessary but not sufficient.

## Source Paths

- `owner_4571_2`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01_RESULTS.json`
- `owner_3982_0`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-3982-0-20260509-01_RESULTS.json`
- `owner_4488_3_aka_default_run`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01_RESULTS.json`
- `owner_2484_4`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2484-4-20260509-01_RESULTS.json`
- `owner_2932_1`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-2932-1-20260509-01_RESULTS.json`
- `owner_4404_0`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4404-0-20260509-01_RESULTS.json`

