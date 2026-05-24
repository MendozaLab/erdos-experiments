# EXP-MATH-EHP114-N15-CELL-02-03-WS01-RESIDUAL-DEEP-DIVE-20260509-01

**Status:** `INTERVAL_NEWTON_DOES_NOT_FIX`

**Claim ceiling:** Internal diagnostic only. One-box decomposition of the WS-01 wall test on CELL-02-03 owner family 4571:2 worst-failing box. Python numerical demonstration. The R-reduction-after-interval-Newton number is an UNCERTIFIED estimate based on geometry, not a certified bound.

## Selected box

- Owner key: `4571:2`
- src_idx: `3852`
- split_path: `root/y0/y0/y1/y0`
- x_interval: `[0.8957386, 0.8965909]`
- y_interval: `[-0.1466974, -0.1424893]`
- Rationale: smallest `required_factor_after_2R = 4.43e-03` among 6 STILL_FAILING_AFTER_WS01 boxes in family 4571:2. Wall LHS misses the threshold by ~226× at the conservative 2R bound.

## Wall-test inequality being checked (WS-01 rewrite)

```
LHS = 0.5 * |F_tt|_upper * R^2 + (T_3/6) * R^3 + R_n
RHS = S * lower(|F_n|)
PASS if LHS < RHS
```

## Per-term decomposition (evaluated at 2R, the worst-case z* placement)

| Term | Formula | Value | Share of LHS |
|---|---|---|---|
| Quadratic | `0.5 * |F_tt|_upper * (2R)^2` | `1.408e-02` | 77.46% |
| Cubic | `(T_3/6) * (2R)^3` | `3.258e-04` | 1.79% |
| Remainder | `R_n` | `3.772e-03` | 20.75% |
| **LHS total** | | **`1.818e-02`** | 100% |
| RHS | `S * lower(|F_n|)` | `8.063e-05` | |
| **LHS / RHS** | | **`225.5×`** | |

Inputs: `R = 1.822e-03`, `|F_tt|_upper = 2120.3`, `T_3_upper = 4.038e+04`, `R_n = 3.772e-03`.

## Per-term decomposition at tight R (best-case z* on box center)

| Term | Value | Share |
|---|---|---|
| Quadratic | `3.521e-03` | 48.0% |
| Cubic | `4.07e-05` | 0.6% |
| Remainder | `3.772e-03` | 51.4% |
| **LHS total** | **`7.334e-03`** | |
| **LHS / RHS** | **`91.0×`** | |

Even at the tight-R bound, the remainder term alone (`3.77e-03`) is 47× above RHS. Quadratic and remainder are roughly co-equal contributors at R; at 2R the quadratic dominates because it scales as R² while remainder is R-independent.

## Comparison with a closing box (4488:3)

Closest-R_n closing box from owner family `4488:3`:

- src_idx: `3577`, split_path: `root/y1/y0/y0/y1`
- Outcome: `CLOSED_BY_WS01`

| Term | Failing (4571:2) | Closing (4488:3, comparable R_n) |
|---|---|---|
| Quad (2R) | `1.41e-02` | comparable |
| Cubic (2R) | `3.26e-04` | comparable |
| R_n | `3.77e-03` | similar (~1.4–3.3e-03 typical) |
| **LHS** | **`1.82e-02`** | similar magnitude |
| **RHS** | **`8.06e-05`** | **`~2e-02`** |
| LHS / RHS | 225× over | 0.09–0.25× under |

**The decisive difference is the wall RHS, not the LHS.** Across all 4571:2 failing boxes, `S*lower(|F_n|)` ranges 8.06e-05 to 8.15e-03 (median 2.62e-03). Across all 4488:3 closing boxes, it ranges 1.93e-02 to 2.88e-02 (median 2.34e-02). That's a ~10–240× systematic gap in wall threshold, set by cell geometry near the |F| zero-curve.

## What does interval-Newton certified z* actually do?

Interval-Newton tightens R from `tangent_radius_upper` (worst-case off-center placement of z* in the box) to a smaller "true off-center radius." Estimated effect:

- Quadratic term scales as R² → drops by `k²` if R drops by factor `k`.
- Cubic term scales as R³ → drops by `k³`.
- Remainder R_n is a Taylor remainder bound on the box; it does **not** depend on z*'s placement within the box. R-shrinkage does not touch it.

### Projection at an optimistic 10× R-reduction (UNCERTIFIED estimate)

| Term | Before (2R) | After (R/10) |
|---|---|---|
| Quadratic | `1.408e-02` | `1.408e-04` |
| Cubic | `3.258e-04` | `3.258e-07` |
| Remainder R_n | `3.772e-03` | `3.772e-03` (unchanged) |
| **LHS total** | `1.818e-02` | `3.913e-03` |
| RHS | `8.06e-05` | `8.06e-05` |
| LHS / RHS | 225× | **48.6×** |

Even with a 10× R-reduction, the box still misses by ~49×, dominated entirely by R_n.

### Theoretical floor at R = 0 (perfect interval-Newton)

LHS floor = R_n = `3.77e-03`, still 46.79× above wall_rhs = `8.06e-05`.
**R-shrinkage cannot close this box, no matter how tight.**

## All six STILL_FAILING boxes in 4571:2: R_n vs wall_rhs

| src_idx | split_path | R_n | wall_rhs | R_n / RHS | IN floor clears? |
|---|---|---|---|---|---|
| 3852 | root/y0/y0/y0/y0 | 3.139e-03 | 8.610e-04 | 3.65× | No |
| 3852 | root/y0/y0/y0/y1 | 3.475e-03 | 3.842e-04 | 9.05× | No |
| 3852 | root/y0/y0/y1/y0 | 3.772e-03 | 8.063e-05 | **46.79×** | No |
| 3852 | root/y1/y0/y0/y0 | 4.081e-03 | 5.001e-04 | 8.16× | No |
| 3852 | root/y1/y0/y0/y1 | 3.813e-03 | 1.475e-03 | 2.58× | No |
| 3852 | root/y1/y0/y1/y0 | 3.471e-03 | 2.622e-03 | 1.32× | No |

**For every single one, R_n already exceeds wall_rhs.** None can be saved by R-shrinkage alone.

## Verdict

**Dominant LHS term:** quadratic at 2R (77%); roughly tied with remainder at R (48% / 51%).

**Does interval-Newton fix this box?** **No.** The deficit is gated on R_n vs `S * lower(|F_n|)`, which is a cell-geometry property (the cell sits near the |F| zero-curve where lower(|F_n|) is small). Interval-Newton certified z* does not shift either of those quantities. Even at R = 0, LHS floors at R_n = 3.77e-3, which is 46.79× the wall threshold of 8.06e-5.

**Implications for the Modal Tier-2 Rust-port spend.** The Python-vs-Rust gap argument — that interval-Newton's certified z* would close the residual failures — does not hold for this box, and does not hold for any of the 6 STILL_FAILING boxes in 4571:2 (R_n / wall_rhs range: 1.32× to 46.79×). The Tier-2 spend may still be justified on the broader 4488:3 family or for marginal CLOSED_BY_WS01_TIGHT_R_ONLY edges, but it cannot be justified by reference to the 4571:2 residual.

## Next dependency

The 6 STILL_FAILING boxes in 4571:2 require one of:

1. **Cell subdivision** until each sub-cell's `lower(|F_n|)` rises enough to put `S*lower(|F_n|)` above its R_n. Most likely effective: `lower(|F_n|)` typically scales linearly with cell-width near the zero-curve, while R_n scales quadratically/cubically with cell-width. Halving the cell drops R_n by 4–8× and S*lower(|F_n|) by ~2×, so the net wall-test ratio improves by 2–4× per subdivision level. ~3–5 levels of subdivision should suffice for the worst box.
2. **Sharper R_n bound.** A smaller-degree Taylor remainder, or a change-of-frame that absorbs more of the curvature into the linearized term, could cut R_n. Direct analytic work, not a port.
3. **Analytic critical-point exclusion** to remove the cell entirely from the wall-test class. Most expensive but cleanest.

Cell subdivision is the recommended path of least resistance.

## Source artifacts (read-only inputs)

- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-OWNER-4571-2-20260509-01/.../_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260509-01/.../_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/proof_path/EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01/.../_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/build_n15_cell0203_centerline_sign_model.py`

SHA-256 receipts in `_RESULTS.json` under `source_sha_status`.
