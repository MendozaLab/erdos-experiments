# Track A2 c_10 Grid-Evaluation Rigor Check

**Experiment ID:** `EXP-MATH-EHP114-TAO-C10-GRID-EVALUATION-RIGOR-CHECK-20260509-01`
**Date (UTC):** 2026-05-09
**Target script:** `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/build_tao_c10_disk_integral_lower_bound.py`
**Target artifact:** `proof_path/EXP-MATH-EHP114-TAO-DEFECT-PSI-COR-C10-LOWER-BOUND-20260509-01/`
**Concern (third-party reviewer, Perplexity + Gemini):** Endpoint sampling on a polar grid does NOT yield a certified interval lower bound; only true interval-form evaluation (Taylor models, monotonicity, or `mpmath.iv` operator chains over interval inputs) does.

## Verdict

**`RIGOROUS_INTERVAL_FORM`**

The script evaluates `g(r, θ) = (|z-1| + |z| - |2z-1|) / |z-1|` on `mpmath.iv` interval objects spanning the full cell, not on point samples. The reported lower bound `c_10 ≥ 0.006900495524140533` is a certified rigorous lower bound, modulo the standard correctness assumption on the `mpmath.iv` library itself.

## Method as implemented

For each of the 200 × 720 = 144,000 polar cells `[r_lo, r_hi] × [θ_lo, θ_hi]`, the script does the following (file: `build_tao_c10_disk_integral_lower_bound.py`):

```
# Cell loop — line 178-196
for i in range(n_r):
    r_iv  = iv.mpf([r_lo, r_hi])              # line 180  — interval over WHOLE radial subinterval
    dr_iv = iv.mpf([r_hi - r_lo, r_hi - r_lo]) # line 181  — exact width
    for j in range(n_theta):
        t_iv  = iv.mpf([t_lo, t_hi])           # line 185  — interval over WHOLE angular subinterval
        dt_iv = iv.mpf([t_hi - t_lo, t_hi - t_lo])
        g_iv  = integrand_interval(r_iv, t_iv) # line 188  — interval evaluation of g on the cell
        g_lower = mpf(g_iv.a)                  # line 191  — LOWER endpoint of the enclosure
        cell_lower = iv.mpf([g_lower, g_lower]) * dr_iv * dt_iv   # line 192
        total_lo += cell_lower                 # line 193
```

The `integrand_interval` function (lines 117-151) is composed entirely of `mpmath.iv` operators on interval inputs:

```
# integrand_interval body — line 117-151
cos_t   = iv.cos(theta_iv)                                          # line 127
abs_z   = r_iv                                                      # line 130
sq_zm1  = r_iv * r_iv - 2 * r_iv * cos_t + 1                        # line 133
sq_zm1  = iv.mpf([max(sq_zm1.a, (20/21)^2), sq_zm1.b])               # line 138 — P2 tightening
abs_zm1 = iv.sqrt(sq_zm1)                                           # line 139
sq_2zm1 = 4 * r_iv * r_iv - 4 * r_iv * cos_t + 1                    # line 142
sq_2zm1 = iv.mpf([max(sq_2zm1.a, (19/21)^2), sq_2zm1.b])             # line 146 — P3 tightening
abs_2zm1 = iv.sqrt(sq_2zm1)                                         # line 147
g       = (abs_zm1 + abs_z - abs_2zm1) / abs_zm1                    # line 149-150
return g
```

Every operator (`+`, `-`, `*`, `/`, `iv.cos`, `iv.sqrt`) is the `mpmath.iv` interval-arithmetic version, which by contract returns an enclosure of the operator's range over the input intervals, with directed (outward) rounding at `dps=60`.

The endpoint-only sampler exists in the file but is the **non-rigorous** sanity routine, clearly labeled as such:

- `integrand_float(r, theta)` (lines 154-160) — point evaluation, mpf precision
- `nonrigorous_midpoint_estimate(...)` (lines 206-217) — midpoint Riemann sum
- `spot_check_nonnegativity(...)` (lines 223-240) — point-grid scan to confirm P1

These three are NOT in the lower-bound computation path. The function `rigorous_lower_bound` (lines 166-203) is the one whose output is recorded as `rigorous_lower_bound` in `_RESULTS.json`.

### One subtlety — the `lo_safe` clamps (lines 138, 146)

The script tightens the lower endpoint of `|z-1|^2` and `|2z-1|^2` to known analytic floors:

- `|z-1|^2 ≥ (20/21)^2` for z ∈ D(0, 1/21), since |z-1| ≥ 1 - |z| > 20/21 (precondition P2, file lines 22-26, 134-137).
- `|2z-1|^2 ≥ (19/21)^2` for z ∈ D(0, 1/21), since |2z-1| ≥ 1 - 2|z| > 19/21 (precondition P3, lines 143-146).

These are mathematical facts about the closed disk, not heuristics. Intersecting an interval-arithmetic enclosure with a proven-true lower bound is a sound *tightening* — it can only shrink the enclosure, never invalidate it. So this hybrid step does not break rigor; it patches the dependency-problem widening that would otherwise make `iv.sqrt` produce a needlessly wide (or nominally negative) interval at cells near θ = 0, r = 1/21.

## Rigor assessment

The bound `c_10 ≥ 0.006900495524140533` is a **certified rigorous lower bound**, assuming standard correctness of `mpmath.iv` (interval operators produce range enclosures with outward rounding). Specifically:

- For every cell, `g_iv.a` is a sound lower bound of `min g` over that cell — not a sample.
- Cell contribution = `g_iv.a · (r_hi - r_lo) · (θ_hi - θ_lo)`, computed in interval arithmetic with outward rounding (line 192). Note: when `g_iv.a` is negative, this contributes negatively to the running sum, which is conservative (it under-counts the cell, never over-counts).
- The summation `total_lo += cell_lower` (line 193) is interval-arithmetic accumulation, again with outward rounding.
- `total_lo.a` (line 203) is the lower endpoint of the running interval, hence a rigorous lower bound on `c_10`.

The 13,288 cells with negative `g_lower` (RESULTS.json field `negative_interval_min_cells`, ~9.2% of cells) are an honest reflection of the dependency-problem widening on a 200×720 grid — not a bug. The script lets them drag the bound down rather than clamping (lines 60-65, 367-369). Even with this honest accounting, the final bound stays positive and below the fine-grid midpoint estimate (0.006900 < 0.007069), as expected for a tight interval lower bound.

The result `c_10 ≥ 0.006900495524140533` is filing-grade rigorous. It can be quoted as a certified bound in proof artifacts.

## Recommended remediation

**None required.** The implementation matches the rigorous-interval-form pattern that Perplexity and Gemini specified. If at any future step a tighter bound is wanted (current gap to non-rigorous estimate is ~2.4%), refinement options in priority order:

1. **Increase grid resolution** (cheapest): bump `N_R_RIGOROUS = 200` and `N_THETA_RIGOROUS = 720` (lines 106-107) to 400 × 1440. Linear cost increase, quadratic narrowing of dependency widening. No code structure change.
2. **Cell-wise centered-form / mean-value-form evaluation**: instead of evaluating `g` on the whole cell interval, evaluate at the cell midpoint plus a Lipschitz / gradient-bound correction. Reduces dependency widening; needs a verified gradient bound on g.
3. **Adaptive subdivision** of the 13,288 negative-`g_lower` cells only, recursing until each yields `g_lower ≥ 0`. Targets the widening-affected cells without globally refining.

None of these are corrections — they would tighten an already-rigorous bound.

## Honest caveats

- Trust assumption: `mpmath.iv` correctly implements outward rounding for its operators at `dps=60`. This is standard mpmath behavior but is taken on faith — not independently proven in this artifact.
- The 9.2% negative-cell rate suggests the bound is not as tight as it could be; another implementation (e.g., centered-form or finer grid) might give 0.0070 instead of 0.0069. But that's a tightness question, not a rigor question. The 0.0069 figure is sound.
- The c_10 number itself is one numeric input. It is NOT a proof of Erdős #114 or an N₀ candidate — the script's `claim_ceiling` field already states this honestly.

## Files

- `EXP-MATH-EHP114-TAO-C10-GRID-EVALUATION-RIGOR-CHECK-20260509-01_REPORT.md` (this file)
- `EXP-MATH-EHP114-TAO-C10-GRID-EVALUATION-RIGOR-CHECK-20260509-01_RESULTS.json`
- `EXP-MATH-EHP114-TAO-C10-GRID-EVALUATION-RIGOR-CHECK-20260509-01_RESULTS.sha256`
