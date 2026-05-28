# Round 23 — Next Local Gate

After Round 23's Mode 1 literature survey, the global_reduction route has three putative morphisms with explicit maps and falsifiable tests, none of which are verified. The route status advances from "named open thread" to "three candidate morphisms with explicit maps," but remains SUMMIT-LEVEL OPEN.

Route state snapshot:
- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL, Chebyshev-rescaled working basis (unchanged from Round 20 pivot)
- `global_reduction`: SUMMIT-LEVEL OPEN — three putative morphisms named (A: hyperelliptic Jacobian, B: NPS theorem, C: KKT optimality), none tested
- `kkt_strict_slack`: SUMMIT-LEVEL OPEN (unchanged; C3 f64 confound still open)
- Six receipts still absent
- FT-02A methodology gap surfaced: corrected test requires complex b-period contour integration

---

## Recommended local actions (priority order)

### Action 1 — Run FT-05A (NPS component-count sweep) — INLINE, CHEAP

**What:** For monic real-rooted polynomials of degrees `n ∈ {10, 20, 50, 100}` with equally spaced roots in `[-1, 1]`, compute:
1. The number of connected components of `{|f_n(x)| < 1} ∩ [-1, 1]`
2. The NPS doubling exponent `β*(D, log|f_n|)` numerically
3. Whether component count grows with n or stays bounded
4. Comparison of `m({|f_n| < 1})` to `c_NPS / log(component_count_n)`

**Why first:** No private receipts, no interval arithmetic, no complex contours. Pure f64 computation with standard quadrature. Can be run in an afternoon. Directly tests the key assumption of Morphism B (bounded gap-component count ≤ 25 as n grows). The result determines whether the NPS lower bound is constant or degrades logarithmically.

**Accept criterion:** Component count ≤ 25 for all tested degrees AND measured `m` ≥ `0.5 c_NPS / log(50)`.
**Reject criterion:** Component count grows with n OR measured `m` falls below `c_NPS / log(2n)` for some n.

**Route consequence if pass:** Morphism B advances from PUTATIVE to TOY-SCALE-POSITIVE. Component-count assumption has toy-scale evidence.
**Route consequence if fail:** Morphism B's constant-lower-bound form is falsified; log(n)-degrading form may still hold.

---

### Action 2 — Run FT-06A (KKT finite-difference Jacobian) — INLINE, CHEAP

**What:** For a degree-10 polynomial with 10 equally spaced roots in `[-1, 1]`:
1. Perturb each root by ε = 0.01
2. Measure the change in `m({|f| < 1})` for each perturbation
3. Compute finite-difference Jacobian `∂m/∂a_i ≈ (m(a_i+ε) - m(a_i-ε)) / (2ε)`
4. Compare to the analytic gap-period matrix entries from the co-area formula

**Why second:** Also cheap and inline; no private receipts, no complex contours. Directly tests whether the gap-period matrix entries coincide with the KKT Jacobian of the measure functional. This is the foundational verification for Morphism C.

**Accept criterion:** Finite-difference Jacobian agrees with analytic co-area formula to 1% relative error.
**Reject criterion:** Disagreement above 1% (implies co-area formula not applicable or singularity in measure gradient at this configuration).

**Route consequence if pass:** Morphism C advances from PUTATIVE to TOY-SCALE-POSITIVE. KKT-as-Jacobian identification has toy-scale evidence.
**Route consequence if fail:** Morphism C's identification is broken at this configuration; investigate singularity.

---

### Action 3 — Implement the corrected FT-02A (Morphism A) — SUBSTANTIAL

The FT-02A test as delivered by PC tests `M_a` (the real a-period matrix) for symmetry and positive-definiteness. This is the wrong test — `M_a` is generically not symmetric; symmetry belongs to the full symplectic period matrix `Omega`.

**Corrected FT-02A protocol:**

1. Fix toy genus-2 hyperelliptic curve: `y² = ∏_{j=1}^{2} (x-a_j)(x-b_j)` for some gap endpoints, e.g., `a_1 = -1, b_1 = -0.3, a_2 = 0.2, b_2 = 0.8`.

2. Compute 2×2 **a-period matrix** `M_a`:
   ```
   M_a[i,j] = integral_{[a_j, b_j]} x^{i-1} / sqrt(|(x-a_1)(x-b_1)(x-a_2)(x-b_2)|) dx
   ```
   Standard Gaussian quadrature with endpoint-safe substitution. Real-valued.

3. Compute 2×2 **b-period matrix** `M_b`:
   The b-cycles for a real hyperelliptic curve go from `b_j` to `a_{j+1}` along the real axis, crossing to the other sheet. For a real curve, `M_b` is purely imaginary.
   Concretely: integrate `x^{i-1} / y` along the real segment `[b_1, a_2]` on both sheets (the two-sheeted cover); the result is `2i · integral_{[b_1, a_2]} x^{i-1} / sqrt(|(x-a_1)(x-b_1)(x-a_2)(x-b_2)|) dx`.
   Reference implementation: Molin–Neurohr (arXiv:1707.07249), or the SageMath `hyperelliptic_curve.period_matrix()` function.

4. Form the symplectic period matrix: `Omega = M_a^{-1} M_b`.

5. Verify:
   - `Omega - Omega^T` has entries below 1e-10 in absolute value (symmetry)
   - `Im(Omega)` is positive definite: `det(Im(Omega)) > 0` and `Im(Omega)[0,0] > 0`

**Why this is the correct test:** The Riemann bilinear relations are proved theorems for the symplectic period matrix, not for `M_a` alone. The morphism A claim is "the #1038 gap-period matrix identifies with `Omega`," and `Omega` is computed from both a-periods and b-periods.

**Effort:** Moderate. Step 3 requires setting up complex contour integration or using an existing period-matrix library. The toy g=2 infrastructure verified above (cond(M_a) = 2.625, det = 13.2) confirms the a-period computation is non-degenerate. The b-period computation adds complexity but is well-studied for g=2.

---

### Action 4 — Future PC round: formalize A+C composability — DEFERRED

The composability observation from PC's PUTATIVE_MORPHISMS_TO_1038.md is the most intellectually interesting direction to formalize: if Morphism A (M = Omega) and Morphism C (M = KKT Jacobian) both hold simultaneously, then KKT stationarity ⇔ Omega ∈ H_{24} (the extremal polynomial's optimality conditions are equivalent to the period matrix lying in the Siegel modular variety with gap-endpoint constraints).

This is worth a future PC round to:
- Determine the conditions under which both identifications hold simultaneously
- Check whether the composability is coherent (the two "M =" identifications must refer to the same matrix with the same normalization)
- Identify what geometric structure on the Siegel modular variety corresponds to the #1038 extremal measure condition
- Assess whether the resulting global reduction is stronger than the individual morphisms

**Why deferred:** Actions 1–3 should complete first to determine whether Morphisms B and C have any toy-scale evidence. If FT-06A fails, the composability discussion is premature. If FT-06A passes, a well-grounded formalization request to PC becomes feasible.

---

## What stays absent (unchanged from prior rounds)

- Six receipts (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS) — still absent; none of the three new morphisms circumvent this blocker on the dependent-Vieta consumer path
- G3 (interval-arithmetic re-implementation) — still local-only; required for any certified seed on the canonical hyperelliptic basis route
- No theorem advance is implied by completing Actions 1–3 above; passing all three tests would advance the global_reduction route from "three putative morphisms" to "three putative morphisms with toy-scale evidence," raising confidence but not closing any summit-level gate

---

## Round 24 dispatch question

Three reasonable Round 24 PC-shaped scopes:

**Option α — FT-05A + FT-06A only (if not run locally).** Both tests are cheap and well-specified. If local bandwidth is limited, PC can execute them. Outcome directly gates the composability formalization.

**Option β — Corrected FT-02A (b-period contour integration).** If the local agent does not have a period-matrix library set up, PC can implement the Molin–Neurohr algorithm for the g=2 case. More involved than FT-05A/FT-06A, but the methodology is published and PC has demonstrated capability with period-matrix computation from prior rounds.

**Option γ — A+C composability formalization.** After Actions 1–3 are complete (or partially complete), a fresh PC round on the composability question. This is the highest-leverage research direction but requires the prior tests as prerequisites.

Option α is the most conservative continuation. Option γ is the most strategically interesting if Actions 1–3 are done locally before the next dispatch.
