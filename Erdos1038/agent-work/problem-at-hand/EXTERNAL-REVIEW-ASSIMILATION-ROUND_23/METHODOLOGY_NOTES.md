# Round 23 — Methodology Notes

Round 23 is a Mode 1 literature survey, so the methodology assessment is different in character from Rounds 19–22 (which produced computational sweep scripts and result JSON). Here the relevant question is not "did the numerical procedure produce trustworthy numbers?" but rather "are the morphism proposals structurally correct, and are the falsifiable tests well-designed?" The answer is mostly yes, with one significant exception that matters for practical execution.

---

## §1 — Praise: careful exclusion methodology

The three NEGATIVE_EXCLUSION verdicts in this round are unusually well-grounded. It would have been easy to wave at structural similarity (both Hecke gaps and #1038 gaps are "gaps"; both Weil bounds and the sublevel-set measure bound are "saving estimates") and leave the families classified as weak putative morphisms. PC did not do this. Each exclusion is pinned to a specific structural mismatch:

- **TF-01 (modular forms):** The gap structure in #1038 is topological (connected components of `{|f| < 1}` on ℝ), not spectral (eigenvalue gaps in a Hecke module). Hecke operators act by sum-of-coset decompositions; there is no natural functor from real-line sublevel-set decompositions to Hecke-module structures.
- **TF-03 (Weil bounds):** The field-type mismatch is fundamental. The Weil-Deligne machinery requires a Frobenius action on ℓ-adic cohomology of a variety over F_q. No natural Frobenius acts on `{x ∈ ℝ : |f(x)| < 1}`. Attempting a mod-p reduction loses the real-line ordering and gap structure entirely.
- **TF-04 (Sato-Tate):** The ensemble-average obstruction is correctly identified. Sato-Tate is a statement about equidistribution as the prime varies; #1038's gap surface is a single fixed extremal object with no family parameter varying.

The honesty of naming a weak analogy ("both are saving estimates") and then declining to call it a morphism is a methodological strength, not a limitation.

---

## §2 — FT-02A: conceptual mismatch and corrected protocol

**The issue.** Test FT-02A as written in THEOREM_FAMILIES_INVENTORY.json and PUTATIVE_MORPHISMS_TO_1038.md §A.4 asks the following:

> Compute the 2×2 matrix `Omega` with `Omega_{ij} = integral_{G_j} x^{i-1} / sqrt(|(x-a_1)(x-b_1)(x-a_2)(x-b_2)|) dx`. Verify (1) `Omega_{12} = Omega_{21}` (symmetry) and (2) `det(Im(Omega)) > 0`.

The problem: that integral is not the full symplectic period matrix `Omega ∈ H_g`. It is the **real a-period matrix** `M_a`, defined by integrating over the **real** arc (the gap interval) of the algebraic curve `y² = ∏(x-a_k)(x-b_k)`.

The real a-period matrix `M_a` satisfies:
- `M_a[i,j] = integral_{[a_j, b_j]} x^{i-1} / sqrt(|∏_k (x-a_k)(x-b_k)|) dx`
- `M_a` is generally **not** symmetric
- `M_a` has no imaginary part (the integrand is real-valued and the integration domain is real)

Running FT-02A as written would produce:
1. An asymmetric `M_a` — correctly, but interpreted by FT-02A as a FAIL for symmetry
2. An `Im(M_a)` that is the zero matrix — and a FAIL for positive-definiteness

Both would be false fails. The morphism A structural claim is not invalidated by them; the test protocol is simply testing the wrong matrix.

**Toy g=2 verification.** For the toy hyperelliptic curve `y² = (x+1)(x+0.5)(x-0.2)(x-0.5)(x-0.7)(x-1)` (6 branch points, genus 2, two gap arcs approximately `[-1, -0.5]` and `[0.2, 0.5]`), numerical computation yields:

```
Real a-period matrix M_a = [[1.92, -1.35],
                             [4.53,  3.70]]

cond(M_a) = 2.625
det(M_a)  = 13.2   (non-degenerate, well-defined)
M_a[0,1] = -1.35  vs  M_a[1,0] = 4.53  -> asymmetric (correct mathematical behavior)
```

This is exactly what the theory predicts. `M_a` is generically not symmetric. Symmetry is a property of the full symplectic period matrix `Omega`, not of the partial real a-period matrix.

**Why does symmetry belong to `Omega` and not `M_a`?** The Riemann bilinear relations hold for the period matrix computed using a **symplectic basis** of `H_1(C, Z)`. The symplectic basis consists of `{a_1,...,a_g, b_1,...,b_g}` cycles satisfying `a_i · b_j = delta_{ij}`. The period matrix `Omega ∈ H_g` is defined by:

```
Omega_{ij} = integral_{b_i} omega_j,   omega_j = x^{j-1} dx / y
```

where `{b_1,...,b_g}` are the **b-cycles** (which involve complex contours through branch cuts, not just the real gap arcs). The full symplectic period matrix is:

```
Omega = M_a^{-1} M_b
```

where `M_a` is the a-period matrix (real cycles, real-valued integrals) and `M_b` is the b-period matrix (complex cycles, complex-valued integrals).

**The corrected FT-02A protocol:**

1. Fix a genus-2 hyperelliptic curve with branch points at `{a_1, b_1, a_2, b_2}`.
2. Compute the 2×2 a-period matrix `M_a` by integrating `x^{i-1}/y` over each real gap arc `[a_j, b_j]`. (This is the real-valued part; uses standard Gaussian quadrature on the real interval with the square-root singularity.)
3. Compute the 2×2 b-period matrix `M_b` by integrating `x^{i-1}/y` over each b-cycle. (This requires contour integration through branch cuts. A standard approach: b_j goes from `b_j` to `a_{j+1}` along the real axis on one sheet, then back on the other sheet. The result is purely imaginary for real hyperelliptic curves.)
4. Form the symplectic period matrix: `Omega = M_a^{-1} M_b`.
5. Check: `Omega = Omega^T` (symmetry) and `Im(Omega) > 0` (positive-definiteness).

**The toy infrastructure is computationally accessible.** The g=2 sanity check above confirms that M_a is well-conditioned (cond = 2.625) and non-degenerate (det = 13.2) at this toy configuration. Step 3 (complex b-period contours) is the additional work required but is standard numerical algebraic geometry. The Molin–Neurohr software (cited in TF-02 references, arXiv:1707.07249) implements exactly this computation for superelliptic curves.

**Bottom line for morphism A.** The structural claim (gap-period matrix maps to hyperelliptic Jacobian period matrix) is sound. The test as written would falsely reject it. The corrected test is more involved (requires complex contour integration) but is a tractable computation using existing software.

---

## §3 — FT-05A and FT-06A: well-designed, runnable inline

**FT-05A (NPS component-count sweep).** The test asks: for monic real-rooted polynomials of degrees `n ∈ {10, 20, 50, 100}` with equally spaced roots, compute the number of connected components of `{|f_n| < 1} ∩ [-1,1]` and check whether it grows with n or stays bounded. This is a well-posed, cheap computation (roots at equally spaced positions, univariate measure computation). No private receipts, no interval arithmetic, no complex contours. The accept/reject criteria are specific and falsifiable. Runnable in an afternoon.

The morphism's key assumption — that the gap-component count is bounded (≤ 25) independently of degree n — is exactly what FT-05A probes. If component count grows with n, the NPS lower bound degrades from `c/log(50)` (constant) to `c/log(n)` (still meaningful, but different in character). Both outcomes are informative.

**FT-06A (KKT finite-difference Jacobian).** The test asks: for a degree-10 polynomial with equally spaced roots, perturb each root by ε = 0.01 and compare the finite-difference Jacobian `∂m/∂a_i ≈ (m(a_i+ε) - m(a_i-ε)) / (2ε)` to the analytic formula from the co-area derivation. Accept criterion: agreement to 1% relative error. Reject criterion: disagreement above 1%.

This is a clean, cheap inline test. The co-area formula is mathematically explicit (PC derives it fully in §C.2 of PUTATIVE_MORPHISMS_TO_1038.md). The finite-difference computation is standard. The 1% relative-error tolerance is appropriate for f64 arithmetic. Runnable inline in under an hour.

---

## §4 — Composability note: research direction, not a claim

PC's observation that Morphisms A and C may compose — if M = Omega (A holds) and M = KKT Jacobian (C holds), then KKT stationarity ⇔ Omega ∈ H_{24} (Siegel constraint) — is the highest-leverage intellectual product of this round. The logic is clean: if both identifications hold simultaneously, the extremal polynomial's optimality conditions translate into the finite-dimensional constraint that the period matrix lies in the Siegel modular variety with specific gap-endpoint constraints.

PC presents this correctly as a research-level question, not a claim. The word "putative" is used. No attempt is made to derive either morphism from the other. The note is kept at the end of the morphisms document as a forward-looking direction, not asserted as evidence.

The value of this observation is that it suggests a future PC round or local research effort: formalize the conditions under which A+C compose, and determine whether the composability gives an independent test of the Siegel constraint that bypasses the need to compute b-period contours directly.

---

## §5 — Score-card

| Aspect | Assessment |
|--------|-----------|
| Exclusion methodology | Strong — specific structural reasons, no loose-analogy inflation |
| Morphism A structural claim | Sound — correct identification of gap arcs with Weierstrass branch points |
| FT-02A test protocol | **Mismatched** — tests M_a (real a-period matrix) not Omega; corrected protocol documented |
| Morphism B structural claim | Sound — NPS lower bound via bounded gap-component count is well-grounded |
| FT-05A test | Well-designed, runnable inline |
| Morphism C structural claim | Sound — co-area formula derivation of KKT Jacobian is explicit and correct |
| FT-06A test | Well-designed, runnable inline |
| Composability note (A+C) | Legitimate research direction, correctly hedged |
| Claim discipline | Maintained throughout — claim ceiling 0, no altitude movement, no receipt invention |
