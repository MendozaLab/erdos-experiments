# Closure Plan — `lindstrom_sieve` axiom (Erdős #755 sidecar)

**Date:** 2026-05-02
**Author:** scoping pass (no `.lean` modifications)
**Toolchain:** `leanprover/lean4:v4.27.0`, Mathlib pin `v4.27.0`
**Source axiom:** `lean/Erdos755_Lindstrom.lean:60`
**Companion axiom in #30 main:** `bfr_core_bound` (`lean/Erdos30_BFR.lean:345`) — same proof technology
**Closure tier:** **T3 — research-tier, Mathlib-blocked**
**Verdict:** **Multi-month Mathlib4 contribution opportunity** (preferred), *not* an indefinite axiom

---

## 1. Axiom statement (verbatim)

```lean
/-- **AXIOM — Lindström 1969 sieve bound** ... -/
axiom lindstrom_sieve (A : Finset ℕ) (N g : ℕ)
    (hS : IsB2GSet A g)
    (hA : A ⊆ Finset.range (N + 1)) :
    ∃ c : ℕ, A.card ^ 2 ≤ g * (2 * N + 1) + c * (g * N + 1)
```

`IsB2GSet A g` (from `Erdos755_BhG.lean`): every sum `s ∈ ℕ` has at most `2g` ordered representations `(a,b) ∈ A × A` with `a + b = s`.

**Quantitative statement target.** The literature form is `|A|² ≤ g(2N+1) + O((gN)^{3/4})`. The Lean statement uses the weaker existential `c · (gN+1)` placeholder so the formalization is decoupled from the exact `3/4` exponent — closing the axiom requires either (a) producing such a `c` explicitly (probably linear in `g`) or (b) sharpening the statement to the literature `(gN)^{3/4}` form once the analytic floor exists.

---

## 2. Proof outline (Lindström 1969, J. Combin. Theory 6, 211–212)

The classical character-sum argument has four legs.

**(L1) Set up the trigonometric polynomial.** Let `f : ℝ/ℤ → ℂ`, `f(θ) = Σ_{a ∈ A} e(aθ)`, with `e(x) = exp(2πix)`. Then `|f|² = Σ_{a,b ∈ A} e((a-b)θ)`.

**(L2) Parseval for the L²-norm (gives `|A|`).** `∫₀¹ |f(θ)|² dθ = |A|` since the off-diagonal characters integrate to zero.

**(L3) Quadruple count for the L⁴-norm.** Writing `|f|⁴ = |f²|²` and expanding,
  `∫₀¹ |f(θ)|⁴ dθ = #{(a,b,c,d) ∈ A⁴ : a + b = c + d} =: E(A)` — the **additive energy** of `A`. By the `B_2[g]` hypothesis, every sum-fiber has size `≤ 2g`, so `E(A) = Σ_s r(s)² ≤ 2g · Σ_s r(s) = 2g · |A|²` where `r(s) = #{(a,b) : a+b=s}`. *(This is the elementary bound. Lindström sharpens it.)*

**(L4) Sieve refinement via Cauchy–Schwarz on the representation function.** The improvement to `g(2N+1) + O((gN)^{3/4})` (factor-of-2 in the leading term) comes from:
- Restricting integration over a *minor arc*: split `[0,1] = M ∪ m` (major + minor arcs around rationals with small denominator).
- On `M`, `f` concentrates and contributes the `g(2N+1)` term (Parseval + Bessel).
- On `m`, Cauchy–Schwarz between the L²-norm and the L⁴-norm bounds the contribution by `(2g · |A|²)^{1/2} · meas(m)^{1/2} · (something with N^{1/4})`. Tracking through the optimization gives the `(gN)^{3/4}` correction.

**(L5) Conversion to the existential `c` form.** Given the analytic bound `|A|² ≤ g(2N+1) + C·(gN)^{3/4}`, choose `c = ⌈C⌉` and bound `(gN)^{3/4} ≤ gN+1` for `gN ≥ C^4` (with bookkeeping for small `gN`). This step is purely elementary and would be ~20 lines of `Nat.sqrt` / `omega` once the analytic bound exists.

**Reference proofs in the literature.** Cilleruelo, Ruzsa, Trujillo, "Upper and lower bounds for finite B_h[g] sequences" (J. Number Theory 97, 2002, 26–34) reprises the argument with the constants explicit and is the most Lean-friendly version. Balogh–Füredi–Roy 2023 (cited for `bfr_core_bound`) sharpens further but reuses the same character-sum core.

---

## 3. Mathlib v4.27 inventory — what exists, what is missing

Surveyed `.lake/packages/mathlib/Mathlib/Analysis/Fourier/`, `Mathlib/Combinatorics/Additive/`, `Mathlib/MeasureTheory/Integral/`, `Mathlib/NumberTheory/`.

### Available (reusable as-is)

| Mathlib component | File | Use in Lindström |
|-------------------|------|------------------|
| `AddChar`, `AddChar.linearIndependent`, character orthogonality | `Mathlib/Analysis/Fourier/FiniteAbelian/Orthogonality.lean` | Discrete Plancherel scaffold |
| `AddChar.zmod`, Pontryagin duality on `ZMod n` | `Mathlib/Analysis/Fourier/FiniteAbelian/PontryaginDuality.lean` | Identifies `ZMod`-characters with the dual group |
| `ZMod.dft` (linear equiv), `dft_dft`, `dft_apply_zero`, multiplicative properties | `Mathlib/Analysis/Fourier/ZMod.lean` (225 lines, 25 lemmas) | The discrete Fourier transform itself |
| `Finset.addEnergy s t` = #{(a₁,a₂,b₁,b₂) ∈ s²×t² : a₁+b₁ = a₂+b₂}; `addEnergy_eq_sum_sq` | `Mathlib/Combinatorics/Additive/Energy.lean` | Exactly the L⁴-norm count `∫ \|f\|⁴ = E(A)`, ready to use |
| `Finset.addConvolution`, basic convolution algebra | `Mathlib/Combinatorics/Additive/Convolution.lean` | Representation function `r_A(s) = (1_A * 1_A)(s)` |
| Continuous Plancherel on `AddCircle T`, `hasSum_sq_fourierCoeff` | `Mathlib/Analysis/Fourier/AddCircle.lean:414` | Available but in *continuous* form — bridge to discrete needed |
| `sum_mul_sq_le_sq_mul_sq` (Cauchy–Schwarz over `Finset` in `ℕ`) | used in `Energy.lean:147` | Discrete C-S, ready |

### Missing (the **upstream Mathlib gap**)

| Required result | Status in Mathlib v4.27 | Estimated formalization effort |
|-----------------|-------------------------|-------------------------------|
| **Discrete Parseval** for `ZMod.dft`: `Σ_k \|𝓕Φ k\|² = N · Σ_j \|Φ j\|²` (or the L²-isometry version with normalization) | **Absent.** The `ZMod.dft` file proves inversion (`dft_dft`) but not the L²-norm identity. | ~1 week. One lemma; uses `dft_dft` + character orthogonality from `FiniteAbelian/Orthogonality.lean`. |
| **Discrete L⁴ identity** linking `Σ_k \|𝓕Φ k\|⁴` to `Σ_s \|(Φ ⋆ Φ)(s)\|²` (additive-energy form for general `Φ : ZMod N → ℂ`) | **Absent.** Specialized to indicator functions this is `Σ \|𝓕(1_A)\|⁴ = N² · E(A)`, but the convolution-Fourier link isn't proved at the `ZMod.dft` level. | ~2 weeks. Combines `dft_const_mul`, `addConvolution`, and a discrete `Parseval`. |
| **No `IsSidon` / `B_h[g]` predicate in Mathlib at all** (`grep -r "Sidon\|IsB2"` returns zero hits) | The project's own `IsB2GSet` is the only formal definition. | Could be upstreamed as part of the contribution; ~1 week. |
| **Major-arc / minor-arc decomposition** on `ℝ/ℤ` (the Hardy–Littlewood circle method scaffold) | **Absent at any usable level.** Mathlib has `AddCircle` and `fourierCoeff`, but no Farey-arc partition, no Diophantine-approximation lemmas in the form needed (Dirichlet's approximation theorem exists in some form, but not packaged for arc partitions). | ~6–8 weeks. This is the single biggest gap; it is general-purpose analytic NT infrastructure used for Vinogradov, Waring, and many Erdős-style problems. |
| **Weyl-type bound / Selberg sieve framework** | **Absent** (Mathlib has Selberg sieve density bounds at the level of definitions only — same gap flagged in `H² CLAUDE.md` example for `selberg_sieve_density`). | ~4 weeks if narrowly scoped to additive-energy applications. |
| **Bessel/Cauchy–Schwarz for L² over `AddCircle T`** specialized to indicator-of-Finset functions | Available in continuous form; needs translation lemmas to discrete `ZMod.dft`. | ~1 week. |

### Half-available (close but not directly usable)

| Component | What's there | What's missing |
|-----------|--------------|----------------|
| `Mathlib/Analysis/Fourier/AddCircleMulti.lean` Parseval (line 226, 235) | Inner-product + norm Parseval for L² on the *multidimensional* torus | Single-`ZMod` discrete version not derived from this in Mathlib; a thin wrapper would do it |
| `Mathlib/Analysis/Fourier/LpSpace.lean:86` Plancherel for L² | Continuous Plancherel on `ℝᵈ` | Discrete-finite-group version not obtained from this |

---

## 4. Sub-lemma decomposition (for any closure attempt — internal or upstream)

Closing `lindstrom_sieve` against the existing `Erdos755_BhG.lean` foundation needs roughly the following Lean lemmas. Each is named with a working namespace; lines are estimates assuming Mathlib gaps are filled.

| # | Lemma | Lines | Depends on |
|---|-------|-------|------------|
| L1 | `Erdos.B2G.charSum (A : Finset ℕ) (N : ℕ) : ZMod (2*N+1) → ℂ` (the `f(θ)` analog as a discrete trig polynomial) | ~30 | `ZMod.dft`, `AddChar.zmod` |
| L2 | `Erdos.B2G.charSum_norm_sq_eq_card`: `Σ_θ \|charSum A N θ\|² = (2N+1) · A.card` | ~40 | **Discrete Parseval (Mathlib gap)** |
| L3 | `Erdos.B2G.charSum_norm_4th_eq_addEnergy`: `Σ_θ \|charSum A N θ\|⁴ = (2N+1) · addEnergy A` | ~80 | **Discrete L⁴ identity (Mathlib gap)**, `Finset.addEnergy` |
| L4 | `Erdos.B2G.b2g_addEnergy_le`: `addEnergy A ≤ 2g · A.card²` (already implicit in `b2g_sum_count`; just needs to be stated against `addEnergy`) | ~50 | `Finset.addEnergy_eq_sum_sq`, `IsB2GSet` |
| L5 | `Erdos.B2G.major_arc_contribution`: bound on the Σ over `θ` near 0 (and rationals with small denominator) | ~150 | **Major-arc framework (Mathlib gap)** |
| L6 | `Erdos.B2G.minor_arc_contribution`: Cauchy–Schwarz bound on the rest | ~100 | L4, `sum_mul_sq_le_sq_mul_sq` |
| L7 | `Erdos.B2G.lindstrom_analytic`: combine L5 + L6 to get `\|A\|² ≤ g(2N+1) + C·(gN)^{3/4}` for explicit `C` | ~80 | L5, L6 |
| L8 | `Erdos.B2G.cubic_root_bound_to_existential`: convert `(gN)^{3/4}` → `c · (gN+1)` form | ~30 | `Nat.sqrt`, `omega` |
| L9 | `lindstrom_sieve` (replacing axiom): cite L7 and L8 | ~20 | L7, L8 |

**Total internal effort (assuming Mathlib gaps closed):** ~580 lines of Lean across 9 lemmas.

**Total effort if Mathlib gaps must also be closed by us:** add 8–12 weeks of upstream PR work — primarily the major-arc partition, discrete Parseval/Plancherel for `ZMod.dft`, and the discrete L⁴ ↔ additive-energy bridge.

---

## 5. Effort estimate (agent-cycles, parallelizable)

Three tracks, parallelizable up to the join at L9:

**Track A — Mathlib upstream (critical path):** discrete Parseval for `ZMod.dft` + discrete L⁴-additive-energy identity. ~2 PRs, ~4–6 weeks each through Mathlib review. Independent of the other tracks; usable elsewhere (Cohn–Elkies sphere packing, Tao additive-combinatorics formalizations, Hardy–Littlewood applications).

**Track B — Mathlib upstream (critical path):** major-arc / minor-arc partition + supporting Diophantine-approximation lemmas. ~6–8 weeks. This is the single largest piece of upstream infrastructure and the gating factor.

**Track C — Internal sidecar (joins after A + B):** assemble L1–L9. Once Tracks A and B land, ~2 agent-cycles (2–3 weeks) of focused Lean work to write and tactic-debug the nine lemmas above.

**Total elapsed time, optimistic:** ~10 weeks if Tracks A and B both run in parallel and Mathlib reviewers cooperate.
**Total elapsed time, realistic:** 4–6 months — Mathlib analytic NT contributions historically take 2–3 review cycles.

**Critical-path observation.** This is the *same* Mathlib gap that blocks `bfr_core_bound` in the #30 main package. Any investment in Tracks A and B simultaneously closes both `lindstrom_sieve` (sidecar) and `bfr_core_bound` (main). That double-payoff makes the upstream contribution route the obviously correct economic call.

---

## 6. Verdict and recommendation

**Verdict: T3, Mathlib-blocked, but a multi-month Mathlib4 PR opportunity — *not* an indefinite axiom.**

The axiom should remain in `Erdos755_Lindstrom.lean` until upstream Mathlib lands either (a) discrete Parseval + L⁴ identity for `ZMod.dft`, or (b) the major-arc / minor-arc analytic-NT scaffold. Both are likely multi-month efforts. Until then, the existing `lindstrom_sieve` axiom is the honest accounting.

**Recommendation: pursue this as a Mathlib4 contribution rather than internal-only formalization.**

Three reasons:

1. **Double-payoff.** The same upstream work closes `bfr_core_bound` in the #30 main package. Erdős #30 publication scope already declares both as T3; this is the natural unification.

2. **Community visibility.** Discrete Parseval + additive-energy ↔ L⁴ identity are *general* tools — they unlock formalization of Plünnecke–Ruzsa quantitative variants, Heath-Brown's sieve, Bourgain restriction estimates, and any future Erdős-corpus problem where character sums show up. This compounds the H² portfolio's Lean credibility per Math/CLAUDE.md "Articulate" priority. The internal-only formalization buries the work.

3. **Match to current Mathlib momentum.** Yaël Dillies (author of `Combinatorics/Additive/Energy.lean` and `Convolution.lean`) is actively expanding additive-combinatorics infrastructure in Mathlib. A discrete-Fourier–energy bridge PR would land in a maintained area with a willing reviewer.

**Action items (priority order, no calendar dates):**
1. Open a Mathlib4 issue scoping discrete Parseval + L⁴-additive-energy identity for `ZMod.dft`. Cross-reference Yaël Dillies and the `Energy.lean` / `Convolution.lean` authors as natural reviewers.
2. Draft the simpler PR first: `dft_norm_sq_eq` (discrete Parseval). Standalone, ~1 file, ~150 lines.
3. Follow with `dft_norm_4th_eq_addEnergy`. Builds on (2). ~200 lines.
4. **Decision point:** with (2) and (3) merged, evaluate whether the major-arc machinery (Track B) is worth pursuing as a third Mathlib PR, or whether a *simpler* Lindström variant — using only Parseval + Cauchy–Schwarz, no arc decomposition, gives the weaker bound `|A|² ≤ g(2N+1) + O((gN)^{1/2})` — would already suffice for the existential-`c` form in `lindstrom_sieve`. This weaker variant might close the axiom without Track B at all.
5. Only after (2)–(4): write the 9-lemma internal closure (Track C).

**Until then:** the axiom stays. Cite this `CLOSURE_PLAN_LINDSTROM_SIEVE.md` from `AXIOM_INVENTORY.md` row 6.

---

## 7. Honest scope note

This document is a *plan*, not a closure. Per H² Formalization Integrity Protocol §7, no `lindstrom_sieve` claim has changed: the file still axiomatizes the result, the axiom still has its citation, and the Erdős #755 sidecar still has honest scope ("proof architecture with explicit axioms"). The upgrade in this session is purely the existence of a documented, prioritized roadmap.
