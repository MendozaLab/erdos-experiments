# Erdős #30 — Axiom Inventory

**Last verified:** 2026-05-02 (post-parallel-subagent attack)
**Lean toolchain:** `leanprover/lean4:v4.27.0`
**Mathlib pin:** `v4.27.0` (per `lakefile.lean require`)
**Build status:** PASS for all 5 main targets (`Erdos30_Sidon_Defs`, `Erdos30_Complete`, `Erdos30_Lindstrom`, `Erdos30_BFR`, `Erdos30_Singer`) — 7894 jobs
**Sorry count (main package):** 0
**Scratch sorries (excluded from package):** 2 in `scratch/Erdos30_SpectralSidon.lean` (research stubs, not part of build)
**Singer primes discharged via `native_decide`:** q ∈ {2, 3, 5, 7, 11, 13} (was {2, 3, 5} before this session)

This file is the **single source of truth** for declared axioms in the Erdős #30 Lean package. Per the H² Formalization Integrity Protocol §7, every axiom has an explicit citation and a closure plan. No hidden tactics smuggling deep results.

---

## Axiom inventory (3 in main package + 2 in #755 sidecar; 2 closed this session)

### Active (main package)

| # | Name | File:Line | Statement (informal) | Citation | Closure tier | Mathlib gap |
|---|------|-----------|----------------------|----------|--------------|-------------|
| 1 | `lindstrom_bound` | `lean/Erdos30_Lindstrom.lean:954` | `\|A\| ≤ √N + ⁴√N + 1` for Sidon A ⊆ [0,N], N>0. | Lindström 1969, J. Combin. Theory 7(1) | **DEFERRED — statement gap** | Closure attempt 2026-05-02 found a **statement-level obstruction**: `lindstrom_quadratic` (proved unconditionally) only delivers `+ 2` slack, not `+ 1`. Worked counterexample at N=15 (lindstrom_quadratic gives `k ≤ 6`, but the axiom claims `k ≤ 5` — the empirical truth `k = 5` is a strictly stronger combinatorial fact). **Decision needed:** weaken statement to `+ 2` (cheapest, ~80 lines via Real.sqrt) or commit to a sharper Lindström 1969 full-form argument. See `CLOSURE_ATTEMPT_LINDSTROM.md`. |
| 2 | `bfr_core_bound` | `lean/Erdos30_BFR.lean:345` | Strengthened upper bound: `1000(k−1) ≤ 1000⌊√N⌋ + 998⌊N^{1/4}⌋` for Sidon A ⊆ [0,N]. | Balogh–Füredi–Roy 2023, "On Sidon sets and Lindström's lower bound" | **T2 (downgraded from T3)** | **Closeable in ~24 agent-cycles** using `Finset.addEnergy`, `addEnergy_eq_sum_sq`, `card_sq_le_card_mul_addEnergy` ALREADY in Mathlib v4.27 — no upstream contribution required. Local `bfr_cauchy_schwarz`, `shifted_inter_card_le_one`, `card_distinctSums_sidon`, `lindstrom_quadratic` already proved. 15 sub-lemmas in 5 blocks (A: rep-function setup; B: §2 discrepancy; C: §3 set-systems; D: §4 Cauchy-Schwarz; E: numeric sqrt rounding at N≥10¹²). Critical-path DAG and risk register in `CLOSURE_PLAN_BFR.md`. |
| 3 | `singer_sidon_exists` | `lean/Erdos30_Singer.lean:256` | For every prime q, ∃ A : Finset ℕ Sidon with `\|A\|=q+1` and `∀a∈A, a ≤ q²+q`. | Singer 1938, Trans. AMS 43(3), 377–385 | T3 (novel formalization, **publishable as Mathlib4 PR**) | GaloisField + Singer cycle on PG(2,q) + perfect difference set extraction. **Nowhere formalized in any prover** (Perplexity 2026-03-29). 9–14 agent-cycles total per `CLOSURE_PLAN_SINGER.md`; recommended path is upstream Mathlib4 PR (publishable artifact). **Discharged for q ∈ {2, 3, 5, 7, 11, 13}** via `native_decide` (was {2,3,5} pre-session). q ∈ {17, 19, 23} queued. |

### Sidecar package (Erdős #755, B_h[g] generalization — not part of #30 publication scope)

| # | Name | File:Line | Statement (informal) | Citation | Closure tier |
|---|------|-----------|----------------------|----------|--------------|
| 4 | `singer_b2g_exists` | `lean/Erdos755_Singer_BhG.lean:69` | Singer/Bose–Chowla lower bound for B₂[g] sets: ∃ A : Finset ℕ B₂[g] with `\|A\|² ≥ g(2N+1)/4`. | Bose–Chowla 1962 | T3 — **defer indefinitely** (same `GaloisField` chokepoint as Singer; lowest priority; conditional on Singer closing). See `CLOSURE_PLAN_SINGER_B2G.md`. |
| 5 | `lindstrom_sieve` | `lean/Erdos755_Lindstrom.lean:60` | Lindström sieve refinement: `\|A\|² ≤ g(2N+1) + c·(gN)^{3/4}` for B₂[g]. | Lindström 1969, "An inequality for B₂ sequences", J. Combin. Theory 6, 211–212 | T3 — **3-PR Mathlib4 contribution** (4–6 months). Same upstream gap as `bfr_core_bound`'s alternative path: discrete Parseval + L⁴ ↔ `addEnergy` bridge. **Double-pays**: closing this also closes BFR via the discrete-Fourier route. Yaël Dillies is natural reviewer. See `CLOSURE_PLAN_LINDSTROM_SIEVE.md`. |

---

## Closed in this session (2026-05-02)

| Name | Was at | Closure mechanism |
|------|--------|-------------------|
| `sidon_elem_bound` | `lean/Erdos30_Lindstrom.lean:122` (axiom) | Replaced with `theorem sidon_elem_bound ... := sidon_difference_count A M hS hA` after fixing 7 Mathlib v4.27.0 API-drift errors in `Erdos30_Complete.lean` (impossible-case `absurd` proofs, `simp` → `rw` for `mem_filter`/`mem_product`, `change` for beta-reducing `Finset.card_bij` heq, `Finset.pair_comm` explicit args, `ring` for odd-card witness). |
| `order_diff_counting` | `lean/Erdos30_Lindstrom.lean:431` (axiom) | **Replaced with full `theorem order_diff_counting`** at line 832 via parallel subagent. +454/−19 lines diff. 18 new private declarations: 3 noncomputable defs (`oList`, `oGet`, `orderDiffs`), 15 private lemmas covering sorted enumeration + monotonicity + sigma-bijection cardinality + row-decomposition telescoping sum bound. Inlined sorted-enum helpers rather than patching `Erdos30_OrderedElements.lean` for v4.27. Tactics used: 22 `omega` + 6 `linarith` (all on legitimate Presburger / linear ℤ goals) + `ring` + `Finset` API + `zify`. Build PASS. |

---

## Tier definitions

- **T1 — closeable now:** Mathlib has the needed lemmas; only assembly required.
- **T2 — closeable next round:** 50–150 lines of Lean, no new Mathlib infrastructure needed. Ordered-element enumeration + telescoping arguments fall here.
- **T3 — research-tier:** Requires Mathlib infrastructure that doesn't exist yet (analytic NT sieves, GaloisField primitive root + trace, Singer cycle construction). Multi-week formalization effort. Acceptable to remain axiomatized indefinitely if upstream Mathlib never lands the supporting machinery.

## Closure roadmap (priority order, post-2026-05-02 parallel attack)

1. **`bfr_core_bound`** (downgraded T3 → T2) — newly highest-ROI target. ~24 agent-cycles, no Mathlib contribution required. `Finset.addEnergy` machinery is already there. See `CLOSURE_PLAN_BFR.md` for the 15-lemma DAG and 4-gotcha risk register.
2. **`lindstrom_bound`** — DECISION FORK before any further work. Either (a) weaken statement to `+ 2` slack and close via Real.sqrt chain (~80 lines), or (b) commit to a sharper Lindström 1969 full-form combinatorial argument that recovers the `+ 1` form. See `CLOSURE_ATTEMPT_LINDSTROM.md` for the algebraic obstruction.
3. **`singer_sidon_exists`** (T3) — pursue as Mathlib4 PR (Singer's theorem unformalized in any prover; the upstream PR is itself a publishable artifact). Optionally ladder more native_decide primes (q ∈ {17, 19, 23} already queued) for "Extended Singer construction verified for q ≤ 23" Zenodo companion DOI.
4. **`lindstrom_sieve`** (sidecar, T3) — 3-PR Mathlib4 contribution shared with `bfr_core_bound` via the discrete-Fourier alternative route. Yaël Dillies is the natural reviewer; coordinate sequencing.
5. **`singer_b2g_exists`** (sidecar, T3) — defer indefinitely; conditional on `singer_sidon_exists` Mathlib PR landing first.

---

## Honest scope statement (per H² Formalization Integrity Protocol §7)

The Erdős #30 Lean package is a **proof architecture** in the sense that:
- All theorems closed by tactics or term-level proofs are **fully verified** (0 sorry, lake build PASS, axiom-print shows only Mathlib `propext`/`Classical.choice`/`Quot.sound` plus the 6 declared axioms above).
- The 6 declared axioms are **explicit** with full bibliographic citations — not hidden behind `nlinarith`/`omega`/`norm_num` smuggling deep results.
- Closure of T2 axioms (`order_diff_counting`, `lindstrom_bound`) is roadmapped as next-round work.
- T3 axioms (`bfr_core_bound`, `singer_sidon_exists`, `singer_b2g_exists`, `lindstrom_sieve`) are stable axiomatization of published 1938/1969/2023 results; their closure is not blocking publication of the formalization itself.

This file replaces the scattered axiom-table copies that previously appeared in `README.md`, paper drafts, and per-session reports. Keep this file canonical and reference it from elsewhere.
