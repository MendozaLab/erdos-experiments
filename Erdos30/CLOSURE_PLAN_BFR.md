# Closure Plan — `bfr_core_bound` (Erdős #30 / BFR 2023)

**Document type:** Scoping plan, not a closure attempt.
**Status:** ROADMAP — axiom remains untouched in `lean/Erdos30_BFR.lean:345`.
**Created:** 2026-05-02
**Audience:** Future-Ken / collaborator picking up the BFR formalization.
**Build state at end of this scoping pass:** unchanged (last verified COMPILED 2026-05-02 per `AUDIT_2026-05-02.md`; `lake build Erdos30_BFR` PASS).
**Lean toolchain:** `leanprover/lean4:v4.27.0`. **Mathlib pin:** `v4.27.0`.

---

## 1. Axiom under scope (verbatim)

From `lean/Erdos30_BFR.lean:345`:

```
axiom bfr_core_bound (A : Finset ℕ) (N : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (hN : N ≥ 10^12) :
    1000 * (A.card - 1) ≤ 1000 * Nat.sqrt N + 998 * Nat.sqrt (Nat.sqrt N)
```

Informal mathematical content (Balogh-Füredi-Roy 2023, Theorem 1.1):

> For a Sidon set A ⊆ {0,...,N} with N sufficiently large,
> |A| ≤ √N + (63/64) · ⁴√N + 1.

The integer form multiplies by 1000 to dodge floats and uses 998/1000 = 0.998
as slack over the paper's 63/64 ≈ 0.984, absorbing `Nat.sqrt` rounding losses
at N ≥ 10¹² where ⁴√N ≥ 1000.

**Reference:** Balogh, J., Füredi, Z., Roy, S. (2023). *On Sidon sets and
Lindström's lower bound.* arXiv:2103.15850. Amer. Math. Monthly 130(5), 437-445.

---

## 2. Existing infrastructure (REUSE — no rebuild needed)

Already proven and importable from `Erdos30_BFR.lean` and `Erdos30_Lindstrom.lean`:

| Lemma | File:Line | Role |
|-------|-----------|------|
| `card_distinctSums_sidon` | `Erdos30_BFR.lean:69` | k(k+1)/2 distinct ordered sums (BFR §2 setup) |
| `distinctSums_subset_range` | `Erdos30_BFR.lean:175` | Sums lie in {0,...,2N} |
| `erdos_turan_counting_bound` | `Erdos30_BFR.lean:190` | k(k+1)/2 ≤ 2N+1 (the trivial bound) |
| `card_shifted` | `Erdos30_BFR.lean:208` | Translates preserve cardinality |
| `shifted_subset_range` | `Erdos30_BFR.lean:217` | A+i ⊆ {0,...,N+m} |
| `shifted_inter_card_le_one` | `Erdos30_BFR.lean:229` | Distinct shifts of a Sidon set meet in ≤1 point — **the Sidon-key intersection bound** |
| `bfr_cauchy_schwarz` | `Erdos30_BFR.lean:273` | Variance decomposition: ∑(d-y)² ≥ |X|(d-d_X)² + ∑y² - |X|d_X² |
| `lindstrom_quadratic` | `Erdos30_Lindstrom.lean:455` | ℓ(2k-ℓ-1)² ≤ 4(ℓ+1)N — the Section 2 Lindström quadratic |
| `sidon_distinct_differences` | `Erdos30_Complete.lean` | Sidon ↔ distinct ordered differences |

**Mathlib v4.27 infrastructure available** (verified by inspection):

| Mathlib lemma | Path | Why it matters |
|----|----|----|
| `Finset.addEnergy` | `Mathlib/Combinatorics/Additive/Energy.lean:65` | E[s,t] = |{(a,b,c,d) ∈ s×s×t×t : a+b=c+d}| — the BFR §4 quadruple count |
| `Finset.addEnergy_eq_sum_sq` | same file | E[s] = ∑_{x ∈ s+s} r_s(x)² where r_s(x) = #{(a,b) ∈ s×s : a+b=x} — Cauchy-Schwarz second moment |
| `Finset.le_card_add_mul_addEnergy` | same file | |s|² · |t|² ≤ |s+t| · E[s,t] — Plünnecke-Ruzsa-style energy bound |
| `Finset.card_sq_le_card_mul_addEnergy` | same file | (∑_{c∈u} r(c))² ≤ |u| · E[s,t] — the Cauchy-Schwarz on representation function |
| `Nat.le_sqrt`, `Nat.sqrt_le`, `Nat.sqrt_lt`, `Nat.sqrt_le_sqrt` | `Mathlib/Data/Nat/Sqrt.lean` | All ℕ-sqrt rounding bounds for the final numerical step |
| `Finset.sum_le_sum`, `Finset.sum_lt_sum`, `Finset.mul_sum`, `Finset.sum_const` | `Mathlib.BigOperators` | Standard BigOperators arithmetic (already used heavily in this package) |

**Critical observation:** Mathlib has the Cauchy-Schwarz on the representation
function (`addEnergy_eq_sum_sq` plus `card_sq_le_card_mul_addEnergy`), but it's
stated for additive groups — not for `Finset ℕ` viewed as a subset of an
ambient `ℤ` or `ℕ`-module. A small bridge lemma is needed (Tier-1 work,
~30-60 lines).

**Critical gap:** Mathlib has no `Sidon` definition or `addEnergy A = A.card²`
specialization for Sidon sets. We have `IsSidonSet` locally; need to prove
`E[A] = |A|²` (or `≤ 2|A|² - |A|` if counting trivial quadruples) as a
local bridge lemma.

---

## 3. The fifteen sub-lemmas (BFR §2-4 decomposition)

Numbering follows the proof flow in the BFR paper. Each entry: informal
statement, target Mathlib API, tier, estimated effort. Names follow the
existing `bfr_*` snake_case convention.

### Block A — Representation function setup (5 lemmas)

**A1. `bfr_repFunction` (def)** — *T1, ~10 lines.*
```
def repFunction (A : Finset ℕ) (s : ℕ) : ℕ :=
  ((A ×ˢ A).filter (fun p => p.1 + p.2 = s)).card
```
Number of ordered representations s = a + b with a,b ∈ A.
Target API: `Finset.filter`, `Finset.card_filter`.

**A2. `bfr_repFunction_sidon_le_two` (lemma)** — *T1, ~15 lines.*
For Sidon A and s ∈ ℕ: `repFunction A s ≤ 2`.
Equality is 2 iff s has a non-diagonal representation; 1 iff s = 2a.
Target API: `Finset.card_le_two`, the existing Sidon `IsSidonSet` unfolding.

**A3. `bfr_sum_repFunction_eq_card_sq` (lemma)** — *T1, ~20 lines.*
`∑_{s ∈ Finset.range (2N+1)} repFunction A s = A.card * A.card`
(if A ⊆ {0,...,N}).
This is just |A × A| via fiberwise counting.
Target API: `Finset.sum_card_fiberwise_eq_card_filter`,
`Finset.card_eq_sum_card_fiberwise` (already used in `pigeonhole_residue`).

**A4. `bfr_sum_repFunction_sq_eq_addEnergy` (lemma)** — *T1, ~30 lines.*
`∑_{s ∈ Finset.range (2N+1)} (repFunction A s)² = Finset.addEnergy A A`
(viewing A as a `Finset ℕ`; ℕ has the trivial additive group structure).
Target API: `Finset.addEnergy_eq_sum_sq` from Mathlib.
**Subtlety:** Mathlib's `addEnergy_eq_sum_sq` sums over the sumset s+s; we
sum over the larger `range (2N+1)`. Bridge via the fact that `repFunction A s
= 0` outside `s+s`, plus `Finset.sum_subset` with the zero filler.

**A5. `bfr_addEnergy_sidon` (lemma)** — *T1, ~40 lines.*
For Sidon A: `Finset.addEnergy A A ≤ 2 * A.card * A.card - A.card`
(or equivalent: `addEnergy A ≤ 2 * |A|² - |A|`).
Mechanism: a quadruple (a,b,c,d) with a+b=c+d in a Sidon set forces
{a,b}={c,d}; counting these gives 2|A|²-|A|.
Target API: `IsSidonSet`, `Finset.addEnergy`, `Finset.card_eq_sum_card_fiberwise`.

### Block B — Section 2 (Erdős-Turán with slack) — discrepancy (3 lemmas)

**B1. `bfr_interval_partition` (def + lemma)** — *T1, ~25 lines.*
For a parameter v ≥ 1, partition `Finset.range (2N+1)` into
`⌈(2N+1)/v⌉` consecutive intervals `I_j = {j·v, ..., min((j+1)·v - 1, 2N)}`.
Target API: `Finset.Ioo`, `Finset.range_eq_Ico`,
or just disjoint unions of `Finset.Icc`.

**B2. `bfr_interval_sumcount` (def)** — *T1, ~10 lines.*
`y_j := ∑_{s ∈ I_j} repFunction A s` — number of ordered (a,b) ∈ A×A
landing in interval I_j.
`d := |A|² / v_total` — the average number per interval (over reals or via `Nat.div`).
Target API: `Finset.sum`.

**B3. `bfr_section_2_discrepancy` (theorem, BFR Eq. 2.4)** — *T2, ~80-120 lines.*
For Sidon A ⊆ {0,...,N} and parameter v:
```
∑_j (y_j - d)² ≤ K(A) - some_slack_term(v, N, |A|)
```
where K(A) is a "concentration defect" measuring excess collisions
of A's distinct sums in any interval. This is the **slack form** of the
Section 2 bound.
Target API: combine A4 (sum-of-squares = addEnergy) with B2 and a
manual algebraic identity for `∑(y - d)²`.
**Status note:** This is the most algebra-heavy lemma in Block B but no new
Mathlib infrastructure required. The variance algebra is essentially
`bfr_cauchy_schwarz` applied with `X = univ`.

### Block C — Section 3 (set-systems / Lindström-style) — translates bound (3 lemmas)

**C1. `bfr_translate_disjoint_union_card` (lemma)** — *T1, ~50 lines.*
For Sidon A and 0 ≤ i < j < m: `(shifted A i ∩ shifted A j).card ≤ 1`.
Total: `(⋃_{i<m} shifted A i).card ≥ m·|A| - C(m,2)`.
Target API: `shifted_inter_card_le_one` (already proved!), inclusion-exclusion
via `Finset.card_biUnion_le_card_mul` or pairwise computation.

**C2. `bfr_translate_set_systems` (theorem, BFR Theorem 3.1)** — *T2, ~80 lines.*
For Sidon A ⊆ {0,...,N} and parameter m ≥ 1:
```
A.card² · m ≤ (2*N + m) * (m + A.card - 1) + C(m,2)
```
Or in BFR's form: `k²·m ≤ (2N+m-1)·(m+k-1)`.
Mechanism: count pairs (i, x) with x ∈ shifted A i ∩ shifted A j for some j ≠ i;
LHS counts all such pairs ≥ k²·m, RHS uses the union bound from C1 plus
the pigeonhole on `range (N+m+1)`.
Target API: `shifted_subset_range` (proved), `shifted_inter_card_le_one`
(proved), `Finset.card_biUnion_le`.

**C3. `bfr_section_3_quadratic` (corollary)** — *T1 once C2 closed, ~30 lines.*
Algebraic rearrangement of C2 to: `k² ≤ 2N + 2k·(m-1)/m + (m-1)`.
This isolates the role of m as a tunable parameter.

### Block D — Section 4 (Cauchy-Schwarz combination) — Claims 4.1-4.3 (3 lemmas)

**D1. `bfr_claim_4_1_low_K` (theorem)** — *T2, ~100 lines.*
**Case: K(A) ≤ ε·N^{3/4} (low concentration defect).**
Combining B3 (low K) with the trivial bound k(k+1)/2 ≤ 2N+1 yields
`k ≤ √N + (63/64)·⁴√N + lower-order`.
Target API: `bfr_cauchy_schwarz` (already proved!), `Nat.le_sqrt`, `Nat.sqrt_le_sqrt`.
**Subtlety:** The case-split parameter ε must be chosen explicitly to make
the casework match BFR's 63/64. This is bookkeeping, not new mathematics.

**D2. `bfr_claim_4_2_high_K` (theorem)** — *T2, ~120 lines.*
**Case: K(A) > ε·N^{3/4} (high concentration defect).**
Combining C2 (set-systems with parameter m chosen as ⌊⁴√N⌋) with the
high-K assumption forces k ≤ √N + (63/64)·⁴√N + lower-order via a
contradiction argument.
Target API: `bfr_translate_set_systems` (C2), `Nat.sqrt`, `Nat.lt_succ_sqrt`.

**D3. `bfr_claim_4_3_combine` (theorem)** — *T1 once D1+D2 closed, ~60 lines.*
Casework on K(A): D1 handles low K, D2 handles high K. Conclude:
`1000·(k-1) ≤ 1000·⌊√N⌋ + 998·⌊⁴√N⌋` for N ≥ 10¹².
Target API: pure case-split + `omega` for the integer rounding.

### Block E — Numerical resolution (1 lemma)

**E1. `bfr_sqrt_rounding_at_1e12` (lemma)** — *T1, ~50 lines.*
For N ≥ 10¹²: the conversion from BFR's real-arithmetic 63/64 to the
integer 998/1000 absorbs all `Nat.sqrt` rounding losses (which are ≤ 2 units
on each of `⌊√N⌋` and `⌊⁴√N⌋`).
Target API: `Nat.le_sqrt`, `Nat.lt_succ_sqrt`, `Nat.sqrt_le_self`, `omega`.
**Why it works:** at N = 10¹², ⁴√N ≥ 1000, so the gap (998 vs 1000·63/64
≈ 984) is 16 units of slack on a quantity of size ≥ 1000 — comfortably
absorbing the ≤ 4 units of `Nat.sqrt` rounding.

---

## 4. Mathlib gaps (what is NOT in v4.27)

The following are NOT directly available and would need either a Mathlib
contribution or a local re-proof:

| Need | Status in Mathlib v4.27 | Workaround |
|------|--------------------------|------------|
| `IsSidonSet` definition | **Not in Mathlib** | Local definition exists in `Erdos30_Sidon_Defs.lean`. No upstream contribution needed for this closure. |
| `Finset.addEnergy A A = |A|²` for Sidon | **Not in Mathlib** (no Sidon machinery at all) | Prove locally as A5 above; ~40 lines. |
| Selberg sieve / large sieve | **Not in Mathlib** | **Not needed for BFR.** BFR is elementary; no sieves. |
| Character sums / Weyl inequality | **Not in Mathlib** | **Not needed for BFR.** BFR uses no characters. |
| Real-valued `Nat.sqrt` bounds tying ℕ-sqrt to ℝ-sqrt | Partial: `Nat.cast_sqrt_le_sqrt`, `Real.nat_sqrt_le_real_sqrt` exist but with different name conventions in v4.27 | E1 stays in pure ℕ; no Real bridge needed. |
| Multi-dimensional Cauchy-Schwarz (Plancherel) | `Finset.inner_mul_le_norm_mul_norm` exists | We don't need it; BFR's Cauchy-Schwarz is the variance form already proved as `bfr_cauchy_schwarz`. |

**Headline:** **No Mathlib contribution is required for this closure.**
Every sub-lemma is provable inside the existing Mathlib + local-package
toolset. This is the key reason BFR is Tier-3 long-effort but not Tier-3
"requires upstream PR."

---

## 5. Tier-1 attack path (no new Mathlib needed, closeable next 1-2 sessions)

If a future session wants to *start* closing this axiom without doing the
whole 15-lemma chain, the following are individually closeable today:

| Lemma | Effort (agent-cycles) | Why Tier-1 |
|-------|----------------------|------------|
| **A1** `bfr_repFunction` (def) | 0.1 | Pure definition, ~5 lines. |
| **A2** `bfr_repFunction_sidon_le_two` | 0.5 | Direct from `IsSidonSet` unfolding + case split on |{representations}|. |
| **A3** `bfr_sum_repFunction_eq_card_sq` | 0.5 | Existing `card_eq_sum_card_fiberwise` pattern, mirrors `pigeonhole_residue` proof structure. |
| **A4** `bfr_sum_repFunction_sq_eq_addEnergy` | 1.0 | Bridge to Mathlib's `addEnergy_eq_sum_sq` over the larger range. |
| **A5** `bfr_addEnergy_sidon` | 1.5 | Quadruple count for Sidon A — analogous to existing `card_distinctSums_sidon` |
| **C1** `bfr_translate_disjoint_union_card` | 1.5 | Reuses `shifted_inter_card_le_one` + Bonferroni-style inclusion-exclusion. |
| **E1** `bfr_sqrt_rounding_at_1e12` | 1.0 | Pure `Nat.sqrt` arithmetic + `omega`. Independent of A-D. |

**Total Tier-1 effort:** ~6 agent-cycles. Output: 7 of 15 sub-lemmas closed,
the BFR §2/§3 representation-function and translate-counting machinery in
place. Axiom **still open** but supporting infrastructure is paid down.

## 6. Tier-2 attack path (no new Mathlib needed, but algebra-heavy)

After Tier-1 lands:

| Lemma | Effort | Notes |
|-------|--------|-------|
| **B1** `bfr_interval_partition` | 1.0 | Bookkeeping. |
| **B2** `bfr_interval_sumcount` | 0.5 | Definitional. |
| **B3** `bfr_section_2_discrepancy` | 3.0 | The variance computation; reuses `bfr_cauchy_schwarz`. |
| **C2** `bfr_translate_set_systems` | 3.0 | The set-systems quadratic. |
| **C3** `bfr_section_3_quadratic` | 1.0 | Algebraic rearrangement. |
| **D1** `bfr_claim_4_1_low_K` | 4.0 | The hardest case — full BFR §4.1 |
| **D2** `bfr_claim_4_2_high_K` | 4.0 | The other hard case — full BFR §4.2 |
| **D3** `bfr_claim_4_3_combine` | 1.5 | Final case split. |

**Total Tier-2 effort:** ~18 agent-cycles. Output: axiom **closed** as a
proper theorem, `lake build Erdos30_BFR` continues to PASS, axiom inventory
drops from 4 main-package items to 3.

**Combined Tier-1 + Tier-2 effort:** ~24 agent-cycles total to close the axiom.
This is consistent with AXIOM_INVENTORY's "T3 research-tier, ~weeks" estimate
(at solo-developer pace) but is **achievable** without any upstream Mathlib PR.

## 7. Tier-3 (would require Mathlib contribution)

**None for this axiom.** Every sub-lemma is closeable inside the current
Mathlib v4.27 + local-package toolset. BFR is elementary combinatorics.

---

## 8. Critical-path graph (DAG, parallelizable)

```
                ┌─ A1 ─ A2 ─ A3 ─┐
                │                ├─ A5 ─┐
                └─ A4 ───────────┘      │
                                        ├─ B3 ─┐
                                        │      │
              B1 ─ B2 ───────────────────┘      │
                                                ├─ D1 ─┐
                                                │      │
              C1 ─ C2 ─ C3 ─────────────────────┘      ├─ D3 ─ E1 ─ axiom_closed
                                                │      │
                                                D2 ────┘
```

- **Critical path length:** 6 nodes (A1 → A2 → A3 → A5 → B3 → D1 → D3).
- **Total work:** 15 sub-lemmas across 4 parallelizable tracks.
- **Fan-out at A:** A1/A4 independent → A2/A3 in parallel after A1.
- **Fan-out at root:** A-track / B-track / C-track / E-track all independent.
- **Fan-in at D3:** D1, D2, B3, C3 all merge into the final case-split.

If 4 agents work in parallel: ~6 critical-path cycles + sync overhead ≈
**8-10 wall-clock cycles** to close the axiom.

---

## 9. Honest risk register

1. **The K(A) defect parameter** in B3 / D1 / D2 is the trickiest
   bookkeeping. BFR's paper is terse on the explicit definition; expect to
   spend 1-2 cycles reading the paper carefully before formalizing B3.
2. **`addEnergy` for `Finset ℕ` in an additive monoid context** — Mathlib
   defines additive energy for `[AddCommGroup α]` (or weaker). ℕ is `[AddCommMonoid]`
   not a group, so we may need to embed in ℤ or use a custom definition. Check
   Mathlib's typeclass requirements before assuming A4 is plug-and-play.
   (Risk: A4 jumps from T1 to T2, ~60-100 lines instead of 30.)
3. **The 63/64 vs 998/1000 slack** in E1 — must verify by hand at N = 10¹²
   that the rounding doesn't accidentally tip the inequality. Numerically
   safe at N ≥ 10¹², but the formal proof needs to encode this carefully.
4. **The current `bfr_cauchy_schwarz` lemma** (BFR §4 Lemma 4.1) is stated for
   real-valued y, real d, real d_X. Block D's case analysis requires
   instantiating this with rational d_X = (∑ y) / |X|, which means division.
   Either: (a) keep working in ℝ throughout D and convert only at E1, or
   (b) clear denominators by multiplying through by |X|² before instantiation.
   Recommend (b) for ℕ purity.

---

## 10. What this scoping pass does NOT do

- Does **not** modify `Erdos30_BFR.lean` or any other file in the package.
- Does **not** introduce a `sorry` or `admit` anywhere.
- Does **not** attempt a partial closure (no Tier-1 lemmas added in this session).
- Does **not** require a new Mathlib contribution.
- Does **not** change the build state — `lake build Erdos30_BFR` continues
  to PASS as of last verification (2026-05-02).

The axiom remains untouched. This document is a roadmap.

---

## 11. Next concrete action

When a future session wants to start closure:

1. Pick **A1, A2, A3** as the first batch (lowest risk, ~2 cycles).
2. Run `lake build Erdos30_BFR` between additions to ensure no regression.
3. Update `AXIOM_INVENTORY.md` with progress notes (do not yet downgrade
   the axiom's tier — only when the full chain closes).
4. Commit incrementally; this work spans many sessions and benefits from
   small, reviewable diffs.

**Do not** attempt the Section-4 case analysis (D1/D2) without first having
Block A and one of {B3, C2} fully proved. The case analysis depends on the
representation-function and set-systems infrastructure being airtight.

---

## 12. Cross-references

- `AXIOM_INVENTORY.md` — canonical axiom-status registry (reference, not duplicated here).
- `lean/Erdos30_BFR.lean:345` — the axiom itself.
- `lean/Erdos30_Lindstrom.lean:455` — `lindstrom_quadratic` (the §2 quadratic this builds on).
- BFR 2023 (arXiv:2103.15850) §2-4 — the paper this formalization tracks.
- `Mathlib/Combinatorics/Additive/Energy.lean` — Mathlib's additive-energy library, the key external dependency.
