import Mathlib
import Erdos30_Sidon_Defs
import Erdos30_OrderedElements
import Erdos30_Lindstrom

open Finset Nat

namespace Erdos.Sidon

/-!
  Erdős #30 — Dense finite Sidon ordered-element interface
  =======================================================

  This file is a **research target**, not part of the imported core package.
  It records the current honest formal state of the dense-finite rigidity
  program after checking the literature scale.

  Local, axiom-free combinatorics previously lived here and has been lifted to
  `Erdos30_OrderedElements.lean`. What remains in this scratch file is exactly
  the content that depends on an external-theorem interface: the
  Balasubramanian–Dutta axiom and the theorems derived from it.

  Important:
  - The external ordered-element theorem is not proved locally.
  - The earlier sharp prefix-discrepancy axiom was removed because it was
    stronger than the published Balasubramanian–Dutta scale.
  - The file is kept in `scratch/` until a real local proof route is identified
    for the external interface.

  **Style note.** The axiom and its derived consequences are stated in ℝ-valued
  form (`Real.sqrt`, `Real.rpow`) to match the published paper. The local
  combinatorics imported from `Erdos30_OrderedElements.lean` are ℕ-valued; casts
  are performed at the point of use.
-/

/-!
### Perplexity Pre-Submission Gate — 2026-04-20

**Query:** Verify Balasubramanian–Dutta, *The m-th Element of a Sidon Set*
(`arXiv:2409.01986`, J. Number Theory Vol. 279, DOI
`10.1016/j.jnt.2025.07.007`) Theorem 3 statement, error exponents,
deficiency dependence, and hypotheses.

**Result: PASS — external theorem interface matches Theorem 3 up to O-constant
absorption.**

- Published form: `a_m = m · n^{1/2} + O(n^{7/8}) + O(L^{1/2} · n^{3/4})`,
  where `L = max(0, n^{1/2} - |A|)` and `m` ranges over `{1, ..., |A|}`.
- No macroscopic-index restriction; bound is uniform in `m`.
- Hypotheses: Sidon-ness and the deficiency relation.
- No known counterexamples, gaps, or sharper follow-ups found by the gate as of
  2026-04-20.

**Absorbed discrepancies:**
- Paper: `A ⊆ {1, ..., n}`; this file: `A ⊆ Finset.range (n+1)` = `{0, ..., n}`.
  Shift by `1` in the ambient range is absorbed into the O-constant.
- Paper: `L = max(0, √n - |A|)` with real square root; this file uses the
  integer-floor deficiency relation `A.card + L = Nat.sqrt n`.
  The floor difference is `< 1`, again absorbed into the O-constant.

**Future-proof note.** One proof route in the paper assumes `L ≤ n^(21/80)`,
but the main theorem statement does not. Any attempt to replace the external
interface with a local proof must handle both deficiency regimes.

**Lean-side pitfalls flagged by gate:**
- `n^(7/8)` and `n^(3/4)` are expressed here with `Real.rpow`.
- O-notation is non-constructive; this file packages the hidden constant as an
  existential `C`.
-/

/-- External theorem interface for Balasubramanian–Dutta,
*The m-th Element of a Sidon Set*.

Verified against the open arXiv version `arXiv:2409.01986` and the Journal of
Number Theory publication (Vol. 279, February 2026, DOI
`10.1016/j.jnt.2025.07.007`).

This is not a local proof in the present file. It records the literature-scale
ordered-element estimate that is currently justified for dense Sidon sets.
-/
axiom dense_sidon_ordered_element_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (_hDense : DenseSidonAtScale A n L)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)

/-- General external theorem interface matching the literature deficiency parameter.

Unlike `dense_sidon_ordered_element_external`, this version does not force the
set to live below `floor (sqrt n)`. The truncated real deficiency
`max(0, sqrt n - |A|)` lets the same statement see actual maximizers and
super-floor dense sets as well.
-/
axiom sidon_in_range_ordered_element_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (_hA : SidonInRange A n)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ)

/-! ### Extremal-surface vocabulary

These definitions are deliberately lightweight. The April 24 exact packets
suggest that the prefix/mass compatibility signal is attached to the exact
extremal surface, not to arbitrary dense Sidon layers. The predicates below let
future statements say that directly without yet proving an optimizer-selection
theorem.
-/

/-- `A` is an exact cardinality-maximal Sidon set inside `[0,n]`. -/
def IsMaximalSidonInRange (A : Finset ℕ) (n : ℕ) : Prop :=
  SidonInRange A n ∧ ∀ B : Finset ℕ, SidonInRange B n → B.card ≤ A.card

/-- `A` is within `δ` elements of the maximal Sidon cardinality inside `[0,n]`. -/
def NearExtremalSidonInRange (A : Finset ℕ) (n δ : ℕ) : Prop :=
  SidonInRange A n ∧ ∀ B : Finset ℕ, SidonInRange B n → B.card ≤ A.card + δ

theorem IsMaximalSidonInRange.sidonInRange {A : Finset ℕ} {n : ℕ}
    (hA : IsMaximalSidonInRange A n) : SidonInRange A n :=
  hA.1

theorem NearExtremalSidonInRange.sidonInRange {A : Finset ℕ} {n δ : ℕ}
    (hA : NearExtremalSidonInRange A n δ) : SidonInRange A n :=
  hA.1

theorem IsMaximalSidonInRange.nearExtremal_zero {A : Finset ℕ} {n : ℕ}
    (hA : IsMaximalSidonInRange A n) : NearExtremalSidonInRange A n 0 := by
  refine ⟨hA.sidonInRange, ?_⟩
  intro B hB
  simpa using hA.2 B hB

/-- Prefix residual after the maximizer-friendly terminal drift has been subtracted. -/
noncomputable def prefixResidualAfterGeneralDrift
    (A : Finset ℕ) (n t : ℕ) : ℝ :=
  max 0
    (|(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)| -
      max (realGapFromSqrt A n) 1 * Real.sqrt (n : ℝ))

/-- Density-adjusted mass deviation from the center `n (|A| + 1) / 2`. -/
noncomputable def densityAdjustedMassDeviation (A : Finset ℕ) (n : ℕ) : ℝ :=
  |(∑ a ∈ A, (a : ℝ)) - ((n : ℝ) * ((A.card : ℝ) + 1) / 2)|

/-- General cutpoint consequence of the literature theorem.

This is the maximizer-friendly analogue of `dense_sidon_prefix_cutpoint_external`.
It uses the real deficiency parameter, so it remains meaningful even when
`A.card > floor (sqrt n)`.
-/
theorem sidon_in_range_prefix_cutpoint_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) -
            ((intervalSlice A 0 (orderedElement A i)).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA i
  have hmain := hC hA i
  have hcard :
      (((i.1 + 1 : ℕ) : ℝ)) =
        ((intervalSlice A 0 (orderedElement A i)).card : ℝ) := by
    norm_num [ordered_prefix_card_target A i]
  simpa [hcard, mul_comm, mul_left_comm, mul_assoc] using hmain

/-- Nearby-prefix theorem for the maximizer-friendly interface. -/
theorem sidon_in_range_nearby_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n)
        (i : Fin A.card) (hi : i.1 + 1 < A.card) {t : ℕ},
        orderedElement A i ≤ t →
        t < orderedElement A ⟨i.1 + 1, hi⟩ →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA i hi t hleft hright
  let j : Fin A.card := ⟨i.1 + 1, hi⟩
  let s : ℝ := Real.sqrt (n : ℝ)
  let E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (realDeficiencyFromSqrt A n) *
        Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hcount : (intervalSlice A 0 t).card = i.1 + 1 :=
    Erdos.Sidon.ordered_prefix_card_between_consecutive A i hi hleft hright
  have hpref : ((intervalSlice A 0 t).card : ℝ) = ((i.1 + 1 : ℕ) : ℝ) := by
    norm_num [hcount]
  have hi_abs : |((orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s)| ≤ E := by
    simpa [s, E] using hC hA i
  have hj_abs : |((orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s)| ≤ E := by
    simpa [j, s, E, Nat.cast_add, add_assoc, add_comm, add_left_comm] using hC hA j
  have hleft_real : (orderedElement A i : ℝ) ≤ t := by exact_mod_cast hleft
  have hright_real : (t : ℝ) < orderedElement A j := by exact_mod_cast hright
  have hlowE : -E ≤ (orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
    exact (abs_le.mp hi_abs).1
  have huppE : (orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s ≤ E := by
    exact (abs_le.mp hj_abs).2
  have hlow :
      -(s + E) ≤ (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
    have hmono :
        (orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤
          (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
      nlinarith
    have hstep : -(s + E) ≤ -E := by
      nlinarith [Real.sqrt_nonneg (n : ℝ)]
    exact le_trans hstep (le_trans hlowE hmono)
  have hupp :
      (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤ s + E := by
    have hupp_lt :
        (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s < s + E := by
      have hmid :
          (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s <
            (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
        nlinarith
      have htop :
          (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤ s + E := by
        have hcastsucc : ((i.1 + 2 : ℕ) : ℝ) = ((i.1 + 1 : ℕ) : ℝ) + 1 := by
          push_cast
          ring
        have hrewrite :
            (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s =
              ((orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s) + s := by
          rw [hcastsucc]
          ring_nf
        rw [hrewrite]
        nlinarith
      exact lt_of_lt_of_le hmid htop
    exact le_of_lt hupp_lt
  have habs :
      |(t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ s + E := by
    exact abs_le.mpr ⟨hlow, hupp⟩
  simpa [hpref, s, E, add_assoc, add_left_comm, add_comm, j] using habs

/-- Index-free internal-prefix theorem for the maximizer-friendly interface. -/
theorem sidon_in_range_internal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) {t : ℕ},
        0 < (intervalSlice A 0 t).card →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_nearby_prefix_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA t hs0 hslt
  let i : Fin A.card := ⟨(intervalSlice A 0 t).card - 1, by omega⟩
  have hbracket := Erdos.Sidon.ordered_prefix_bracketing_of_internal A hs0 hslt
  dsimp [i] at hbracket
  rcases hbracket with ⟨hleft, hright⟩
  have hi : i.1 + 1 < A.card := by
    dsimp [i]
    omega
  have hindex : (⟨i.1 + 1, hi⟩ : Fin A.card) = ⟨(intervalSlice A 0 t).card, hslt⟩ := by
    ext
    dsimp [i]
    omega
  have hleft' : orderedElement A i ≤ t := by
    simpa [i] using hleft
  have hright' : t < orderedElement A ⟨i.1 + 1, hi⟩ := by
    simpa [hindex] using hright
  simpa [i] using hC hA i hi hleft' hright'

/-- Empty-prefix boundary theorem for the maximizer-friendly interface. -/
theorem sidon_in_range_empty_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) {t : ℕ},
        (intervalSlice A 0 t).card = 0 →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA t hzero hslt
  have hApos : 0 < A.card := by
    simpa [hzero] using hslt
  let i0 : Fin A.card := ⟨0, hApos⟩
  have hfirst_gt : t < orderedElement A i0 := by
    by_contra hnot
    have hle : orderedElement A i0 ≤ t := Nat.le_of_not_gt hnot
    have hmem : orderedElement A i0 ∈ intervalSlice A 0 t := by
      simp [intervalSlice, orderedElement_mem, hle]
    have hpos : 0 < (intervalSlice A 0 t).card := Finset.card_pos.mpr ⟨orderedElement A i0, hmem⟩
    omega
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hfirst_bd : |(orderedElement A i0 : ℝ) - s| ≤ E := by
    simpa [i0, s, E] using hC hA i0
  have hfirst_le : (orderedElement A i0 : ℝ) ≤ s + E := by
    have hupp := (abs_le.mp hfirst_bd).2
    nlinarith
  have ht_le : (t : ℝ) ≤ s + E := by
    have hfirst_gt_real : (t : ℝ) < orderedElement A i0 := by
      exact_mod_cast hfirst_gt
    nlinarith
  have ht_nonneg : 0 ≤ (t : ℝ) := by
    exact_mod_cast Nat.zero_le t
  have habs :
      |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * s| = (t : ℝ) := by
    simp [hzero, s, abs_of_nonneg ht_nonneg]
  have hmain : |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * s| ≤ s + E := by
    rw [habs]
    exact ht_le
  simpa [s, E, add_assoc, add_left_comm, add_comm] using hmain

/-- Nonterminal-prefix theorem for the maximizer-friendly interface. -/
theorem sidon_in_range_nonterminal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) {t : ℕ},
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_internal_prefix_external with ⟨Cint, hCint, hInt⟩
  rcases sidon_in_range_empty_prefix_external with ⟨Cempty, hCempty, hEmpty⟩
  let C : ℝ := Cint + Cempty
  refine ⟨C, add_nonneg hCint hCempty, ?_⟩
  intro A n hA t hslt
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hmono_int :
      Real.sqrt (n : ℝ) + Cint * u + Cint * v ≤
        Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    have hu : 0 ≤ u := by
      dsimp [u]
      exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
    have hv : 0 ≤ v := by
      dsimp [v]
      exact mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _)
    nlinarith
  have hmono_empty :
      Real.sqrt (n : ℝ) + Cempty * u + Cempty * v ≤
        Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    have hu : 0 ≤ u := by
      dsimp [u]
      exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
    have hv : 0 ≤ v := by
      dsimp [v]
      exact mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _)
    nlinarith
  by_cases hs0 : 0 < (intervalSlice A 0 t).card
  · have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cint * u + Cint * v := by
      simpa [u, v, mul_assoc] using hInt hA hs0 hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_int
    simpa [u, v, mul_assoc] using hmain
  · have hzero : (intervalSlice A 0 t).card = 0 := Nat.eq_zero_of_not_pos hs0
    have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cempty * u + Cempty * v := by
      simpa [u, v, mul_assoc] using hEmpty hA hzero hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_empty
    simpa [u, v, mul_assoc] using hmain

/-- Terminal-prefix theorem for the maximizer-friendly interface. -/
theorem sidon_in_range_terminal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) {t : ℕ},
        0 < A.card →
        (intervalSlice A 0 t).card = A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (Erdos.Sidon.realGapFromSqrt A n) * Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA t hApos hfull ht
  let iLast : Fin A.card := ⟨A.card - 1, by omega⟩
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hEnonneg : 0 ≤ E := by
    dsimp [E]
    exact add_nonneg
      (mul_nonneg hCnonneg (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _))
      (mul_nonneg (mul_nonneg hCnonneg (Real.sqrt_nonneg _))
        (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _))
  have hgap_nonneg : 0 ≤ Erdos.Sidon.realGapFromSqrt A n :=
    Erdos.Sidon.realGapFromSqrt_nonneg A n
  have hdrift_nonneg : 0 ≤ (Erdos.Sidon.realGapFromSqrt A n) * s := by
    exact mul_nonneg hgap_nonneg hs_nonneg
  have hlast_raw := hC hA iLast
  have hlast_est : |(orderedElement A iLast : ℝ) - (A.card : ℝ) * s| ≤ E := by
    have hcard_nat : iLast.1 + 1 = A.card := by
      dsimp [iLast]
      omega
    have hcard : (((iLast.1 + 1 : ℕ) : ℝ)) = (A.card : ℝ) := by
      exact_mod_cast hcard_nat
    simpa [s, E, hcard, mul_assoc] using hlast_raw
  have hlast_le_t : orderedElement A iLast ≤ t := by
    by_contra hnot
    have hgt : t < orderedElement A iLast := Nat.lt_of_not_ge hnot
    have hsubset : intervalSlice A 0 t ⊆ A := by
      intro x hx
      simp [intervalSlice] at hx
      exact hx.1
    have hmemLast : orderedElement A iLast ∈ A := orderedElement_mem A iLast
    have hnotmemLast : orderedElement A iLast ∉ intervalSlice A 0 t := by
      simp [intervalSlice, orderedElement_mem, Nat.not_le_of_lt hgt]
    have hssub : intervalSlice A 0 t ⊂ A := by
      exact Finset.ssubset_iff_subset_ne.mpr ⟨hsubset, by
        intro heq
        exact hnotmemLast (heq.symm ▸ hmemLast)⟩
    have hcard_lt := Finset.card_lt_card hssub
    omega
  have ht_le_real : (t : ℝ) ≤ n := by
    exact_mod_cast ht
  have hsq : (n : ℝ) = s ^ 2 := by
    dsimp [s]
    symm
    exact Real.sq_sqrt (by positivity)
  have hfac :
      s - (A.card : ℝ) ≤ Erdos.Sidon.realGapFromSqrt A n := by
    dsimp [Erdos.Sidon.realGapFromSqrt, s]
    nlinarith [neg_le_abs ((A.card : ℝ) - Real.sqrt (n : ℝ))]
  have hupper0 :
      (t : ℝ) - (A.card : ℝ) * s ≤ (Erdos.Sidon.realGapFromSqrt A n) * s := by
    have hstep : (t : ℝ) - (A.card : ℝ) * s ≤ (n : ℝ) - (A.card : ℝ) * s := by
      nlinarith
    have hdom : (n : ℝ) - (A.card : ℝ) * s ≤ (Erdos.Sidon.realGapFromSqrt A n) * s := by
      rw [hsq]
      nlinarith
    exact le_trans hstep hdom
  have hlast_le_real : (orderedElement A iLast : ℝ) ≤ t := by
    exact_mod_cast hlast_le_t
  have hmono :
      (orderedElement A iLast : ℝ) - (A.card : ℝ) * s ≤
        (t : ℝ) - (A.card : ℝ) * s := by
    nlinarith
  have hlowE : -E ≤ (orderedElement A iLast : ℝ) - (A.card : ℝ) * s := by
    exact (abs_le.mp hlast_est).1
  have hlow :
      -((Erdos.Sidon.realGapFromSqrt A n) * s + E) ≤ (t : ℝ) - (A.card : ℝ) * s := by
    nlinarith [hdrift_nonneg, hlowE, hmono]
  have hupp :
      (t : ℝ) - (A.card : ℝ) * s ≤ (Erdos.Sidon.realGapFromSqrt A n) * s + E := by
    nlinarith [hupper0, hEnonneg]
  have habs :
      |(t : ℝ) - (A.card : ℝ) * s| ≤ (Erdos.Sidon.realGapFromSqrt A n) * s + E := by
    exact abs_le.mpr ⟨hlow, hupp⟩
  simpa [hfull, s, E, add_assoc, add_left_comm, add_comm, mul_assoc] using habs

/-- Honest all-prefix theorem for positive-card Sidon sets in range. -/
theorem sidon_in_range_positive_card_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) {t : ℕ},
        0 < A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (max (Erdos.Sidon.realGapFromSqrt A n) 1) * Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_nonterminal_prefix_external with ⟨Cnon, hCnon, hNon⟩
  rcases sidon_in_range_terminal_prefix_external with ⟨Cterm, hCterm, hTerm⟩
  let C : ℝ := Cnon + Cterm
  refine ⟨C, add_nonneg hCnon hCterm, ?_⟩
  intro A n hA t hApos ht
  set g : ℝ := Erdos.Sidon.realGapFromSqrt A n
  set m : ℝ := max g 1
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hg : 0 ≤ g := by
    dsimp [g]
    exact Erdos.Sidon.realGapFromSqrt_nonneg A n
  have hm1 : (1 : ℝ) ≤ m := by
    dsimp [m]
    exact le_max_right _ _
  have hmg : g ≤ m := by
    dsimp [m]
    exact le_max_left _ _
  have hs_nonneg : 0 ≤ Real.sqrt (n : ℝ) := Real.sqrt_nonneg _
  have hu : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hv : 0 ≤ v := by
    dsimp [v]
    exact mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _)
  have hmono_non :
      Real.sqrt (n : ℝ) + Cnon * u + Cnon * v ≤
        m * Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    nlinarith
  have hmono_term :
      g * Real.sqrt (n : ℝ) + Cterm * u + Cterm * v ≤
        m * Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    nlinarith
  have hsubset : intervalSlice A 0 t ⊆ A := by
    intro x hx
    simp [intervalSlice] at hx
    exact hx.1
  have hcard_le : (intervalSlice A 0 t).card ≤ A.card := Finset.card_le_card hsubset
  by_cases hslt : (intervalSlice A 0 t).card < A.card
  · have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cnon * u + Cnon * v := by
      simpa [u, v, mul_assoc] using hNon hA hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ m * Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_non
    simpa [m, u, v, mul_assoc, add_assoc, add_left_comm, add_comm] using hmain
  · have hfull : (intervalSlice A 0 t).card = A.card := by
      omega
    have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ g * Real.sqrt (n : ℝ) + Cterm * u + Cterm * v := by
      simpa [g, u, v, mul_assoc] using hTerm hA hApos hfull ht
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ m * Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_term
    simpa [m, u, v, mul_assoc, add_assoc, add_left_comm, add_comm] using hmain

/-- In the super-floor regime `|A| ≥ floor (sqrt n)`, the symmetric density gap is
already small on the lower side, and Lindström bounds the upper side by the
fourth-root correction. This is the first honest place where the terminal drift
can be replaced by an `n`-only quantity rather than something still depending on
the particular set. -/
theorem realGapFromSqrt_le_fourthRoot_of_superfloor
    {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n) (hn : 0 < n)
    (hfloor : Nat.sqrt n ≤ A.card) :
    Erdos.Sidon.realGapFromSqrt A n ≤ (Nat.sqrt (Nat.sqrt n) : ℝ) + 1 := by
  dsimp [Erdos.Sidon.realGapFromSqrt]
  have hs_lo : (Nat.sqrt n : ℝ) ≤ Real.sqrt (n : ℝ) :=
    Real.nat_sqrt_le_real_sqrt
  have hs_hi : Real.sqrt (n : ℝ) ≤ (Nat.sqrt n : ℝ) + 1 := by
    simpa using Real.real_sqrt_le_nat_sqrt_succ (a := n)
  have hcard_lo : (Nat.sqrt n : ℝ) ≤ (A.card : ℝ) := by
    exact_mod_cast hfloor
  have hL := lindstrom_bound A n hA.1 hA.2 hn
  have hcard_hi :
      (A.card : ℝ) ≤ (Nat.sqrt n : ℝ) + (Nat.sqrt (Nat.sqrt n) : ℝ) + 1 := by
    exact_mod_cast hL
  have hlower :
      -((Nat.sqrt (Nat.sqrt n) : ℝ) + 1) ≤ (A.card : ℝ) - Real.sqrt (n : ℝ) := by
    nlinarith
  have hupper :
      (A.card : ℝ) - Real.sqrt (n : ℝ) ≤ (Nat.sqrt (Nat.sqrt n) : ℝ) + 1 := by
    nlinarith
  exact abs_le.mpr ⟨hlower, hupper⟩

/-- In the super-floor regime `|A| ≥ floor (sqrt n)`, the truncated real
deficiency is at most `1`, since the only possible deficit below the true
square-root scale comes from the floor gap. -/
theorem realDeficiencyFromSqrt_le_one_of_superfloor
    {A : Finset ℕ} {n : ℕ} (_hA : SidonInRange A n)
    (hfloor : Nat.sqrt n ≤ A.card) :
    realDeficiencyFromSqrt A n ≤ 1 := by
  dsimp [realDeficiencyFromSqrt]
  refine max_le ?_ ?_
  · linarith
  · have hs_hi : Real.sqrt (n : ℝ) ≤ (Nat.sqrt n : ℝ) + 1 := by
      simpa using Real.real_sqrt_le_nat_sqrt_succ (a := n)
    have hcard_lo : (Nat.sqrt n : ℝ) ≤ (A.card : ℝ) := by
      exact_mod_cast hfloor
    nlinarith

/-- Prefix theorem with the endpoint drift absorbed into an `n`-only fourth-root
correction on the super-floor corridor.

This does not yet remove the terminal correction entirely, but it no longer
depends on the particular set `A`. That makes it the honest next formal target
after the maximizer calibration packet: shrink the deterministic endpoint term
from a set-dependent gap to a pure `n`-scale correction, and then test whether
that correction can itself be folded into the literature error.
-/
theorem sidon_in_range_superfloor_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        0 < A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((Nat.sqrt (Nat.sqrt n) : ℕ) : ℝ) + 1) * Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_positive_card_prefix_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n t hA hn hfloor hApos ht
  set s : ℝ := Real.sqrt (n : ℝ)
  set g : ℝ := Erdos.Sidon.realGapFromSqrt A n
  set d : ℝ := realDeficiencyFromSqrt A n
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.rpow (n : ℝ) (3 / 4 : ℝ)
  set q : ℝ := (Nat.sqrt (Nat.sqrt n) : ℝ) + 1
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hu_nonneg : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hv_nonneg : 0 ≤ v := by
    dsimp [v]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hq_one : (1 : ℝ) ≤ q := by
    dsimp [q]
    nlinarith
  have hgap : g ≤ q := by
    dsimp [g, q]
    exact realGapFromSqrt_le_fourthRoot_of_superfloor hA hn hfloor
  have hmax : max g 1 ≤ q := max_le hgap hq_one
  have hdef : d ≤ 1 := by
    dsimp [d]
    exact realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt d ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · have hd_nonneg : 0 ≤ d := by
        dsimp [d]
        exact realDeficiencyFromSqrt_nonneg A n
      simpa using hdef
  have hbase := hC hA hApos ht
  have hterm1 : (max g 1) * s ≤ q * s := by
    exact mul_le_mul_of_nonneg_right hmax hs_nonneg
  have hterm2 : C * Real.sqrt d * v ≤ C * v := by
    have hmul : C * Real.sqrt d ≤ C * 1 := by
      gcongr
    simpa [one_mul, mul_assoc] using mul_le_mul_of_nonneg_right hmul hv_nonneg
  have hmain :
      (max g 1) * s + C * u + C * Real.sqrt d * v ≤ q * s + C * u + C * v := by
    nlinarith
  exact le_trans hbase (by
    simpa [s, g, d, u, v, q, mul_assoc, add_assoc, add_left_comm, add_comm] using hmain)

/-- The pure ambient correction from `sidon_in_range_superfloor_prefix_external`
already lives below the `n^(7/8)` scale for every positive `n`. -/
theorem superfloor_ambient_le_two_sevenEighths {n : ℕ} (hn : 0 < n) :
    (((Nat.sqrt (Nat.sqrt n) : ℕ) : ℝ) + 1) * Real.sqrt (n : ℝ) ≤
      2 * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  have hn0 : 0 ≤ (n : ℝ) := by exact_mod_cast Nat.zero_le n
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hs_le :
      Real.sqrt (n : ℝ) ≤ Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hquarter_le :
      (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
    have h1 : (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.sqrt (Nat.sqrt n : ℝ) :=
      Real.nat_sqrt_le_real_sqrt
    have h2 : Real.sqrt (Nat.sqrt n : ℝ) ≤ Real.sqrt (Real.sqrt (n : ℝ)) := by
      apply Real.sqrt_le_sqrt
      exact_mod_cast Real.nat_sqrt_le_real_sqrt
    have h3 : Real.sqrt (Real.sqrt (n : ℝ)) = Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
      rw [Real.sqrt_eq_rpow, Real.sqrt_eq_rpow, ← Real.rpow_mul hn0]
      norm_num
    exact le_trans h1 (h2.trans_eq h3)
  have hquarter_mul :
      (Nat.sqrt (Nat.sqrt n) : ℝ) * Real.sqrt (n : ℝ) ≤
        Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
    calc
      (Nat.sqrt (Nat.sqrt n) : ℝ) * Real.sqrt (n : ℝ)
          ≤ Real.rpow (n : ℝ) (1 / 4 : ℝ) * Real.sqrt (n : ℝ) := by
            gcongr
      _ = Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
            rw [Real.sqrt_eq_rpow]
            have htmp :
                Real.rpow (n : ℝ) ((1 / 2 : ℝ) + (1 / 4 : ℝ)) =
                  Real.rpow (n : ℝ) (1 / 2 : ℝ) * Real.rpow (n : ℝ) (1 / 4 : ℝ) :=
              Real.rpow_add_of_nonneg hn0 (by norm_num) (by norm_num)
            calc
              Real.rpow (n : ℝ) (1 / 4 : ℝ) * Real.rpow (n : ℝ) (1 / 2 : ℝ)
                  = Real.rpow (n : ℝ) (1 / 2 : ℝ) * Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
                      ring
              _ = Real.rpow (n : ℝ) ((1 / 2 : ℝ) + (1 / 4 : ℝ)) := by
                    exact htmp.symm
              _ = Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
                    congr 1
                    norm_num
      _ ≤ Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
            exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hexpand :
      (((Nat.sqrt (Nat.sqrt n) : ℕ) : ℝ) + 1) * Real.sqrt (n : ℝ) =
        (Nat.sqrt (Nat.sqrt n) : ℝ) * Real.sqrt (n : ℝ) + Real.sqrt (n : ℝ) := by
    ring
  rw [hexpand]
  have hsum := add_le_add hquarter_mul hs_le
  simpa [two_mul] using hsum

/-- Coarse super-floor prefix theorem with the remaining deterministic correction
absorbed into a single `n^(7/8)` scale.

This is the first compiled theorem in the current lane where the prefix bound is
stated with no explicit set-dependent drift and no extra fourth-root ambient
term. The price is a coarser constant, but the scale is now unified.
-/
theorem sidon_in_range_superfloor_prefix_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        0 < A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_prefix_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 2 + 2 * C
  refine ⟨C', by nlinarith, ?_⟩
  intro A n t hA hn hfloor hApos ht
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  have hn0 : 0 ≤ (n : ℝ) := by exact_mod_cast Nat.zero_le n
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hthreeFour_le :
      Real.rpow (n : ℝ) (3 / 4 : ℝ) ≤ u := by
    dsimp [u]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hCu :
      C * Real.rpow (n : ℝ) (3 / 4 : ℝ) ≤ C * u := by
    exact mul_le_mul_of_nonneg_left hthreeFour_le hCnonneg
  calc
    |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
        ≤ (((Nat.sqrt (Nat.sqrt n) : ℕ) : ℝ) + 1) * Real.sqrt (n : ℝ) +
            C * u + C * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
              simpa [u] using hC hA hn hfloor hApos ht
    _ ≤ 2 * u + C * u + C * u := by
          nlinarith [superfloor_ambient_le_two_sevenEighths hn, hCu]
    _ = C' * u := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
          rfl

/-- Cleaner super-floor coarse theorem: once `0 < n` and `floor (sqrt n) ≤ |A|`,
the positivity of `A.card` is automatic, so the extra cardinality hypothesis can
be dropped from the user-facing statement. -/
theorem sidon_in_range_superfloor_prefix_coarse_external' :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_prefix_coarse_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n t hA hn hfloor ht
  have hsqrt_pos : 0 < Nat.sqrt n := Nat.sqrt_pos.mpr hn
  have hApos : 0 < A.card := by
    omega
  exact hC hA hn hfloor hApos ht

/-- Maximizer-friendly index-difference theorem.

This is the `SidonInRange` analogue of `dense_sidon_index_difference_external`.
It upgrades the single-index ordered-element theorem to a two-index displacement
statement and is the natural cutpoint-to-interval bridge in the wider regime.
-/
theorem sidon_in_range_index_difference_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ))
            - ((j.1 - i.1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA i j hij
  have hi_bd := hC hA i
  have hj_bd := hC hA j
  set a_i : ℝ := (orderedElement A i : ℝ)
  set a_j : ℝ := (orderedElement A j : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hi_abs : |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ E := hi_bd
  have hj_abs : |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| ≤ E := hj_bd
  have hdiff :
      ((j.1 - i.1 : ℕ) : ℝ) = ((j.1 + 1 : ℕ) : ℝ) - ((i.1 + 1 : ℕ) : ℝ) := by
    rw [Nat.cast_sub hij]
    push_cast
    ring
  have hkey :
      (a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s =
        (a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
          (a_i - ((i.1 + 1 : ℕ) : ℝ) * s) := by
    rw [hdiff]
    ring
  calc
    |(a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s|
        = |(a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
            (a_i - ((i.1 + 1 : ℕ) : ℝ) * s)| := by rw [hkey]
    _ ≤ |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| +
          |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| := abs_sub _ _
    _ ≤ E + E := add_le_add hj_abs hi_abs
    _ = 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
          2 * C * Real.sqrt (realDeficiencyFromSqrt A n) *
            Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
        show (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) +
              (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) = _
        ring

/-- Cutpoint-interval theorem in the maximizer-friendly `SidonInRange` lane. -/
theorem sidon_in_range_cutpoint_interval_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ)) -
            (((intervalSlice A 0 (orderedElement A j)).card -
                (intervalSlice A 0 (orderedElement A i)).card : ℕ) : ℝ) *
              Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases sidon_in_range_index_difference_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA i j hij
  have hmain := hC hA i j hij
  have hcount :
      (((intervalSlice A 0 (orderedElement A j)).card -
          (intervalSlice A 0 (orderedElement A i)).card : ℕ) : ℝ) =
        ((j.1 - i.1 : ℕ) : ℝ) := by
    have hi_card : (intervalSlice A 0 (orderedElement A i)).card = i.1 + 1 :=
      ordered_prefix_card_target A i
    have hj_card : (intervalSlice A 0 (orderedElement A j)).card = j.1 + 1 :=
      ordered_prefix_card_target A j
    have hnat :
        (intervalSlice A 0 (orderedElement A j)).card -
            (intervalSlice A 0 (orderedElement A i)).card = j.1 - i.1 := by
      omega
    exact_mod_cast hnat
  simpa [hcount, mul_comm, mul_left_comm, mul_assoc] using hmain

/-- Super-floor coarse cutpoint-interval theorem.

On the super-floor corridor, the maximizer-friendly cutpoint-interval theorem
collapses to the single `n^(7/8)` scale with no remaining deficiency term.
-/
theorem sidon_in_range_superfloor_cutpoint_interval_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        0 < n →
        Nat.sqrt n ≤ A.card →
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ)) -
            (((intervalSlice A 0 (orderedElement A j)).card -
                (intervalSlice A 0 (orderedElement A i)).card : ℕ) : ℝ) *
              Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_cutpoint_interval_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 4 * C
  refine ⟨C', by positivity, ?_⟩
  intro A n hA i j hij hn hfloor
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hdef : realDeficiencyFromSqrt A n ≤ 1 :=
    realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt (realDeficiencyFromSqrt A n) ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · simpa using hdef
  have hthreeFour_le :
      Real.rpow (n : ℝ) (3 / 4 : ℝ) ≤ u := by
    dsimp [u]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hthreeFour_nonneg : 0 ≤ Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hmix :
      2 * C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
        ≤ 2 * C * u := by
    have hsmall :
        Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ) ≤ u := by
      calc
        Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
            ≤ 1 * u := by
              gcongr
        _ = u := by ring
    have hcoeff_nonneg : 0 ≤ 2 * C := by positivity
    simpa [mul_assoc] using mul_le_mul_of_nonneg_left hsmall hcoeff_nonneg
  calc
    |((orderedElement A j : ℝ) - (orderedElement A i : ℝ)) -
        (((intervalSlice A 0 (orderedElement A j)).card -
            (intervalSlice A 0 (orderedElement A i)).card : ℕ) : ℝ) *
          Real.sqrt (n : ℝ)|
        ≤ 2 * C * u +
            2 * C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
              simpa [u] using hC hA i j hij
    _ ≤ 2 * C * u + 2 * C * u := by
          gcongr
    _ = C' * u := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
          rfl

/-- Coarse internal-prefix theorem on the super-floor corridor.

For a prefix that is genuinely internal (`0 < s(t) < |A|`), the maximizer-friendly
nearby-prefix theorem collapses to a single `n^(7/8)` scale with no endpoint drift
and no extra ambient correction. This is the first index-free nearby-prefix theorem
in the current lane that talks directly about arbitrary internal cut points rather
than only about ordered cutpoints or explicit consecutive brackets.
-/
theorem sidon_in_range_superfloor_internal_prefix_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        0 < (intervalSlice A 0 t).card →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_internal_prefix_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 1 + 2 * C
  refine ⟨C', by nlinarith, ?_⟩
  intro A n t hA hn hfloor hs0 hslt
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by
    exact_mod_cast hn
  have hsqrt_le :
      Real.sqrt (n : ℝ) ≤ u := by
    dsimp [u]
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hdef : realDeficiencyFromSqrt A n ≤ 1 :=
    realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt (realDeficiencyFromSqrt A n) ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · simpa using hdef
  have hthreeFour_le :
      v ≤ u := by
    dsimp [u, v]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hv_nonneg : 0 ≤ v := by
    dsimp [v]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hmix :
      C * Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ C * u := by
    have hsmall :
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ u := by
      calc
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ 1 * u := by
          gcongr
        _ = u := by ring
    simpa [mul_assoc] using mul_le_mul_of_nonneg_left hsmall hCnonneg
  calc
    |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
        ≤ Real.sqrt (n : ℝ) + C * u +
            C * Real.sqrt (realDeficiencyFromSqrt A n) * v := by
              simpa [u, v, mul_assoc] using hC hA hs0 hslt
    _ ≤ u + C * u + C * u := by
          nlinarith [hsqrt_le, hmix]
    _ = C' * u := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
          rfl

/-- Coarse empty-prefix theorem on the super-floor corridor.

This isolates the left boundary case of the prefix package. Once `0 < n` and
`floor (sqrt n) ≤ |A|`, the first ordered element is already controlled at the
literature scale, so a prefix with zero points also sits at the same single
`n^(7/8)` scale.
-/
theorem sidon_in_range_superfloor_empty_prefix_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        (intervalSlice A 0 t).card = 0 →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_empty_prefix_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 1 + 2 * C
  refine ⟨C', by nlinarith, ?_⟩
  intro A n t hA hn hfloor hzero
  have hsqrt_pos : 0 < Nat.sqrt n := Nat.sqrt_pos.mpr hn
  have hApos : 0 < A.card := by
    omega
  have hslt : (intervalSlice A 0 t).card < A.card := by
    rw [hzero]
    omega
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by
    exact_mod_cast hn
  have hsqrt_le :
      Real.sqrt (n : ℝ) ≤ u := by
    dsimp [u]
    rw [Real.sqrt_eq_rpow]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hdef : realDeficiencyFromSqrt A n ≤ 1 :=
    realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt (realDeficiencyFromSqrt A n) ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · simpa using hdef
  have hthreeFour_le :
      v ≤ u := by
    dsimp [u, v]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hv_nonneg : 0 ≤ v := by
    dsimp [v]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hmix :
      C * Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ C * u := by
    have hsmall :
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ u := by
      calc
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ 1 * u := by
          gcongr
        _ = u := by ring
    simpa [mul_assoc] using mul_le_mul_of_nonneg_left hsmall hCnonneg
  calc
    |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
        ≤ Real.sqrt (n : ℝ) + C * u +
            C * Real.sqrt (realDeficiencyFromSqrt A n) * v := by
              simpa [u, v, mul_assoc] using hC hA hzero hslt
    _ ≤ u + C * u + C * u := by
          nlinarith [hsqrt_le, hmix]
    _ = C' * u := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
          rfl

/-- Coarse nonterminal-prefix theorem on the super-floor corridor.

This packages the empty-prefix and internal-prefix branches into one clean
statement with no endpoint drift. It is the natural local bulk theorem: every
prefix that has not yet swallowed all of `A` sits at the single `n^(7/8)` scale.
-/
theorem sidon_in_range_superfloor_nonterminal_prefix_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_internal_prefix_coarse_external with
    ⟨Cint, hCint, hInt⟩
  rcases sidon_in_range_superfloor_empty_prefix_coarse_external with
    ⟨Cempty, hCempty, hEmpty⟩
  let C : ℝ := Cint + Cempty
  refine ⟨C, add_nonneg hCint hCempty, ?_⟩
  intro A n t hA hn hfloor hslt
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  have hu : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  by_cases hs0 : 0 < (intervalSlice A 0 t).card
  · have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Cint * u := by
      simpa [u] using hInt hA hn hfloor hs0 hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * u := by
      dsimp [C]
      nlinarith
    simpa [u] using hmain
  · have hzero : (intervalSlice A 0 t).card = 0 := Nat.eq_zero_of_not_pos hs0
    have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Cempty * u := by
      simpa [u] using hEmpty hA hn hfloor hzero
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * u := by
      dsimp [C]
      nlinarith
    simpa [u] using hmain

/-- Coarse terminal-prefix theorem on the super-floor corridor.

This isolates the right boundary case. Once the density gap and truncated
deficiency are both absorbed by the super-floor estimates, the terminal branch
also collapses to the single `n^(7/8)` scale.
-/
theorem sidon_in_range_superfloor_terminal_prefix_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        (intervalSlice A 0 t).card = A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
  rcases sidon_in_range_terminal_prefix_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 2 + 2 * C
  refine ⟨C', by nlinarith, ?_⟩
  intro A n t hA hn hfloor hfull ht
  have hsqrt_pos : 0 < Nat.sqrt n := Nat.sqrt_pos.mpr hn
  have hApos : 0 < A.card := by
    omega
  set s : ℝ := Real.sqrt (n : ℝ)
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.rpow (n : ℝ) (3 / 4 : ℝ)
  set g : ℝ := Erdos.Sidon.realGapFromSqrt A n
  set q : ℝ := (Nat.sqrt (Nat.sqrt n) : ℝ) + 1
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by
    exact_mod_cast hn
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hgap : g ≤ q := by
    dsimp [g, q]
    exact realGapFromSqrt_le_fourthRoot_of_superfloor hA hn hfloor
  have hambient :
      q * s ≤ 2 * u := by
    simpa [q, s, u] using superfloor_ambient_le_two_sevenEighths hn
  have hgapterm :
      g * s ≤ 2 * u := by
    exact le_trans (mul_le_mul_of_nonneg_right hgap hs_nonneg) hambient
  have hdef : realDeficiencyFromSqrt A n ≤ 1 :=
    realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt (realDeficiencyFromSqrt A n) ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · simpa using hdef
  have hthreeFour_le :
      v ≤ u := by
    dsimp [u, v]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hv_nonneg : 0 ≤ v := by
    dsimp [v]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hmix :
      C * Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ C * u := by
    have hsmall :
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ u := by
      calc
        Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ 1 * u := by
          gcongr
        _ = u := by ring
    simpa [mul_assoc] using mul_le_mul_of_nonneg_left hsmall hCnonneg
  calc
    |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
        ≤ g * s + C * u +
            C * Real.sqrt (realDeficiencyFromSqrt A n) * v := by
              simpa [g, s, u, v, mul_assoc] using hC hA hApos hfull ht
    _ ≤ 2 * u + C * u + C * u := by
          nlinarith [hgapterm, hmix]
    _ = C' * u := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (7 / 8 : ℝ) := by
          rfl

/-- Literature consequence at the actual cutpoints `t = a_i`.

Combining the external Balasubramanian–Dutta theorem with the exact local
identity `|(A ∩ [0, a_i])| = i+1` gives a genuine prefix-count statement at the
ordered-element cutpoints. This is weaker than a uniform prefix discrepancy
theorem, but it is an honest consequence of the verified literature interface.
-/
theorem dense_sidon_prefix_cutpoint_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) -
            ((intervalSlice A 0 (orderedElement A i)).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i
  have hmain := hC hDense i
  have hcard :
      (((i.1 + 1 : ℕ) : ℝ)) =
        ((intervalSlice A 0 (orderedElement A i)).card : ℝ) := by
    norm_num [ordered_prefix_card_target A i]
  simpa [hcard, mul_comm, mul_left_comm, mul_assoc] using hmain

/-- Consecutive-gap corollary of the external ordered-element theorem.

Two applications of Balasubramanian–Dutta at indices `i` and `i+1`, stitched
by the triangle inequality, give that consecutive ordered elements differ by
`√n` up to twice the Balasubramanian–Dutta error. This is the stepping stone
to a nearby-prefix theorem for general `t`: once the gap is controlled, the
prefix count cannot jump more than `√n + error` across any single step, so
moving `t` off a cutpoint perturbs `|t − s(t)·√n|` by at most one gap.

No new axioms used; this is pure triangle inequality on
`dense_sidon_ordered_element_external`.
-/
theorem dense_sidon_consecutive_gap_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i : Fin A.card) (hi : i.1 + 1 < A.card),
        |((orderedElement A ⟨i.1 + 1, hi⟩ : ℝ) - (orderedElement A i : ℝ))
            - Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i hi
  set j : Fin A.card := ⟨i.1 + 1, hi⟩ with hj
  have hi_est := hC hDense i
  have hj_est := hC hDense j
  set a_i : ℝ := (orderedElement A i : ℝ)
  set a_j : ℝ := (orderedElement A j : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hi_bd : |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ E := hi_est
  have hj_bd : |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| ≤ E := hj_est
  have hjval : ((j.1 + 1 : ℕ) : ℝ) = ((i.1 + 1 : ℕ) : ℝ) + 1 := by
    push_cast [hj]
    ring
  have hj_bd' : |a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s| ≤ E := by
    simpa [hjval] using hj_bd
  have hkey :
      (a_j - a_i) - s =
        (a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s) -
          (a_i - ((i.1 + 1 : ℕ) : ℝ) * s) := by ring
  calc
    |(a_j - a_i) - s|
        = |(a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s) -
            (a_i - ((i.1 + 1 : ℕ) : ℝ) * s)| := by rw [hkey]
    _ ≤ |a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s| +
          |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| := abs_sub _ _
    _ ≤ E + E := add_le_add hj_bd' hi_bd
    _ = 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
          2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
        show (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) +
              (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) = _
        ring

/-- Index-difference bound: multi-step generalization of the consecutive gap.

For any two ordered-element indices `i ≤ j`, the displacement `a_j − a_i` is
within `2·E` of `(j − i)·√n`, where `E` is the Balasubramanian–Dutta error.
Two applications of the axiom at `i` and `j` plus triangle inequality.

This subsumes `dense_sidon_consecutive_gap_external` (`j = i + 1` case) and
is the natural tool for bounding `t` that sits between two ordered elements:
once `t` is bracketed by `a_i ≤ t ≤ a_j`, interval control follows from this
lemma plus monotonicity of the prefix count.

The gap between two cutpoints `(j − i) · √n` grows linearly with index
separation, while the error stays at B–D scale — so widely separated cutpoints
give a sharper relative bound than consecutive ones.

No new axioms used.
-/
theorem dense_sidon_index_difference_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ))
            - ((j.1 - i.1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i j hij
  have hi_bd := hC hDense i
  have hj_bd := hC hDense j
  set a_i : ℝ := (orderedElement A i : ℝ)
  set a_j : ℝ := (orderedElement A j : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hi_abs : |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ E := hi_bd
  have hj_abs : |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| ≤ E := hj_bd
  have hdiff :
      ((j.1 - i.1 : ℕ) : ℝ) = ((j.1 + 1 : ℕ) : ℝ) - ((i.1 + 1 : ℕ) : ℝ) := by
    rw [Nat.cast_sub hij]
    push_cast
    ring
  have hkey :
      (a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s =
        (a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
          (a_i - ((i.1 + 1 : ℕ) : ℝ) * s) := by
    rw [hdiff]; ring
  calc
    |(a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s|
        = |(a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
            (a_i - ((i.1 + 1 : ℕ) : ℝ) * s)| := by rw [hkey]
    _ ≤ |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| +
          |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| := abs_sub _ _
    _ ≤ E + E := add_le_add hj_abs hi_abs
    _ = 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
          2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
        show (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) +
              (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) = _
        ring

/-- Cutpoint-interval version of the external index-difference theorem.

For two ordered cutpoints `a_i ≤ a_j`, the displacement `a_j - a_i` is close to
the prefix-count difference times `√n`. This is the exact prefix-count rewrite
of `dense_sidon_index_difference_external`, obtained by substituting the local
identity `|(A ∩ [0, a_k])| = k+1` at both endpoints.

This is still only a theorem at cutpoints, not for arbitrary `t`, but it is the
right bridge toward a future nearby-prefix theorem.
-/
theorem dense_sidon_cutpoint_interval_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ)) -
            (((intervalSlice A 0 (orderedElement A j)).card -
                (intervalSlice A 0 (orderedElement A i)).card : ℕ) : ℝ) *
              Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_index_difference_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i j hij
  have hmain := hC hDense i j hij
  have hcount :
      (intervalSlice A 0 (orderedElement A j)).card -
          (intervalSlice A 0 (orderedElement A i)).card = j.1 - i.1 := by
    rw [ordered_prefix_card_target A j, ordered_prefix_card_target A i]
    omega
  simpa [hcount, mul_comm, mul_left_comm, mul_assoc] using hmain

/-- Nearby-prefix theorem between consecutive ordered cutpoints.

If `t` lies between `a_i` and `a_{i+1}`, then the prefix count is exactly
`i+1`, and `t` is within one `√n` step plus the Balasubramanian–Dutta error of
that prefix count times `√n`.

This is the first honest theorem in the file for arbitrary `t` rather than only
for the cutpoints `t = a_i`.
-/
theorem dense_sidon_nearby_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i : Fin A.card) (hi : i.1 + 1 < A.card) {t : ℕ},
        orderedElement A i ≤ t →
        t < orderedElement A ⟨i.1 + 1, hi⟩ →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i hi t hleft hright
  let j : Fin A.card := ⟨i.1 + 1, hi⟩
  let s : ℝ := Real.sqrt (n : ℝ)
  let E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hcount : (intervalSlice A 0 t).card = i.1 + 1 :=
    Erdos.Sidon.ordered_prefix_card_between_consecutive A i hi hleft hright
  have hpref : ((intervalSlice A 0 t).card : ℝ) = ((i.1 + 1 : ℕ) : ℝ) := by
    norm_num [hcount]
  have hi_abs : |((orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s)| ≤ E := by
    simpa [s, E] using hC hDense i
  have hj_abs : |((orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s)| ≤ E := by
    simpa [j, s, E, Nat.cast_add, add_assoc, add_comm, add_left_comm] using hC hDense j
  have hs_nonneg : 0 ≤ s := by
    exact Real.sqrt_nonneg _
  have hleft_real : (orderedElement A i : ℝ) ≤ t := by exact_mod_cast hleft
  have hright_real : (t : ℝ) < orderedElement A j := by exact_mod_cast hright
  have hlowE : -E ≤ (orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
    exact (abs_le.mp hi_abs).1
  have huppE : (orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s ≤ E := by
    exact (abs_le.mp hj_abs).2
  have hlow :
      - (s + E) ≤ (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
    have hmono :
        (orderedElement A i : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤
          (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
      nlinarith
    have hstep : -(s + E) ≤ -E := by
      nlinarith
    exact le_trans hstep (le_trans hlowE hmono)
  have hupp :
      (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤ s + E := by
    have hupp_lt :
        (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s < s + E := by
      have hmid :
          (t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s <
            (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s := by
        nlinarith
      have htop :
          (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s ≤ s + E := by
        have hcastsucc : ((i.1 + 2 : ℕ) : ℝ) = ((i.1 + 1 : ℕ) : ℝ) + 1 := by
          push_cast
          ring
        have hrewrite :
            (orderedElement A j : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s =
              ((orderedElement A j : ℝ) - ((i.1 + 2 : ℕ) : ℝ) * s) + s := by
          rw [hcastsucc]
          ring_nf
        rw [hrewrite]
        nlinarith
      exact lt_of_lt_of_le hmid htop
    exact le_of_lt hupp_lt
  have habs :
      |(t : ℝ) - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ s + E := by
    exact abs_le.mpr ⟨hlow, hupp⟩
  simpa [hpref, s, E, add_assoc, add_left_comm, add_comm, j] using habs

/-- Index-free nearby-prefix theorem.

If the internal prefix count `s(t) = |A ∩ [0,t]|` satisfies `0 < s(t) < |A|`,
then `t` is within one `√n` step plus Balasubramanian–Dutta error of
`s(t)·√n`. This packages `dense_sidon_nearby_prefix_external` without exposing
the bracketing index.
-/
theorem dense_sidon_internal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L) {t : ℕ},
        0 < (intervalSlice A 0 t).card →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_nearby_prefix_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense t hs0 hslt
  let i : Fin A.card := ⟨(intervalSlice A 0 t).card - 1, by omega⟩
  have hbracket := Erdos.Sidon.ordered_prefix_bracketing_of_internal A hs0 hslt
  dsimp [i] at hbracket
  rcases hbracket with ⟨hleft, hright⟩
  have hi : i.1 + 1 < A.card := by
    dsimp [i]
    omega
  have hs : i.1 + 1 = (intervalSlice A 0 t).card := by
    dsimp [i]
    omega
  have hindex : (⟨i.1 + 1, hi⟩ : Fin A.card) = ⟨(intervalSlice A 0 t).card, hslt⟩ := by
    ext
    exact hs
  have hleft' : orderedElement A i ≤ t := by
    simpa [i] using hleft
  have hright' : t < orderedElement A ⟨i.1 + 1, hi⟩ := by
    simpa [hindex] using hright
  simpa [i] using hC hDense i hi hleft' hright'

/-- Empty-prefix boundary case.

If `|A ∩ [0,t]| = 0`, then `t` lies before the first ordered element of `A`.
Balasubramanian–Dutta at the first element therefore bounds `t` by one `√n`
step plus the published error term.
-/
theorem dense_sidon_empty_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L) {t : ℕ},
        (intervalSlice A 0 t).card = 0 →
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense t hzero hslt
  have hApos : 0 < A.card := by
    simpa [hzero] using hslt
  let i0 : Fin A.card := ⟨0, hApos⟩
  have hfirst_gt : t < orderedElement A i0 := by
    by_contra hnot
    have hle : orderedElement A i0 ≤ t := Nat.le_of_not_gt hnot
    have hmem : orderedElement A i0 ∈ intervalSlice A 0 t := by
      simp [intervalSlice, orderedElement_mem, hle]
    have hpos : 0 < (intervalSlice A 0 t).card := Finset.card_pos.mpr ⟨orderedElement A i0, hmem⟩
    omega
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hfirst_bd : |(orderedElement A i0 : ℝ) - s| ≤ E := by
    simpa [i0, s, E] using hC hDense i0
  have hfirst_le : (orderedElement A i0 : ℝ) ≤ s + E := by
    have hupp := (abs_le.mp hfirst_bd).2
    nlinarith
  have ht_le : (t : ℝ) ≤ s + E := by
    have hfirst_gt_real : (t : ℝ) < orderedElement A i0 := by
      exact_mod_cast hfirst_gt
    nlinarith
  have ht_nonneg : 0 ≤ (t : ℝ) := by
    exact_mod_cast Nat.zero_le t
  have habs :
      |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * s| = (t : ℝ) := by
    simp [hzero, s, abs_of_nonneg ht_nonneg]
  have hmain : |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * s| ≤ s + E := by
    rw [habs]
    exact ht_le
  simpa [s, E, add_assoc, add_left_comm, add_comm] using hmain

/-- Honest prefix theorem away from the terminal regime.

This packages the internal case `0 < |A ∩ [0,t]| < |A|` together with the
boundary case `|A ∩ [0,t]| = 0`. The only excluded regime is the terminal one
`|A ∩ [0,t]| = |A|`, where the literature-scale `m · √n` profile is not the
right center without further bookkeeping.
-/
theorem dense_sidon_nonterminal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L) {t : ℕ},
        (intervalSlice A 0 t).card < A.card →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_internal_prefix_external with ⟨Cint, hCint, hInt⟩
  rcases dense_sidon_empty_prefix_external with ⟨Cempty, hCempty, hEmpty⟩
  let C : ℝ := Cint + Cempty
  refine ⟨C, add_nonneg hCint hCempty, ?_⟩
  intro A n L hDense t hslt
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hu : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hv : 0 ≤ v := by
    dsimp [v]
    exact mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _)
  have hmono_int :
      Real.sqrt (n : ℝ) + Cint * u + Cint * v ≤
        Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    nlinarith
  have hmono_empty :
      Real.sqrt (n : ℝ) + Cempty * u + Cempty * v ≤
        Real.sqrt (n : ℝ) + C * u + C * v := by
    dsimp [C]
    nlinarith
  by_cases hs0 : 0 < (intervalSlice A 0 t).card
  · have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cint * u + Cint * v := by
      simpa [u, v, mul_assoc] using hInt hDense hs0 hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_int
    simpa [u, v, mul_assoc] using hmain
  · have hzero : (intervalSlice A 0 t).card = 0 := Nat.eq_zero_of_not_pos hs0
    have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cempty * u + Cempty * v := by
      simpa [u, v, mul_assoc] using hEmpty hDense hzero hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + C * u + C * v :=
      le_trans hbase hmono_empty
    simpa [u, v, mul_assoc] using hmain

/-- Terminal-prefix theorem with the correct deficiency drift.

If the prefix has already captured all of `A`, then the naive center
`|A ∩ [0,t]| · √n = |A| · √n` is still usable, but only after paying the
natural terminal drift `((L+1) : ℝ) · √n`. This drift comes from the floor
relation `|A| + L = Nat.sqrt n` together with the fact that `t` can range all
the way up to `n`.
-/
theorem dense_sidon_terminal_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L) {t : ℕ},
        0 < A.card →
        (intervalSlice A 0 t).card = A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense t hApos hfull ht
  have hDense' := hDense
  rcases hDense with ⟨_hSidon, _hRange, hdef⟩
  let iLast : Fin A.card := ⟨A.card - 1, by omega⟩
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hEnonneg : 0 ≤ E := by
    dsimp [E]
    exact add_nonneg
      (mul_nonneg hCnonneg (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _))
      (mul_nonneg (mul_nonneg hCnonneg (Real.sqrt_nonneg _))
        (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _))
  have hlast_raw := hC hDense' iLast
  have hlast_est : |(orderedElement A iLast : ℝ) - (A.card : ℝ) * s| ≤ E := by
    have hcard_nat : iLast.1 + 1 = A.card := by
      dsimp [iLast]
      omega
    have hcard : (((iLast.1 + 1 : ℕ) : ℝ)) = (A.card : ℝ) := by
      exact_mod_cast hcard_nat
    simpa [s, E, hcard, mul_assoc] using hlast_raw
  have hlast_le_t : orderedElement A iLast ≤ t := by
    by_contra hnot
    have hgt : t < orderedElement A iLast := Nat.lt_of_not_ge hnot
    have hsubset : intervalSlice A 0 t ⊆ A := by
      intro x hx
      simp [intervalSlice] at hx
      exact hx.1
    have hmemLast : orderedElement A iLast ∈ A := orderedElement_mem A iLast
    have hnotmemLast : orderedElement A iLast ∉ intervalSlice A 0 t := by
      simp [intervalSlice, orderedElement_mem, Nat.not_le_of_lt hgt]
    have hssub : intervalSlice A 0 t ⊂ A := by
      exact Finset.ssubset_iff_subset_ne.mpr ⟨hsubset, by
        intro heq
        exact hnotmemLast (heq.symm ▸ hmemLast)⟩
    have hcard_lt := Finset.card_lt_card hssub
    omega
  have hdrift :
      (n : ℝ) - (A.card : ℝ) * s ≤ ((L + 1 : ℕ) : ℝ) * s := by
    have hnat' : (n : ℝ) <
        (((Nat.sqrt n).succ : ℕ) : ℝ) * (((Nat.sqrt n).succ : ℕ) : ℝ) := by
      exact_mod_cast Nat.lt_succ_sqrt n
    have hnat : (n : ℝ) < (((Nat.sqrt n + 1 : ℕ) : ℝ)) ^ 2 := by
      simpa [pow_two, Nat.succ_eq_add_one, add_comm, add_left_comm, add_assoc] using hnat'
    have hs_lt : s < ((Nat.sqrt n + 1 : ℕ) : ℝ) := by
      rw [Real.sqrt_lt' (by positivity)]
      simpa [s] using hnat
    have hsq : (n : ℝ) = s ^ 2 := by
      dsimp [s]
      symm
      exact Real.sq_sqrt (by positivity)
    have hmain : (n : ℝ) ≤ ((Nat.sqrt n + 1 : ℕ) : ℝ) * s := by
      rw [hsq]
      nlinarith
    have hdef_cast : (Nat.sqrt n : ℝ) = (A.card : ℝ) + (L : ℝ) := by
      exact_mod_cast hdef.symm
    have hdef_cast' : (((Nat.sqrt n + 1 : ℕ) : ℝ)) = (A.card : ℝ) + (L : ℝ) + 1 := by
      rw [Nat.cast_add, Nat.cast_one, hdef_cast]
    rw [hdef_cast'] at hmain
    have hsub := sub_le_sub_right hmain ((A.card : ℝ) * s)
    simpa [Nat.cast_add, Nat.cast_one, add_mul, mul_add, add_assoc, add_left_comm, add_comm] using hsub
  have hdrift_nonneg : 0 ≤ ((L + 1 : ℕ) : ℝ) * s := by
    exact mul_nonneg (by positivity) hs_nonneg
  have ht_le_real : (t : ℝ) ≤ n := by
    exact_mod_cast ht
  have hupper0 : (t : ℝ) - (A.card : ℝ) * s ≤ ((L + 1 : ℕ) : ℝ) * s := by
    nlinarith [ht_le_real, hdrift]
  have hlast_le_real : (orderedElement A iLast : ℝ) ≤ t := by
    exact_mod_cast hlast_le_t
  have hmono :
      (orderedElement A iLast : ℝ) - (A.card : ℝ) * s ≤
        (t : ℝ) - (A.card : ℝ) * s := by
    nlinarith
  have hlowE : -E ≤ (orderedElement A iLast : ℝ) - (A.card : ℝ) * s := by
    exact (abs_le.mp hlast_est).1
  have hlow :
      -((((L + 1 : ℕ) : ℝ) * s) + E) ≤ (t : ℝ) - (A.card : ℝ) * s := by
    nlinarith [hlowE, hmono, hdrift_nonneg]
  have hupp :
      (t : ℝ) - (A.card : ℝ) * s ≤ (((L + 1 : ℕ) : ℝ) * s) + E := by
    nlinarith [hupper0, hEnonneg]
  have habs :
      |(t : ℝ) - (A.card : ℝ) * s| ≤ (((L + 1 : ℕ) : ℝ) * s) + E := by
    exact abs_le.mpr ⟨hlow, hupp⟩
  simpa [hfull, s, E, add_assoc, add_left_comm, add_comm, mul_assoc] using habs

/-- Honest all-prefix theorem for positive-card dense Sidon sets.

This combines the nonterminal and terminal regimes into a single statement for
all prefixes `[0,t]` with `t ≤ n`, provided `A` is nonempty. The price of
uniformity is exactly the terminal deficiency drift `((L+1) : ℝ) · √n`.
-/
theorem dense_sidon_positive_card_prefix_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L) {t : ℕ},
        0 < A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) +
            C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_nonterminal_prefix_external with ⟨Cnon, hCnon, hNon⟩
  rcases dense_sidon_terminal_prefix_external with ⟨Cterm, hCterm, hTerm⟩
  let C : ℝ := Cnon + Cterm
  refine ⟨C, add_nonneg hCnon hCterm, ?_⟩
  intro A n L hDense t hApos ht
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hs_nonneg : 0 ≤ Real.sqrt (n : ℝ) := Real.sqrt_nonneg _
  have hdrift_dom :
      Real.sqrt (n : ℝ) ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) := by
    have hfac : (1 : ℝ) ≤ ((L + 1 : ℕ) : ℝ) := by
      exact_mod_cast Nat.succ_le_succ (Nat.zero_le L)
    nlinarith
  have hu : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _
  have hv : 0 ≤ v := by
    dsimp [v]
    exact mul_nonneg (Real.sqrt_nonneg _) (Real.rpow_nonneg (by exact_mod_cast Nat.zero_le n) _)
  have hmono_non :
      Real.sqrt (n : ℝ) + Cnon * u + Cnon * v ≤
        (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + C * u + C * v := by
    dsimp [C]
    nlinarith
  have hmono_term :
      (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + Cterm * u + Cterm * v ≤
        (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + C * u + C * v := by
    dsimp [C]
    nlinarith
  have hsubset : intervalSlice A 0 t ⊆ A := by
    intro x hx
    simp [intervalSlice] at hx
    exact hx.1
  have hcard_le : (intervalSlice A 0 t).card ≤ A.card := Finset.card_le_card hsubset
  by_cases hslt : (intervalSlice A 0 t).card < A.card
  · have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Real.sqrt (n : ℝ) + Cnon * u + Cnon * v := by
      simpa [u, v, mul_assoc] using hNon hDense hslt
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + C * u + C * v :=
      le_trans hbase hmono_non
    simpa [u, v, mul_assoc] using hmain
  · have hfull : (intervalSlice A 0 t).card = A.card := by
      omega
    have hbase :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + Cterm * u + Cterm * v := by
      simpa [u, v, mul_assoc] using hTerm hDense hApos hfull ht
    have hmain :
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ (((L + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)) + C * u + C * v :=
      le_trans hbase hmono_term
    simpa [u, v, mul_assoc] using hmain

/-- Ordered mass-balance theorem at the literature scale.

Summing the ordered-element theorem over every index gives global control of the
center of mass of a dense Sidon set relative to the linear `m · √n` profile.
This is the exact ordered-index form of the Balasubramanian–Dutta mass-balance
consequence, before simplifying the profile sum to a closed form.
-/
theorem dense_sidon_ordered_mass_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L),
        |(∑ i : Fin A.card, (orderedElement A i : ℝ)) -
            ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
          ≤ (A.card : ℝ) *
            (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
              C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hrewrite :
      (∑ i : Fin A.card, (orderedElement A i : ℝ)) -
          ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * s) =
        ∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)) := by
    rw [Finset.sum_sub_distrib]
  have habs :
      |∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s))|
        ≤ ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)| := by
    exact Finset.abs_sum_le_sum_abs _ _
  have hsum :
      ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)|
        ≤ ∑ _i : Fin A.card, E := by
    exact Finset.sum_le_sum (fun i _ => by simpa [s, E] using hC hDense i)
  have hcard :
      (∑ _i : Fin A.card, E) = (A.card : ℝ) * E := by
    simp [E, mul_add, mul_assoc]
  calc
    |(∑ i : Fin A.card, (orderedElement A i : ℝ)) -
        ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
      = |∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s))| := by
          rw [← hrewrite]
    _ ≤ ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)| := habs
    _ ≤ ∑ _i : Fin A.card, E := hsum
    _ = (A.card : ℝ) *
          (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
          simpa [E, mul_add, mul_assoc] using hcard

/-- Finset-mass version of the external ordered mass-balance theorem.

Using the local theorem that the ordered enumeration sums to `A.sum id`, the
ordered-index mass statement becomes a genuine statement about the total mass
of the Sidon set itself.
-/
theorem dense_sidon_finset_mass_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L),
        |(∑ a ∈ A, (a : ℝ)) -
            ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
          ≤ (A.card : ℝ) *
            (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
              C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
  rcases dense_sidon_ordered_mass_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense
  have hsum_nat : (∑ i : Fin A.card, orderedElement A i) = A.sum (fun a => a) :=
    sum_orderedElement_eq_sum A
  have hsum :
      (∑ i : Fin A.card, (orderedElement A i : ℝ)) = ∑ a ∈ A, (a : ℝ) := by
    calc
      (∑ i : Fin A.card, (orderedElement A i : ℝ))
          = ((∑ i : Fin A.card, orderedElement A i : ℕ) : ℝ) := by
              rw [← Nat.cast_sum]
      _ = ∑ a ∈ A, (a : ℝ) := by
            rw [hsum_nat, Nat.cast_sum]
  simpa [hsum] using hC hDense

/-- Closed form for the linear profile sum, without the final triangular-number
compression.

This is the readable arithmetic center that follows from `Finset.sum_range_id`
without spending extra proof effort on the identity
`k * (k - 1) / 2 + k = k * (k + 1) / 2`.
-/
theorem dense_sidon_profile_sum_explicit (k : ℕ) (s : ℝ) :
    (∑ i : Fin k, (((i.1 + 1 : ℕ) : ℝ) * s)) =
      (((k * (k - 1) / 2 + k : ℕ) : ℝ) * s) := by
  rw [← Finset.sum_mul]
  congr 1
  have hnat : ∑ x ∈ Finset.range k, (x + 1 : ℕ) = k * (k - 1) / 2 + k := by
    calc
      ∑ x ∈ Finset.range k, (x + 1 : ℕ)
        = (∑ x ∈ Finset.range k, x) + (∑ _x ∈ Finset.range k, (1 : ℕ)) := by
            rw [Finset.sum_add_distrib]
      _ = k * (k - 1) / 2 + k := by simp [Finset.sum_range_id]
  rw [Fin.sum_univ_eq_sum_range (fun m => (((m + 1 : ℕ) : ℝ))) k]
  rw [← Nat.cast_sum]
  exact_mod_cast hnat

/-- Finset-mass theorem with an explicit arithmetic center.

This is the same theorem as `dense_sidon_finset_mass_external`, but the profile
sum has been rewritten to the concrete center
`((|A| * (|A| - 1) / 2) + |A|) · √n`.
-/
theorem dense_sidon_finset_mass_explicit_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L),
        |(∑ a ∈ A, (a : ℝ)) -
            (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))|
          ≤ (A.card : ℝ) *
            (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
              C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
  rcases dense_sidon_finset_mass_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense
  have hprofile := dense_sidon_profile_sum_explicit A.card (Real.sqrt (n : ℝ))
  have hprofile' :
      (∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))) =
        ((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)) := by
    simpa [Nat.cast_add] using hprofile
  have hmain := hC hDense
  rw [hprofile'] at hmain
  exact hmain

/-- Maximizer-friendly finset-mass theorem with the explicit arithmetic center.

This is the same mass-balance consequence as
`dense_sidon_finset_mass_explicit_external`, but parameterized by the
literature deficiency `max(0, sqrt n - |A|)` so it still applies to true
maximizers when they lie above `floor (sqrt n)`.
-/
theorem sidon_in_range_finset_mass_explicit_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n),
        |(∑ a ∈ A, (a : ℝ)) -
            (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))|
          ≤ (A.card : ℝ) *
            (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
              C * Real.sqrt (realDeficiencyFromSqrt A n) *
                Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
  rcases sidon_in_range_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n hA
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
    C * Real.sqrt (realDeficiencyFromSqrt A n) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hrewrite :
      (∑ i : Fin A.card, (orderedElement A i : ℝ)) -
          ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * s) =
        ∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)) := by
    rw [Finset.sum_sub_distrib]
  have habs :
      |∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s))|
        ≤ ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)| := by
    exact Finset.abs_sum_le_sum_abs _ _
  have hsum :
      ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)|
        ≤ ∑ _i : Fin A.card, E := by
    exact Finset.sum_le_sum (fun i _ => by simpa [s, E] using hC hA i)
  have hcard :
      (∑ _i : Fin A.card, E) = (A.card : ℝ) * E := by
    simp [E, mul_add, mul_assoc]
  have hsum_nat : (∑ i : Fin A.card, orderedElement A i) = A.sum (fun a => a) :=
    sum_orderedElement_eq_sum A
  have hsum_real :
      (∑ i : Fin A.card, (orderedElement A i : ℝ)) = ∑ a ∈ A, (a : ℝ) := by
    calc
      (∑ i : Fin A.card, (orderedElement A i : ℝ))
          = ((∑ i : Fin A.card, orderedElement A i : ℕ) : ℝ) := by
              rw [← Nat.cast_sum]
      _ = ∑ a ∈ A, (a : ℝ) := by
            rw [hsum_nat, Nat.cast_sum]
  have hprofile := dense_sidon_profile_sum_explicit A.card (Real.sqrt (n : ℝ))
  have hprofile' :
      (∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))) =
        ((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)) := by
    simpa [Nat.cast_add] using hprofile
  calc
    |(∑ a ∈ A, (a : ℝ)) -
        (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))|
      = |(∑ i : Fin A.card, (orderedElement A i : ℝ)) -
          ∑ i : Fin A.card, (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))| := by
          rw [← hsum_real, hprofile']
    _ = |∑ i : Fin A.card, ((orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s))| := by
          rw [← hrewrite]
    _ ≤ ∑ i : Fin A.card, |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * s)| := habs
    _ ≤ ∑ _i : Fin A.card, E := hsum
    _ = (A.card : ℝ) *
          (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (realDeficiencyFromSqrt A n) *
              Real.rpow (n : ℝ) (3 / 4 : ℝ)) := by
          simpa [E, mul_add, mul_assoc] using hcard

/-- Super-floor coarse mass theorem.

On the super-floor corridor, the maximizer-friendly mass theorem collapses to a
single `n^(11/8)` scale. The extra factor of `n^(1/2)` reflects the fact that
mass sums over `|A| ≍ sqrt(n)` ordered elements rather than controlling one
prefix or one cutpoint.
-/
theorem sidon_in_range_superfloor_finset_mass_coarse_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        |(∑ a ∈ A, (a : ℝ)) -
            (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))|
          ≤ C * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
  rcases sidon_in_range_finset_mass_explicit_external with ⟨C, hCnonneg, hC⟩
  let C' : ℝ := 6 * C
  refine ⟨C', by positivity, ?_⟩
  intro A n hA hn hfloor
  set s : ℝ := Real.sqrt (n : ℝ)
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set v : ℝ := Real.rpow (n : ℝ) (3 / 4 : ℝ)
  set w : ℝ := Real.rpow (n : ℝ) (11 / 8 : ℝ)
  have hn0 : 0 ≤ (n : ℝ) := by exact_mod_cast Nat.zero_le n
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hu_nonneg : 0 ≤ u := by
    dsimp [u]
    exact Real.rpow_nonneg hn0 _
  have hv_nonneg : 0 ≤ v := by
    dsimp [v]
    exact Real.rpow_nonneg hn0 _
  have hw_nonneg : 0 ≤ w := by
    dsimp [w]
    exact Real.rpow_nonneg hn0 _
  have hsqrt_nat :
      (Nat.sqrt n : ℝ) ≤ s := by
    dsimp [s]
    exact Real.nat_sqrt_le_real_sqrt
  have hquarter_nat :
      (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ s := by
    have h1 : (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
      have hnat : (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.sqrt (Nat.sqrt n : ℝ) :=
        Real.nat_sqrt_le_real_sqrt
      have hmono : Real.sqrt (Nat.sqrt n : ℝ) ≤ Real.sqrt (Real.sqrt (n : ℝ)) := by
        apply Real.sqrt_le_sqrt
        exact_mod_cast Real.nat_sqrt_le_real_sqrt
      have hpow : Real.sqrt (Real.sqrt (n : ℝ)) = Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
        rw [Real.sqrt_eq_rpow, Real.sqrt_eq_rpow, ← Real.rpow_mul hn0]
        norm_num
      exact le_trans hnat (hmono.trans_eq hpow)
    have h2 : Real.rpow (n : ℝ) (1 / 4 : ℝ) ≤ s := by
      dsimp [s]
      rw [Real.sqrt_eq_rpow]
      exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
    exact le_trans h1 h2
  have hsqrt_one : (1 : ℝ) ≤ s := by
    have hsq : (n : ℝ) = s ^ 2 := by
      dsimp [s]
      symm
      exact Real.sq_sqrt (by positivity)
    nlinarith
  have hL := lindstrom_bound A n hA.1 hA.2 hn
  have hcard_hi :
      (A.card : ℝ) ≤ (Nat.sqrt n : ℝ) + (Nat.sqrt (Nat.sqrt n) : ℝ) + 1 := by
    exact_mod_cast hL
  have hcard_s :
      (A.card : ℝ) ≤ 3 * s := by
    nlinarith [hcard_hi, hsqrt_nat, hquarter_nat, hsqrt_one]
  have hs_mul_u :
      s * u = w := by
    dsimp [s, u, w]
    rw [Real.sqrt_eq_rpow]
    rw [← Real.rpow_add_of_nonneg hn0 (by positivity : 0 ≤ (1 / 2 : ℝ))
      (by positivity : 0 ≤ (7 / 8 : ℝ))]
    congr 1
    norm_num
  have hv_le_u :
      v ≤ u := by
    dsimp [u, v]
    exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
  have hs_mul_v_le :
      s * v ≤ w := by
    have hmul : s * v ≤ s * u := by
      exact mul_le_mul_of_nonneg_left hv_le_u hs_nonneg
    rw [hs_mul_u] at hmul
    exact hmul
  have hdef : realDeficiencyFromSqrt A n ≤ 1 :=
    realDeficiencyFromSqrt_le_one_of_superfloor hA hfloor
  have hsqrtdef : Real.sqrt (realDeficiencyFromSqrt A n) ≤ 1 := by
    rw [Real.sqrt_le_iff]
    constructor
    · positivity
    · simpa using hdef
  have hAu :
      (A.card : ℝ) * u ≤ 3 * w := by
    calc
      (A.card : ℝ) * u ≤ (3 * s) * u := by
        exact mul_le_mul_of_nonneg_right hcard_s hu_nonneg
      _ = 3 * w := by rw [mul_assoc, hs_mul_u]
  have hAv :
      (A.card : ℝ) * Real.sqrt (realDeficiencyFromSqrt A n) * v ≤ 3 * w := by
    calc
      (A.card : ℝ) * Real.sqrt (realDeficiencyFromSqrt A n) * v
          ≤ (A.card : ℝ) * 1 * v := by
            gcongr
      _ = (A.card : ℝ) * v := by ring
      _ ≤ (3 * s) * v := by
            exact mul_le_mul_of_nonneg_right hcard_s hv_nonneg
      _ ≤ 3 * w := by
            nlinarith [hs_mul_v_le]
  calc
    |(∑ a ∈ A, (a : ℝ)) -
        (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))|
        ≤ (A.card : ℝ) *
            (C * u + C * Real.sqrt (realDeficiencyFromSqrt A n) * v) := by
              simpa [u, v, s, mul_assoc] using hC hA
    _ = C * ((A.card : ℝ) * u) +
          C * ((A.card : ℝ) * Real.sqrt (realDeficiencyFromSqrt A n) * v) := by
          ring
    _ ≤ C * (3 * w) + C * (3 * w) := by
          gcongr
    _ = C' * w := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
          rfl

/-- Super-floor coarse mass theorem with density-adjusted center.

The exact maximizer packet shows that the mass side is better centered by the
ambient slope `n / |A|` than by the raw `sqrt(n)` profile. On the super-floor
corridor, the difference between these two centers is itself absorbed at the
same coarse `n^(11/8)` scale, so the theorem can be repackaged around the
readable density-adjusted center `n (|A| + 1) / 2`.
-/
theorem sidon_in_range_superfloor_finset_mass_density_adjusted_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        |(∑ a ∈ A, (a : ℝ)) - ((n : ℝ) * ((A.card : ℝ) + 1) / 2)|
          ≤ C * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_finset_mass_coarse_external with
    ⟨C, hCnonneg, hC⟩
  let C' : ℝ := C + 4
  refine ⟨C', by linarith, ?_⟩
  intro A n hA hn hfloor
  set k : ℝ := (A.card : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set u : ℝ := Real.rpow (n : ℝ) (7 / 8 : ℝ)
  set w : ℝ := Real.rpow (n : ℝ) (11 / 8 : ℝ)
  set q : ℝ := (Nat.sqrt (Nat.sqrt n) : ℝ) + 1
  set oldCenter : ℝ :=
    (((↑(A.card * (A.card - 1) / 2) + ↑A.card) * Real.sqrt (n : ℝ)))
  set newCenter : ℝ := ((n : ℝ) * ((A.card : ℝ) + 1) / 2)
  have hn0 : 0 ≤ (n : ℝ) := by exact_mod_cast Nat.zero_le n
  have hn1 : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have hs_nonneg : 0 ≤ s := by
    dsimp [s]
    exact Real.sqrt_nonneg _
  have hsqrt_one : (1 : ℝ) ≤ s := by
    have hsq : (n : ℝ) = s ^ 2 := by
      dsimp [s]
      symm
      exact Real.sq_sqrt (by positivity)
    nlinarith
  have hs_mul_u :
      s * u = w := by
    dsimp [s, u, w]
    rw [Real.sqrt_eq_rpow]
    rw [← Real.rpow_add_of_nonneg hn0 (by positivity : 0 ≤ (1 / 2 : ℝ))
      (by positivity : 0 ≤ (7 / 8 : ℝ))]
    congr 1
    norm_num
  have hmain :
      |(∑ a ∈ A, (a : ℝ)) - oldCenter| ≤ C * w := by
    simpa [oldCenter, w] using hC hA hn hfloor
  have hshift :
      |oldCenter - newCenter| ≤ 4 * w := by
    have hL := lindstrom_bound A n hA.1 hA.2 hn
    have hcard_hi :
        k ≤ (Nat.sqrt n : ℝ) + (Nat.sqrt (Nat.sqrt n) : ℝ) + 1 := by
      dsimp [k]
      exact_mod_cast hL
    have hsqrt_nat :
        (Nat.sqrt n : ℝ) ≤ s := by
      dsimp [s]
      exact Real.nat_sqrt_le_real_sqrt
    have hquarter_nat :
        (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ s := by
      have h1 : (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
        have hnat : (Nat.sqrt (Nat.sqrt n) : ℝ) ≤ Real.sqrt (Nat.sqrt n : ℝ) :=
          Real.nat_sqrt_le_real_sqrt
        have hmono : Real.sqrt (Nat.sqrt n : ℝ) ≤ Real.sqrt (Real.sqrt (n : ℝ)) := by
          apply Real.sqrt_le_sqrt
          exact_mod_cast Real.nat_sqrt_le_real_sqrt
        have hpow : Real.sqrt (Real.sqrt (n : ℝ)) = Real.rpow (n : ℝ) (1 / 4 : ℝ) := by
          rw [Real.sqrt_eq_rpow, Real.sqrt_eq_rpow, ← Real.rpow_mul hn0]
          norm_num
        exact le_trans hnat (hmono.trans_eq hpow)
      have h2 : Real.rpow (n : ℝ) (1 / 4 : ℝ) ≤ s := by
        dsimp [s]
        rw [Real.sqrt_eq_rpow]
        exact Real.rpow_le_rpow_of_exponent_le hn1 (by norm_num)
      exact le_trans h1 h2
    have hcard_s :
        k ≤ 3 * s := by
      nlinarith [hcard_hi, hsqrt_nat, hquarter_nat, hsqrt_one]
    have hgap :
        realGapFromSqrt A n ≤ q := by
      dsimp [q]
      exact realGapFromSqrt_le_fourthRoot_of_superfloor hA hn hfloor
    have hambient :
        q * s ≤ 2 * u := by
      simpa [q, s, u] using superfloor_ambient_le_two_sevenEighths hn
    have hk1 :
        k + 1 ≤ 4 * s := by
      nlinarith [hcard_s, hsqrt_one]
    have hold_rewrite :
        oldCenter = ((k * (k + 1) / 2) * s) := by
      dsimp [oldCenter, k, s]
      have htri : A.card * (A.card - 1) / 2 + A.card = A.card * (A.card + 1) / 2 := by
        have hleft_even : 2 ∣ A.card * (A.card - 1) :=
          even_iff_two_dvd.mp (Nat.even_mul_pred_self A.card)
        have hright_even : 2 ∣ A.card * (A.card + 1) :=
          even_iff_two_dvd.mp (Nat.even_mul_succ_self A.card)
        apply Nat.eq_of_mul_eq_mul_left (by norm_num : 0 < 2)
        calc
          2 * (A.card * (A.card - 1) / 2 + A.card)
              = A.card * (A.card - 1) + 2 * A.card := by
                  rw [Nat.mul_add, Nat.mul_div_cancel' hleft_even]
          _ = A.card * (A.card + 1) := by
                cases' A.card with m
                · simp
                · simp
                  ring_nf
          _ = 2 * (A.card * (A.card + 1) / 2) := by
                rw [Nat.mul_div_cancel' hright_even]
      rw [show ((↑(A.card * (A.card - 1) / 2) + ↑A.card) : ℝ) =
          (↑(A.card * (A.card + 1) / 2) : ℝ) by
            rw [← Nat.cast_add, htri]]
      have hdiv : (↑(A.card * (A.card + 1) / 2) : ℝ) = k * (k + 1) / 2 := by
        dsimp [k]
        rw [Nat.cast_div]
        · norm_num
        · exact even_iff_two_dvd.mp (Nat.even_mul_succ_self A.card)
        · norm_num
      rw [hdiv]
    have hnew_rewrite :
        newCenter = ((k + 1) / 2) * s * s := by
      dsimp [newCenter, k]
      have hsq : (n : ℝ) = s * s := by
        dsimp [s]
        rw [← sq]
        symm
        exact Real.sq_sqrt hn0
      rw [hsq]
      ring
    have hshift_formula :
        oldCenter - newCenter = (((k + 1) / 2) * s) * (k - s) := by
      rw [hold_rewrite, hnew_rewrite]
      ring
    rw [hshift_formula, abs_mul]
    have hs_term :
        |((k + 1) / 2) * s| ≤ 2 * s * s := by
      have hnonneg : 0 ≤ ((k + 1) / 2) * s := by
        positivity
      rw [abs_of_nonneg hnonneg]
      nlinarith [hk1, hs_nonneg]
    have hgap' : |k - s| ≤ q := by
      dsimp [k, s]
      simpa [realGapFromSqrt, abs_sub_comm] using hgap
    have hqss : q * (s * s) ≤ 2 * w := by
      have hstep : q * (s * s) = (q * s) * s := by ring
      rw [hstep]
      calc
        (q * s) * s ≤ (2 * u) * s := by
          gcongr
        _ = 2 * w := by
          calc
            (2 * u) * s = 2 * (u * s) := by ring
            _ = 2 * (s * u) := by ring
            _ = 2 * w := by rw [hs_mul_u]
    calc
      |((k + 1) / 2) * s| * |k - s|
          ≤ (2 * s * s) * q := by
            gcongr
      _ = 2 * (q * (s * s)) := by ring
      _ ≤ 2 * (2 * w) := by gcongr
      _ = 4 * w := by ring
  calc
    |(∑ a ∈ A, (a : ℝ)) - newCenter|
        ≤ |(∑ a ∈ A, (a : ℝ)) - oldCenter| + |oldCenter - newCenter| := by
          simpa using abs_sub_le (∑ a ∈ A, (a : ℝ)) oldCenter newCenter
    _ ≤ C * w + 4 * w := by
          nlinarith [hmain, hshift]
    _ = C' * w := by
          dsimp [C']
          ring
    _ = C' * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
          rfl

/-- Joint super-floor envelope theorem for the prefix and density-adjusted mass
observables.

This is the nearest honest formal shadow of the exact observable-split story:
the same Sidon set simultaneously satisfies a coarse prefix envelope at
`n^(7/8)` scale and a density-adjusted mass envelope at `n^(11/8)` scale. It
does not claim a tradeoff or incompatibility between the two observables. It
only packages the fact that the current Lean substrate can see both readouts on
the same object at once.
-/
theorem sidon_in_range_superfloor_prefix_mass_joint_envelope_external :
    ∃ Cpref Cmass : ℝ, 0 ≤ Cpref ∧ 0 ≤ Cmass ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : SidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Cpref * Real.rpow (n : ℝ) (7 / 8 : ℝ) ∧
        |(∑ a ∈ A, (a : ℝ)) - ((n : ℝ) * ((A.card : ℝ) + 1) / 2)|
          ≤ Cmass * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_prefix_coarse_external' with
    ⟨Cpref, hCpref_nonneg, hPref⟩
  rcases sidon_in_range_superfloor_finset_mass_density_adjusted_external with
    ⟨Cmass, hCmass_nonneg, hMass⟩
  refine ⟨Cpref, Cmass, hCpref_nonneg, hCmass_nonneg, ?_⟩
  intro A n t hA hn hfloor ht
  exact ⟨hPref hA hn hfloor ht, hMass hA hn hfloor⟩

/-- Exact-extremal-surface wrapper for the joint prefix/mass envelope.

This theorem does not prove the finite Pareto or compatibility phenomenon. It
only connects the compiled joint envelope to the new exact-extremal vocabulary
so future statements can be conditioned on `IsMaximalSidonInRange` without
reopening the prefix and mass proof packages.
-/
theorem maximal_sidon_in_range_superfloor_prefix_mass_joint_envelope_external :
    ∃ Cpref Cmass : ℝ, 0 ≤ Cpref ∧ 0 ≤ Cmass ∧
      ∀ {A : Finset ℕ} {n t : ℕ} (hA : IsMaximalSidonInRange A n),
        0 < n →
        Nat.sqrt n ≤ A.card →
        t ≤ n →
        |(t : ℝ) - ((intervalSlice A 0 t).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ Cpref * Real.rpow (n : ℝ) (7 / 8 : ℝ) ∧
        densityAdjustedMassDeviation A n
          ≤ Cmass * Real.rpow (n : ℝ) (11 / 8 : ℝ) := by
  rcases sidon_in_range_superfloor_prefix_mass_joint_envelope_external with
    ⟨Cpref, Cmass, hCpref_nonneg, hCmass_nonneg, hJoint⟩
  refine ⟨Cpref, Cmass, hCpref_nonneg, hCmass_nonneg, ?_⟩
  intro A n t hA hn hfloor ht
  simpa [densityAdjustedMassDeviation] using hJoint hA.sidonInRange hn hfloor ht

end Erdos.Sidon
