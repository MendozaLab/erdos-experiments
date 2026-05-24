/-
  Erdős #30 / BFR representation bridge scout
  ===========================================

  Review-only scout file.  This file starts the Block A bridge from
  `CLOSURE_PLAN_BFR.md` without modifying the canonical `Erdos30_BFR.lean`
  artifact or promoting any theorem status.

  Purpose:
  * define the ordered representation function r_A(s),
  * pin down the finite support range used by the BFR interval partition,
  * add small reusable lemmas that are independent of the constant-bearing
    BFR inequality.

  Claim ceiling: local formalization scout only.  It does not close
  `bfr_core_bound`, improve the public Sidon coefficient, or alter any
  registry/publication state.
-/

import Mathlib
import Erdos30_Sidon_Defs

open Finset Nat

namespace Erdos.Sidon

/-! ## Ordered representation function -/

/-- Ordered representation count `r_A(s) = #{(a,b) in A x A | a + b = s}`. -/
def bfrRepFunction (A : Finset ℕ) (s : ℕ) : ℕ :=
  ((A ×ˢ A).filter (fun p : ℕ × ℕ => p.1 + p.2 = s)).card

/-- The BFR representation function is bounded by the ordered square. -/
theorem bfrRepFunction_le_card_product (A : Finset ℕ) (s : ℕ) :
    bfrRepFunction A s ≤ (A ×ˢ A).card := by
  unfold bfrRepFunction
  exact Finset.card_filter_le _ _

/-- The ordered square has cardinality `|A|^2`. -/
theorem bfrRepFunction_le_card_sq (A : Finset ℕ) (s : ℕ) :
    bfrRepFunction A s ≤ A.card * A.card := by
  calc
    bfrRepFunction A s ≤ (A ×ˢ A).card := bfrRepFunction_le_card_product A s
    _ = A.card * A.card := by rw [Finset.card_product]

/-- The ordered representations with `a <= b` are unique for a Sidon set. -/
theorem bfrRepFunction_le_half_card_le_one
    (A : Finset ℕ) (s : ℕ) (hS : IsSidonSet A) :
    ((A ×ˢ A).filter
      (fun p : ℕ × ℕ => p.1 + p.2 = s ∧ p.1 ≤ p.2)).card ≤ 1 := by
  rw [Finset.card_le_one_iff]
  intro p q hp hq
  simp only [Finset.mem_filter, Finset.mem_product] at hp hq
  rcases hp with ⟨⟨hp₁, hp₂⟩, hpsum, hple⟩
  rcases hq with ⟨⟨hq₁, hq₂⟩, hqsum, hqle⟩
  have hsum : p.1 + p.2 = q.1 + q.2 := by omega
  have h := hS p.1 hp₁ p.2 hp₂ q.1 hq₁ q.2 hq₂ hple hqle hsum
  exact Prod.ext h.1 h.2

/-- The ordered representations with `b < a` are unique for a Sidon set. -/
theorem bfrRepFunction_gt_half_card_le_one
    (A : Finset ℕ) (s : ℕ) (hS : IsSidonSet A) :
    ((A ×ˢ A).filter
      (fun p : ℕ × ℕ => p.1 + p.2 = s ∧ p.2 < p.1)).card ≤ 1 := by
  rw [Finset.card_le_one_iff]
  intro p q hp hq
  simp only [Finset.mem_filter, Finset.mem_product] at hp hq
  rcases hp with ⟨⟨hp₁, hp₂⟩, hpsum, hplt⟩
  rcases hq with ⟨⟨hq₁, hq₂⟩, hqsum, hqlt⟩
  have hsum : p.2 + p.1 = q.2 + q.1 := by omega
  have h := hS p.2 hp₂ p.1 hp₁ q.2 hq₂ q.1 hq₁
    (le_of_lt hplt) (le_of_lt hqlt) hsum
  exact Prod.ext h.2 h.1

/-- For a Sidon set, each integer has at most two ordered representations as `a + b`. -/
theorem bfrRepFunction_sidon_le_two
    (A : Finset ℕ) (s : ℕ) (hS : IsSidonSet A) :
    bfrRepFunction A s ≤ 2 := by
  let F := (A ×ˢ A).filter (fun p : ℕ × ℕ => p.1 + p.2 = s)
  let Fle := (A ×ˢ A).filter
    (fun p : ℕ × ℕ => p.1 + p.2 = s ∧ p.1 ≤ p.2)
  let Fgt := (A ×ˢ A).filter
    (fun p : ℕ × ℕ => p.1 + p.2 = s ∧ p.2 < p.1)
  have h_union : F = Fle ∪ Fgt := by
    ext p
    simp only [F, Fle, Fgt, Finset.mem_filter, Finset.mem_product, Finset.mem_union]
    constructor
    · intro hp
      rcases hp with ⟨hpA, hsum⟩
      rcases le_or_gt p.1 p.2 with hle | hgt
      · exact Or.inl ⟨hpA, hsum, hle⟩
      · exact Or.inr ⟨hpA, hsum, hgt⟩
    · intro hp
      rcases hp with hp | hp
      · exact ⟨hp.1, hp.2.1⟩
      · exact ⟨hp.1, hp.2.1⟩
  have h_disj : Disjoint Fle Fgt := by
    rw [Finset.disjoint_iff_ne]
    intro p hp q hq hpq
    subst q
    simp only [Fle, Fgt, Finset.mem_filter, Finset.mem_product] at hp hq
    exact (not_lt.mpr hp.2.2) hq.2.2
  have hle : Fle.card ≤ 1 := bfrRepFunction_le_half_card_le_one A s hS
  have hgt : Fgt.card ≤ 1 := bfrRepFunction_gt_half_card_le_one A s hS
  calc
    bfrRepFunction A s = F.card := rfl
    _ = (Fle ∪ Fgt).card := by rw [h_union]
    _ = Fle.card + Fgt.card := Finset.card_union_of_disjoint h_disj
    _ ≤ 1 + 1 := Nat.add_le_add hle hgt
    _ = 2 := rfl

/-- If `A` is contained in `{0,...,N}`, then every represented sum is below `2N + 1`. -/
theorem bfrRepFunction_sum_mem_range_of_pair
    (A : Finset ℕ) (N s : ℕ)
    (hA : A ⊆ Finset.range (N + 1))
    {p : ℕ × ℕ}
    (hp : p ∈ (A ×ˢ A).filter (fun q : ℕ × ℕ => q.1 + q.2 = s)) :
    s ∈ Finset.range (2 * N + 1) := by
  simp only [Finset.mem_filter, Finset.mem_product] at hp
  rcases hp with ⟨⟨ha, hb⟩, hsum⟩
  have haN : p.1 ≤ N := by
    have := Finset.mem_range.mp (hA ha)
    omega
  have hbN : p.2 ≤ N := by
    have := Finset.mem_range.mp (hA hb)
    omega
  exact Finset.mem_range.mpr (by omega)

/-- No ordered representations exist outside the natural BFR sum range. -/
theorem bfrRepFunction_eq_zero_of_not_mem_range
    (A : Finset ℕ) (N s : ℕ)
    (hA : A ⊆ Finset.range (N + 1))
    (hs : s ∉ Finset.range (2 * N + 1)) :
    bfrRepFunction A s = 0 := by
  unfold bfrRepFunction
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro p hp hsum
  have hmem : p ∈ (A ×ˢ A).filter (fun q : ℕ × ℕ => q.1 + q.2 = s) :=
    Finset.mem_filter.mpr ⟨hp, hsum⟩
  exact hs (bfrRepFunction_sum_mem_range_of_pair A N s hA hmem)

/-! ## Explicit BFR support range -/

/-- The finite support range used for BFR ordered sums when `A ⊆ {0,...,N}`. -/
def bfrSumRange (N : ℕ) : Finset ℕ :=
  Finset.range (2 * N + 1)

theorem mem_bfrSumRange_iff (N s : ℕ) :
    s ∈ bfrSumRange N ↔ s ≤ 2 * N := by
  unfold bfrSumRange
  constructor
  · intro h
    have := Finset.mem_range.mp h
    omega
  · intro h
    exact Finset.mem_range.mpr (by omega)

theorem bfrRepFunction_eq_zero_of_not_mem_bfrSumRange
    (A : Finset ℕ) (N s : ℕ)
    (hA : A ⊆ Finset.range (N + 1))
    (hs : s ∉ bfrSumRange N) :
    bfrRepFunction A s = 0 :=
  bfrRepFunction_eq_zero_of_not_mem_range A N s hA hs

/-! ## Fiberwise counting bridge -/

/-- Summing the ordered representation function over its support range recovers `|A|^2`. -/
theorem bfr_sum_repFunction_eq_card_sq
    (A : Finset ℕ) (N : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    ∑ s ∈ bfrSumRange N, bfrRepFunction A s = A.card * A.card := by
  have h_maps : (↑(A ×ˢ A) : Set (ℕ × ℕ)).MapsTo
      (fun p : ℕ × ℕ => p.1 + p.2) (↑(bfrSumRange N) : Set ℕ) := by
    intro p hp
    simp only [Finset.mem_coe, Finset.mem_product] at hp
    rcases hp with ⟨ha, hb⟩
    have haN : p.1 ≤ N := by
      have := Finset.mem_range.mp (hA ha)
      omega
    have hbN : p.2 ≤ N := by
      have := Finset.mem_range.mp (hA hb)
      omega
    exact (mem_bfrSumRange_iff N (p.1 + p.2)).mpr (by omega)
  calc
    ∑ s ∈ bfrSumRange N, bfrRepFunction A s
        = (A ×ˢ A).card := by
            rw [Finset.card_eq_sum_card_fiberwise h_maps]
            simp [bfrRepFunction]
    _ = A.card * A.card := by rw [Finset.card_product]

end Erdos.Sidon
