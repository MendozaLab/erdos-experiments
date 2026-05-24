/-
  Erdős #30 / BFR interval partition scout
  ========================================

  Review-only scout file. This extends the interval-occupancy work by showing
  that the clipped BFR intervals cover the full support range and are pairwise
  disjoint when the width is positive.

  Claim ceiling: local formalization scout only. It does not close the BFR core
  bound, improve the public Sidon coefficient, or alter registry/publication
  state.
-/

import Erdos30_BFR

open Finset Nat

namespace Erdos.Sidon

/-! ## Finite interval partition of the BFR support range -/

/-- The `i`-th width-`v` interval, clipped to the natural BFR support range. -/
def bfrPartitionInterval (N v i : ℕ) : Finset ℕ :=
  (bfrSumRange N).filter (fun s : ℕ => i * v ≤ s ∧ s < (i + 1) * v)

/-- A finite index set large enough to cover `bfrSumRange N` for every positive width. -/
def bfrPartitionIndexSet (N _v : ℕ) : Finset ℕ :=
  Finset.range (2 * N + 1)

/-- Occupancy of one partition interval by ordered Sidon representations. -/
def bfrPartitionOccupancy (A : Finset ℕ) (N v i : ℕ) : ℕ :=
  ∑ s ∈ bfrPartitionInterval N v i, bfrRepFunction A s

theorem mem_bfrPartitionInterval_iff (N v i s : ℕ) :
    s ∈ bfrPartitionInterval N v i ↔
      s ≤ 2 * N ∧ i * v ≤ s ∧ s < (i + 1) * v := by
  unfold bfrPartitionInterval
  simp [mem_bfrSumRange_iff]

theorem bfrPartitionInterval_subset_sumRange (N v i : ℕ) :
    bfrPartitionInterval N v i ⊆ bfrSumRange N := by
  unfold bfrPartitionInterval
  exact Finset.filter_subset _ _

/-- Every represented sum in the BFR range lies in its quotient-index interval. -/
theorem bfrPartitionInterval_mem_div_index
    (N v s : ℕ) (hv : 0 < v) (hs : s ∈ bfrSumRange N) :
    s ∈ bfrPartitionInterval N v (s / v) := by
  rw [mem_bfrPartitionInterval_iff]
  have hsN : s ≤ 2 * N := (mem_bfrSumRange_iff N s).mp hs
  have hlow : s / v * v ≤ s := Nat.div_mul_le_self s v
  have hupper_mul : s < v * (s / v + 1) := Nat.lt_mul_div_succ s hv
  have hupper : s < (s / v + 1) * v := by
    simpa [Nat.mul_comm] using hupper_mul
  exact ⟨hsN, hlow, hupper⟩

/-- The finite interval family covers the BFR support range. -/
theorem bfrPartitionInterval_cover_sumRange
    (N v : ℕ) (hv : 0 < v) :
    (bfrPartitionIndexSet N v).biUnion (bfrPartitionInterval N v) = bfrSumRange N := by
  ext s
  constructor
  · intro hs
    rw [Finset.mem_biUnion] at hs
    rcases hs with ⟨i, hi, hsi⟩
    exact bfrPartitionInterval_subset_sumRange N v i hsi
  · intro hs
    rw [Finset.mem_biUnion]
    refine ⟨s / v, ?_, bfrPartitionInterval_mem_div_index N v s hv hs⟩
    unfold bfrPartitionIndexSet
    rw [Finset.mem_range]
    have hsN : s ≤ 2 * N := (mem_bfrSumRange_iff N s).mp hs
    have hdiv : s / v ≤ s := Nat.div_le_self s v
    omega

/-- Distinct positive-width partition intervals are disjoint. -/
theorem bfrPartitionInterval_disjoint_of_ne
    (N v i j : ℕ) (hij : i ≠ j) :
    Disjoint (bfrPartitionInterval N v i) (bfrPartitionInterval N v j) := by
  rw [Finset.disjoint_left]
  intro s hsi hsj
  rw [mem_bfrPartitionInterval_iff] at hsi hsj
  rcases hsi with ⟨_, hi_low, hi_high⟩
  rcases hsj with ⟨_, hj_low, hj_high⟩
  rcases lt_or_gt_of_ne hij with hij_lt | hji_lt
  · have hij_succ : i + 1 ≤ j := Nat.succ_le_of_lt hij_lt
    have hgap : (i + 1) * v ≤ j * v := Nat.mul_le_mul_right v hij_succ
    have hs_lt_j : s < j * v := lt_of_lt_of_le hi_high hgap
    exact (not_lt.mpr hj_low) hs_lt_j
  · have hji_succ : j + 1 ≤ i := Nat.succ_le_of_lt hji_lt
    have hgap : (j + 1) * v ≤ i * v := Nat.mul_le_mul_right v hji_succ
    have hs_lt_i : s < i * v := lt_of_lt_of_le hj_high hgap
    exact (not_lt.mpr hi_low) hs_lt_i

/-- The partition intervals are pairwise disjoint over the finite index set. -/
theorem bfrPartitionInterval_pairwiseDisjoint (N v : ℕ) :
    Set.PairwiseDisjoint (↑(bfrPartitionIndexSet N v) : Set ℕ)
      (bfrPartitionInterval N v) := by
  intro i hi j hj hij
  exact bfrPartitionInterval_disjoint_of_ne N v i j hij

/-- Summing interval occupancies over the partition recovers the support-range mass. -/
theorem bfrPartitionOccupancy_sum_eq_total
    (A : Finset ℕ) (N v : ℕ) (hv : 0 < v) :
    ∑ i ∈ bfrPartitionIndexSet N v, bfrPartitionOccupancy A N v i =
      ∑ s ∈ bfrSumRange N, bfrRepFunction A s := by
  unfold bfrPartitionOccupancy
  rw [← Finset.sum_biUnion (s := bfrPartitionIndexSet N v)
    (t := bfrPartitionInterval N v)
    (f := fun s => bfrRepFunction A s)
    (bfrPartitionInterval_pairwiseDisjoint N v)]
  rw [bfrPartitionInterval_cover_sumRange N v hv]

/-- With `A ⊆ {0,...,N}`, the partitioned interval mass is `|A|^2`. -/
theorem bfrPartitionOccupancy_sum_eq_card_sq
    (A : Finset ℕ) (N v : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) (hv : 0 < v) :
    ∑ i ∈ bfrPartitionIndexSet N v, bfrPartitionOccupancy A N v i =
      A.card * A.card := by
  rw [bfrPartitionOccupancy_sum_eq_total A N v hv]
  exact bfr_sum_repFunction_eq_card_sq A N hA

end Erdos.Sidon
