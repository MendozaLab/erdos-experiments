/-
  Erdős #30 / BFR interval occupancy scout
  ========================================

  Review-only scout file. This starts the next bridge after the representation
  function port: finite interval occupancy over the BFR sum range.

  Claim ceiling: local formalization scout only. It does not close the BFR core
  bound, improve the public Sidon coefficient, or alter registry/publication
  state.
-/

import Erdos30_BFR

open Finset Nat

namespace Erdos.Sidon

/-! ## Interval occupancy over the BFR sum range -/

/-- The `i`-th interval of width `v`, clipped to the natural BFR support range. -/
def bfrIntervalInRange (N v i : ℕ) : Finset ℕ :=
  (bfrSumRange N).filter (fun s : ℕ => i * v ≤ s ∧ s < (i + 1) * v)

theorem mem_bfrIntervalInRange_iff (N v i s : ℕ) :
    s ∈ bfrIntervalInRange N v i ↔
      s ≤ 2 * N ∧ i * v ≤ s ∧ s < (i + 1) * v := by
  unfold bfrIntervalInRange
  simp [mem_bfrSumRange_iff]

/-- Occupancy of a clipped BFR interval by ordered Sidon representations. -/
def bfrIntervalOccupancy (A : Finset ℕ) (N v i : ℕ) : ℕ :=
  ∑ s ∈ bfrIntervalInRange N v i, bfrRepFunction A s

theorem bfrIntervalInRange_subset_sumRange (N v i : ℕ) :
    bfrIntervalInRange N v i ⊆ bfrSumRange N := by
  unfold bfrIntervalInRange
  exact Finset.filter_subset _ _

/-- A single interval occupancy is bounded by the total ordered representation mass. -/
theorem bfrIntervalOccupancy_le_total (A : Finset ℕ) (N v i : ℕ) :
    bfrIntervalOccupancy A N v i ≤
      ∑ s ∈ bfrSumRange N, bfrRepFunction A s := by
  unfold bfrIntervalOccupancy
  exact Finset.sum_le_sum_of_subset_of_nonneg
    (bfrIntervalInRange_subset_sumRange N v i)
    (by
      intro s hs hs_not
      exact Nat.zero_le _)

/-- A single interval occupancy is bounded by `|A|^2`. -/
theorem bfrIntervalOccupancy_le_card_sq
    (A : Finset ℕ) (N v i : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    bfrIntervalOccupancy A N v i ≤ A.card * A.card := by
  calc
    bfrIntervalOccupancy A N v i
        ≤ ∑ s ∈ bfrSumRange N, bfrRepFunction A s :=
          bfrIntervalOccupancy_le_total A N v i
    _ = A.card * A.card := bfr_sum_repFunction_eq_card_sq A N hA

end Erdos.Sidon
