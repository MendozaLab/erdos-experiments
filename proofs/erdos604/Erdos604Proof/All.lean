import Erdos604Proof.Target

open OAI.WeakPinned PinnedDistanceAdapter

/-! Aggregate check and degenerate-cardinality sanity checks for the local adapter. -/

example (x : Plane) : pinnedDistances ∅ x = ∅ := by
  classical
  simp [pinnedDistances]

example (x : Plane) : pinnedDistances {x} x = ∅ := by
  classical
  simp [pinnedDistances]

example (s : ℝ) : richPairs ∅ s = ∅ := by
  classical
  simp [richPairs]

example (x : Plane) (s : ℝ) : richPairs {x} s = ∅ := by
  classical
  simp [richPairs]

example (P : Finset Plane) (x : Plane) (r : ℝ) (hr : r ∈ pinnedDistances P x) :
    0 < r := by
  obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hr
  exact dist_pos.2 (Finset.mem_erase.1 hy).1
