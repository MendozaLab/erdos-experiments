import Erdos604Proof.Definitions

open Filter
open scoped Topology
open OAI.WeakPinned

namespace PinnedDistanceAdapter

lemma pairFraction_le_two (P : Finset Plane) (s : ℝ) (hP : 2 ≤ P.card) :
    pairFraction P s ≤ 2 := by
  classical
  have hn : (2 : ℝ) ≤ P.card := by exact_mod_cast hP
  have hden : 0 < (P.card : ℝ) * ((P.card : ℝ) - 1) := mul_pos (by linarith) (by linarith)
  have hcard : (richPairs P s).card ≤ P.card * P.card := by
    simpa [richPairs] using Finset.card_le_card (Finset.filter_subset
      (fun xy : Plane × Plane => xy.1 ≠ xy.2 ∧
        (P.card : ℝ) ^ s ≤ (k P xy.1 xy.2 : ℝ)) (P ×ˢ P))
  have hreal : ((richPairs P s).card : ℝ) ≤ (P.card : ℝ) * P.card := by
    exact_mod_cast hcard
  unfold pairFraction
  apply (div_le_iff₀ hden).2
  nlinarith

lemma pairFraction_le_F (P : Finset Plane) (s : ℝ) (hP : 2 ≤ P.card) :
    pairFraction P s ≤ F P.card s := by
  unfold F
  rw [if_pos hP]
  apply le_csSup
  · refine ⟨2, ?_⟩
    rintro a ⟨Q, hQ, rfl⟩
    exact pairFraction_le_two Q s (hQ ▸ hP)
  · exact ⟨P, rfl, rfl⟩

/-- The endpoint's supremum supplies a single cardinality cutoff for every configuration. -/
theorem uniform_pairFraction
    (endpoint : ∀ s : ℝ, 0 < s → Tendsto (fun n : ℕ => F n s) atTop (𝓝 0))
    (s δ : ℝ) (hs : 0 < s) (hδ : 0 < δ) :
    ∃ N : ℕ, 2 ≤ N ∧ ∀ P : Finset Plane, N ≤ P.card → pairFraction P s < δ := by
  have h := (endpoint s hs).eventually (gt_mem_nhds hδ)
  obtain ⟨N, hN⟩ := eventually_atTop.1 h
  refine ⟨max N 2, le_max_right _ _, ?_⟩
  intro P hP
  exact (pairFraction_le_F P s (le_trans (le_max_right _ _) hP)).trans_lt
    (hN P.card (le_trans (le_max_left _ _) hP))


end PinnedDistanceAdapter
