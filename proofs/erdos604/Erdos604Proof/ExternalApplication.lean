import OAI.Geometry.PinnedDistances.Main
import Erdos604Proof.Target

open Filter
open scoped Topology
open OAI.WeakPinned

/-- Actual endpoint application, uniform in every sufficiently large configuration. -/
theorem almost_all_pins_from_openai (ε δ : ℝ) (hε : 0 < ε) (hδ : 0 < δ) :
    ∃ N : ℕ, 2 ≤ N ∧ ∀ P : Finset Plane, N ≤ P.card →
      (1 - δ) * (P.card : ℝ) <
        ((P.filter (fun x => (P.card : ℝ) ^ (1 - ε) <
          ((PinnedDistanceAdapter.pinnedDistances P x).card : ℝ))).card : ℝ) :=
  PinnedDistanceAdapter.almost_all_pins_good OAI.WeakPinned.main ε δ hε hδ

/-- The weak Erdős 604 statement, using the imported OpenAI endpoint. -/
theorem erdos604_from_openai :
    ∀ ε : ℝ, 0 < ε → ∀ᶠ n : ℕ in atTop,
      (n : ℝ) ^ (1 - ε) ≤ (Erdos604Proof.Target.pinnedDistinctDistances n : ℝ) :=
  PinnedDistanceAdapter.erdos604_of_main OAI.WeakPinned.main
