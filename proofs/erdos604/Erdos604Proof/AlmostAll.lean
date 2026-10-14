import Erdos604Proof.Counting
import Erdos604Proof.Uniform

/-! Conditional almost-all-pins consequence of the endpoint's uniform rich-pair limit.
This file uses no theorem from the external analytic/arithmetic proof closure. -/

open Filter
open scoped Topology
open OAI.WeakPinned

namespace PinnedDistanceAdapter

/-- The small-fiber contribution is eventually at most half of the nonself row. -/
lemma eventually_small_fiber_loss (ε : ℝ) (hε : 0 < ε) :
    ∃ N : ℕ, 4 ≤ N ∧ ∀ n : ℕ, N ≤ n →
      (n : ℝ) ^ (1 - ε) * (n : ℝ) ^ (ε / 2) ≤ ((n : ℝ) - 1) / 2 := by
  have hlim : Tendsto (fun n : ℕ => (n : ℝ) ^ (-(ε / 2))) atTop (𝓝 0) :=
    (tendsto_rpow_neg_atTop (by linarith : 0 < ε / 2)).comp
      tendsto_natCast_atTop_atTop
  obtain ⟨N, hN⟩ := eventually_atTop.1
    (hlim.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 4)))
  refine ⟨max N 4, le_max_right _ _, ?_⟩
  intro n hn
  have hn4 : (4 : ℝ) ≤ n := by exact_mod_cast (le_trans (le_max_right N 4) hn)
  have hn0 : (0 : ℝ) < n := by linarith
  have hsmall := hN n (le_trans (le_max_left N 4) hn)
  have hpow : (n : ℝ) ^ (1 - ε) * (n : ℝ) ^ (ε / 2) =
      (n : ℝ) * (n : ℝ) ^ (-(ε / 2)) := by
    calc
      _ = (n : ℝ) ^ ((1 - ε) + ε / 2) := (Real.rpow_add hn0 _ _).symm
      _ = (n : ℝ) ^ (1 + -(ε / 2)) := by congr 1; ring
      _ = (n : ℝ) * (n : ℝ) ^ (-(ε / 2)) := by rw [Real.rpow_add hn0, Real.rpow_one]
  rw [hpow]
  nlinarith

/-- For every positive exponent loss and exceptional proportion, one cardinality cutoff
works for every planar configuration. Distances to the pin itself are excluded. -/
theorem almost_all_pins
    (endpoint : ∀ s : ℝ, 0 < s → Tendsto (fun n : ℕ => F n s) atTop (𝓝 0))
    (ε δ : ℝ) (hε : 0 < ε) (hδ : 0 < δ) :
    ∃ N : ℕ, 2 ≤ N ∧ ∀ P : Finset Plane, N ≤ P.card →
      ((badPins P ((P.card : ℝ) ^ (1 - ε))).card : ℝ) < δ * (P.card : ℝ) := by
  obtain ⟨N₁, hN₁, hsmall⟩ := eventually_small_fiber_loss ε hε
  obtain ⟨N₂, hN₂, hfrac⟩ := uniform_pairFraction endpoint (ε / 2) (δ / 2)
    (by linarith) (by linarith)
  refine ⟨max N₁ N₂, le_trans hN₂ (le_max_right _ _), ?_⟩
  intro P hP
  have hP₁ : N₁ ≤ P.card := le_trans (le_max_left _ _) hP
  have hP₂ : N₂ ≤ P.card := le_trans (le_max_right _ _) hP
  have hPtwo : 2 ≤ P.card := le_trans hN₂ hP₂
  have hPpos : (0 : ℝ) < P.card := by exact_mod_cast (by omega : 0 < P.card)
  apply (div_lt_iff₀ hPpos).1
  calc
    ((badPins P ((P.card : ℝ) ^ (1 - ε))).card : ℝ) / (P.card : ℝ)
        ≤ 2 * pairFraction P (ε / 2) :=
      badPins_fraction_le_two_pairFraction P _ _ hPtwo (hsmall P.card hP₁)
    _ < δ := by linarith [hfrac P hP₂]

/-- Equivalently, more than a `1 - δ` proportion of pins see strictly more than
`n ^ (1 - ε)` distinct nonzero distances, uniformly over all large configurations. -/
theorem almost_all_pins_good
    (endpoint : ∀ s : ℝ, 0 < s → Tendsto (fun n : ℕ => F n s) atTop (𝓝 0))
    (ε δ : ℝ) (hε : 0 < ε) (hδ : 0 < δ) :
    ∃ N : ℕ, 2 ≤ N ∧ ∀ P : Finset Plane, N ≤ P.card →
      (1 - δ) * (P.card : ℝ) <
        ((P.filter (fun x => (P.card : ℝ) ^ (1 - ε) <
          ((pinnedDistances P x).card : ℝ))).card : ℝ) := by
  classical
  obtain ⟨N, hN, hbad⟩ := almost_all_pins endpoint ε δ hε hδ
  refine ⟨N, hN, ?_⟩
  intro P hP
  have hpartition : (badPins P ((P.card : ℝ) ^ (1 - ε))).card +
      (P.filter (fun x => (P.card : ℝ) ^ (1 - ε) <
        ((pinnedDistances P x).card : ℝ))).card = P.card := by
    simpa only [badPins, not_le] using
      (Finset.card_filter_add_card_filter_not (s := P)
        (fun x => ((pinnedDistances P x).card : ℝ) ≤ (P.card : ℝ) ^ (1 - ε)))
  have hpartition_real : ((badPins P ((P.card : ℝ) ^ (1 - ε))).card : ℝ) +
      ((P.filter (fun x => (P.card : ℝ) ^ (1 - ε) <
        ((pinnedDistances P x).card : ℝ))).card : ℝ) = P.card := by
    exact_mod_cast hpartition
  linarith [hbad P hP]


end PinnedDistanceAdapter
