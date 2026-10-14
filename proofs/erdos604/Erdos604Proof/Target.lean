import Erdos604Proof.AlmostAll
import Erdos604Proof.TargetDefinitions

open Filter
open scoped Topology
open OAI.WeakPinned Erdos604Proof.Target

namespace PinnedDistanceAdapter

/-- The local distance count is exactly the current target count, including erasure of the pin. -/
theorem pinnedDistances_card_eq (P : Finset Plane) (x : Plane) :
    (pinnedDistances P x).card = distinctDistancesFrom P x := by
  classical
  unfold pinnedDistances distinctDistancesFrom
  congr 1
  ext r
  simp

/-- Every cardinality is realized, so the target's natural infimum is attained. -/
theorem pinnedDistinctDistances_attained (n : ℕ) :
    ∃ P : Finset Plane, P.card = n ∧
      P.sup (distinctDistancesFrom P) = Erdos604Proof.Target.pinnedDistinctDistances n := by
  have hne : Set.Nonempty {m : ℕ | ∃ P : Finset Plane,
      P.card = n ∧ P.sup (distinctDistancesFrom P) = m} := by
    obtain ⟨P, hP⟩ := Finset.exists_card_eq (α := Plane) n
    exact ⟨_, P, hP, rfl⟩
  exact Nat.sInf_mem hne

/-- A conditional application to the precise right-hand side of current Erdős 604. -/
theorem erdos604_of_main
    (endpoint : ∀ s : ℝ, 0 < s → Tendsto (fun n : ℕ => F n s) atTop (𝓝 0)) :
    ∀ ε : ℝ, 0 < ε → ∀ᶠ n : ℕ in atTop,
      (n : ℝ) ^ (1 - ε) ≤ (Erdos604Proof.Target.pinnedDistinctDistances n : ℝ) := by
  classical
  intro ε hε
  obtain ⟨N, hN2, hN⟩ := almost_all_pins endpoint ε 1 hε zero_lt_one
  apply eventually_atTop.2
  refine ⟨N, ?_⟩
  intro n hn
  obtain ⟨P, hPn, hmin⟩ := pinnedDistinctDistances_attained n
  have hbad := hN P (hPn ▸ hn)
  have hex : ∃ x ∈ P, (n : ℝ) ^ (1 - ε) < (distinctDistancesFrom P x : ℝ) := by
    by_contra! h
    have heq : badPins P ((P.card : ℝ) ^ (1 - ε)) = P := by
      apply Finset.filter_eq_self.2
      intro x hx
      simpa only [pinnedDistances_card_eq, hPn] using h x hx
    rw [heq, one_mul] at hbad
    exact lt_irrefl _ hbad
  obtain ⟨x, hx, hxd⟩ := hex
  have hsup : distinctDistancesFrom P x ≤ P.sup (distinctDistancesFrom P) :=
    Finset.le_sup hx
  have hsup' : (distinctDistancesFrom P x : ℝ) ≤
      (Erdos604Proof.Target.pinnedDistinctDistances n : ℝ) := by exact_mod_cast (hmin ▸ hsup)
  exact hxd.le.trans hsup'

end PinnedDistanceAdapter
