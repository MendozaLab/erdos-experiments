/-
Bridge: Erdős 123, snug variant for 2, 3, 5 (Formal Conjectures statement), derived from the
Principia Math band theorem `Erdos123.Statement.erdos123_dcomplete_real`.

The mathematics is Principia Math's (Erdős #123 band theorem; see `Challenge.lean` and
`Solution.lean` in the Principia package, a pinned Lake dependency). This file only
specialises it to (a, b, c) = (2, 3, 5) and rewrites the conclusion in the exact vocabulary of
`FormalConjectures/ErdosProblems/123.lean`, theorem `erdos_123.variants.powers_2_3_5_snug`.
The definition `IsSnug` and the statement below are verbatim copies of that file.
-/

import Erdos123
import Solution

set_option autoImplicit false

open Filter
open Submonoid
open scoped Pointwise

namespace SnugBridge

/-- Verbatim copy of `Erdos123.IsSnug` from google-deepmind/formal-conjectures. -/
def IsSnug (ε : ℝ) (A : Finset ℕ) : Prop :=
  ∃ hA : A.Nonempty, ∀ a ∈ A, a < (1 + ε) * A.min' hA

/-- Every `2^k * 3^l * 5^m` lies in the pointwise product of the three power submonoids. -/
lemma smooth_mem_powers {s : ℕ} (hs : s ∈ Erdos123.Statement.Smooth3 2 3 5) :
    s ∈ (↑(powers 2) * ↑(powers 3) * ↑(powers 5) : Set ℕ) := by
  obtain ⟨k, l, m, rfl⟩ := hs
  exact Set.mul_mem_mul (Set.mul_mem_mul ⟨k, rfl⟩ ⟨l, rfl⟩) ⟨m, rfl⟩

/-- **Erdős 123, snug variant for 2, 3, 5**, exactly as stated in Formal Conjectures. -/
theorem erdos_123_powers_2_3_5_snug :
    ∀ ε > 0, ∀ᶠ n in atTop,
      ∃ A : Finset ℕ, (A : Set ℕ) ⊆ ↑(powers 2) * ↑(powers 3) * ↑(powers 5) ∧ IsSnug ε A ∧
        ∑ x ∈ A, x = n := by
  intro ε hε
  -- Band ratio: ρ = min (1 + ε) (3/2), so 1 < ρ ≤ 1 + ε and ρ < 2 = min 2 (min 3 5).
  set ρ : ℝ := min (1 + ε) (3 / 2) with hρ
  have hρ1 : (1 : ℝ) < ρ := lt_min (by linarith) (by norm_num)
  have hρe : ρ ≤ 1 + ε := min_le_left _ _
  have h32 : ρ ≤ 3 / 2 := min_le_right _ _
  have hco : Erdos123.Statement.PairwiseCoprime3 2 3 5 := by
    refine ⟨?_, ?_, ?_⟩ <;> decide
  obtain ⟨N₀, hN₀⟩ := Erdos123.Statement.erdos123_dcomplete_real 2 3 5 ρ
    (by norm_num) (by norm_num) (by norm_num) hco hρ1 (by norm_num; linarith)
  rw [Filter.eventually_atTop]
  refine ⟨max N₀ 1, fun n hn => ?_⟩
  obtain ⟨x, T, hT, -, hsum⟩ := hN₀ n (le_of_max_le_left hn)
  have hn1 : 1 ≤ n := le_of_max_le_right hn
  -- T is nonempty because its sum is positive.
  have hTne : T.Nonempty := by
    rcases T.eq_empty_or_nonempty with h | h
    · subst h; simp at hsum; omega
    · exact h
  refine ⟨T, fun s hs => smooth_mem_powers (hT s hs).1, ⟨hTne, fun a ha => ?_⟩, by simpa using hsum⟩
  -- a < ρ x ≤ ρ (min T) ≤ (1 + ε) (min T).
  have hmin := hT (T.min' hTne) (T.min'_mem hTne)
  have ha' := hT a ha
  have hρpos : (0 : ℝ) ≤ ρ := le_of_lt (lt_trans zero_lt_one hρ1)
  have hminpos : (0 : ℝ) ≤ (T.min' hTne : ℝ) := Nat.cast_nonneg _
  calc (a : ℝ) < ρ * x := ha'.2.2
    _ ≤ ρ * (T.min' hTne : ℝ) := mul_le_mul_of_nonneg_left hmin.2.1 hρpos
    _ ≤ (1 + ε) * (T.min' hTne : ℝ) := mul_le_mul_of_nonneg_right hρe hminpos

end SnugBridge

#print axioms SnugBridge.erdos_123_powers_2_3_5_snug
