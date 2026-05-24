/-
  Erdős Problem #755 — Symmetry Quotient for B_2[g]

  The existing #755 files use the ordered convention: every sum has at most
  2g ordered pair representations.  The native combinatorial statement is
  usually unordered: count only pairs (a,b) with a ≤ b, and allow at most g.

  This file closes that translation gap for h = 2.  The "shadow" term is the
  reversed strict half of the ordered fiber.  It injects into the unordered
  upper-triangle fiber by swapping coordinates, so the ordered fiber has at
  most twice the unordered cardinality.

  Lean version: leanprover/lean4:v4.27.0
  Mathlib version: v4.27.0
-/

import Mathlib
import Erdos755_BhG

open Finset Nat

namespace Erdos.B2G.SymmetryQuotient

open Erdos.B2G

/-- Ordered sum fiber in `A × A`. -/
def orderedFiber (A : Finset ℕ) (s : ℕ) : Finset (ℕ × ℕ) :=
  (A ×ˢ A).filter (fun p : ℕ × ℕ => p.1 + p.2 = s)

/-- Native unordered sum fiber: keep only the upper triangle `a ≤ b`. -/
def unorderedFiber (A : Finset ℕ) (s : ℕ) : Finset (ℕ × ℕ) :=
  (A ×ˢ A).filter (fun p : ℕ × ℕ => p.1 ≤ p.2 ∧ p.1 + p.2 = s)

/-- The shadow half of the ordered fiber: strict lower triangle `b < a`. -/
def reversedStrictFiber (A : Finset ℕ) (s : ℕ) : Finset (ℕ × ℕ) :=
  (A ×ˢ A).filter (fun p : ℕ × ℕ => p.2 < p.1 ∧ p.1 + p.2 = s)

/-- Native unordered B_2[g] predicate. -/
abbrev IsB2GSetUnordered (A : Finset ℕ) (g : ℕ) : Prop :=
  ∀ s : ℕ, (unorderedFiber A s).card ≤ g

/-- The strict reversed half injects into the unordered half by swapping. -/
theorem reversedStrictFiber_card_le_unorderedFiber (A : Finset ℕ) (s : ℕ) :
    (reversedStrictFiber A s).card ≤ (unorderedFiber A s).card := by
  classical
  refine Finset.card_le_card_of_injOn (fun p : ℕ × ℕ => (p.2, p.1)) ?maps ?inj
  · intro p hp
    simp [reversedStrictFiber, unorderedFiber] at hp ⊢
    exact ⟨⟨hp.1.2, hp.1.1⟩, Nat.le_of_lt hp.2.1, by omega⟩
  · intro p _hp q _hq hpq
    cases p with
    | mk p₁ p₂ =>
      cases q with
      | mk q₁ q₂ =>
        simp at hpq ⊢
        exact ⟨hpq.2, hpq.1⟩

/-- Ordered fibers are covered by the unordered half plus the reversed shadow. -/
theorem orderedFiber_card_le_unordered_plus_shadow (A : Finset ℕ) (s : ℕ) :
    (orderedFiber A s).card ≤
      (unorderedFiber A s).card + (reversedStrictFiber A s).card := by
  classical
  have h_subset :
      orderedFiber A s ⊆ unorderedFiber A s ∪ reversedStrictFiber A s := by
    intro p hp
    simp [orderedFiber, unorderedFiber, reversedStrictFiber] at hp ⊢
    rcases hp with ⟨⟨hpA₁, hpA₂⟩, hsum⟩
    by_cases hle : p.1 ≤ p.2
    · exact Or.inl ⟨⟨hpA₁, hpA₂⟩, hle, hsum⟩
    · have hlt : p.2 < p.1 := by omega
      exact Or.inr ⟨⟨hpA₁, hpA₂⟩, hlt, hsum⟩
  exact (Finset.card_le_card h_subset).trans (Finset.card_union_le _ _)

/-- Native unordered B_2[g] implies the ordered convention used by `Erdos755_BhG`. -/
theorem unordered_to_ordered (A : Finset ℕ) (g : ℕ)
    (hU : IsB2GSetUnordered A g) :
    IsB2GSet A g := by
  intro s
  have h_cover := orderedFiber_card_le_unordered_plus_shadow A s
  have h_shadow := reversedStrictFiber_card_le_unorderedFiber A s
  have h_upper := hU s
  change (orderedFiber A s).card ≤ 2 * g
  omega

/-- Sum-counting bound stated against the native unordered B_2[g] predicate. -/
theorem b2g_sum_count_unordered (A : Finset ℕ) (N g : ℕ)
    (hU : IsB2GSetUnordered A g)
    (hA : A ⊆ Finset.range (N + 1)) :
    A.card ^ 2 ≤ 2 * g * (2 * N + 1) := by
  exact b2g_sum_count A N g (unordered_to_ordered A g hU) hA

/-- Square-form corollary under the native unordered predicate. -/
theorem b2g_card_sq_bound_unordered (A : Finset ℕ) (N g : ℕ)
    (hU : IsB2GSetUnordered A g)
    (hA : A ⊆ Finset.range (N + 1)) :
    A.card * A.card ≤ 2 * g * (2 * N + 1) := by
  have h := b2g_sum_count_unordered A N g hU hA
  have : A.card ^ 2 = A.card * A.card := sq A.card
  omega

end Erdos.B2G.SymmetryQuotient
