/-
  Erdős Problem #30 — Difference-Counting Bound for Sidon Sets

  Theorem: For a Sidon set A ⊆ [N], |A| * (|A| - 1) ≤ 2N

  This is the standard textbook bound, NOT the Erdős-Turán conjecture.
  The prize ($1,000) requires proving h(N) = √N + o(√N).
  Best known: Lindström (1969), h(N) ≤ √N + N^(1/4) + 1

  Status: TIER 3 — textbook bound only, not prize-eligible

  Proof sketch (difference counting):
    A Sidon set A has all pairwise sums a+b (a ≤ b) distinct.
    Equivalently, all pairwise differences a-b (a ≠ b) are distinct.
    For ordered pairs (a,b) with a > b, the difference a-b ∈ {1,...,N}.
    There are k(k-1)/2 such pairs (where k = |A|).
    By injectivity, k(k-1)/2 ≤ N, so k(k-1) ≤ 2N.

  Note: The equivalence between distinct sums and distinct differences
  follows from: a₁+b₁ = a₂+b₂ iff a₁-a₂ = b₂-b₁.

  Lean version: leanprover/lean4:v4.24.0
  Mathlib version: f897ebcf72cd16f89ab4577d0c826cd14afaafc7
-/

import Mathlib
import Erdos30_Sidon_Defs

open Finset Nat Erdos.Sidon

namespace Erdos

/-- Key lemma: In a Sidon set, if a₁ > b₁ and a₂ > b₂ and a₁ - b₁ = a₂ - b₂,
    then a₁ = a₂ and b₁ = b₂. This follows from the Sidon property because
    a₁ - b₁ = a₂ - b₂ implies a₁ + b₂ = a₂ + b₁. -/
theorem sidon_diff_injective (A : Finset ℕ)
    (hS : IsSidonSet A)
    (a₁ b₁ a₂ b₂ : ℕ)
    (ha₁ : a₁ ∈ A) (hb₁ : b₁ ∈ A) (ha₂ : a₂ ∈ A) (hb₂ : b₂ ∈ A)
    (hlt₁ : b₁ < a₁) (hlt₂ : b₂ < a₂)
    (heq : a₁ - b₁ = a₂ - b₂) :
    a₁ = a₂ ∧ b₁ = b₂ := by
  -- From a₁ - b₁ = a₂ - b₂ (in ℕ, with a₁ > b₁ and a₂ > b₂),
  -- we get a₁ + b₂ = a₂ + b₁
  have h_sum : a₁ + b₂ = a₂ + b₁ := by omega
  -- Apply Sidon property: need to arrange so first arg ≤ second arg
  -- We know a₁ + b₂ = a₂ + b₁, and we need ordered pairs
  -- Case split on whether b₂ ≤ a₁ and b₁ ≤ a₂ (or the reverse)
  by_cases h1 : b₂ ≤ a₁
  · by_cases h2 : b₁ ≤ a₂
    · -- b₂ ≤ a₁ and b₁ ≤ a₂: apply Sidon to (b₂, a₁) and (b₁, a₂)
      have := hS b₂ hb₂ a₁ ha₁ b₁ hb₁ a₂ ha₂ h1 h2 (by omega)
      exact ⟨this.2, this.1⟩
    · -- b₂ ≤ a₁ and a₂ < b₁: apply Sidon to (b₂, a₁) and (a₂, b₁)
      push_neg at h2
      have := hS b₂ hb₂ a₁ ha₁ a₂ ha₂ b₁ hb₁ h1 (Nat.le_of_lt h2) (by omega)
      exact ⟨this.2, this.1.symm ▸ (by omega)⟩
  · push_neg at h1
    by_cases h2 : b₁ ≤ a₂
    · -- a₁ < b₂ and b₁ ≤ a₂: apply Sidon to (a₁, b₂) and (b₁, a₂)
      have := hS a₁ ha₁ b₂ hb₂ b₁ hb₁ a₂ ha₂ (Nat.le_of_lt h1) h2 (by omega)
      exact ⟨this.1.symm ▸ (by omega), this.2⟩
    · -- a₁ < b₂ and a₂ < b₁: apply Sidon to (a₁, b₂) and (a₂, b₁)
      push_neg at h2
      have := hS a₁ ha₁ b₂ hb₂ a₂ ha₂ b₁ hb₁ (Nat.le_of_lt h1) (Nat.le_of_lt h2) (by omega)
      constructor <;> omega

/--
  **Difference-Counting Bound for Sidon Sets**

  For a Sidon set A ⊆ Finset.range (N + 1), we have:
    A.card * (A.card - 1) ≤ 2 * N

  Proof: Consider ordered pairs (a, b) with a, b ∈ A and b < a.
  There are k(k-1)/2 such pairs (where k = |A|).
  The map (a,b) ↦ a - b is injective by the Sidon property.
  Each difference a - b lies in {1, ..., N} (since 0 ≤ b < a ≤ N).
  So k(k-1)/2 ≤ N, hence k(k-1) ≤ 2N.

  This gives h(N) ≤ √(2N) + O(1), the elementary textbook bound.
  The Erdős-Turán conjecture (prize: $1,000) asks for h(N) = √N + o(√N).
-/
theorem sidon_difference_count (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) :
    A.card * (A.card - 1) ≤ 2 * N := by
  -- Define the set of ordered pairs (a, b) with a, b ∈ A and b < a
  set pairs := Finset.filter (fun p : ℕ × ℕ => p.2 < p.1) (A ×ˢ A) with pairs_def
  -- The difference map (a,b) ↦ a - b
  set diff_map := fun (p : ℕ × ℕ) => p.1 - p.2
  -- Step 1: The image of diff_map on pairs is injective
  have h_inj : Set.InjOn diff_map (↑pairs) := by
    intro ⟨a₁, b₁⟩ h₁ ⟨a₂, b₂⟩ h₂ heq
    simp [pairs_def, Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at h₁ h₂
    have := sidon_diff_injective A hS a₁ b₁ a₂ b₂ h₁.1 h₁.2 h₂.1 h₂.2 h₁.3 h₂.3 heq
    exact Prod.ext this.1 this.2
  -- Step 2: Count the pairs — there are k*(k-1)/2 pairs with b < a
  -- Actually we count ALL ordered pairs with a ≠ b, which gives k*(k-1).
  -- But we use only b < a, which gives k*(k-1)/2.
  -- We'll use the full set of pairs with b < a.
  -- Step 3: The image lies in Finset.range N (differences are in {1,...,N})
  have h_range : Finset.image diff_map pairs ⊆ Finset.Icc 1 N := by
    intro d hd
    simp [Finset.mem_image, pairs_def, Finset.mem_filter, Finset.mem_product] at hd
    obtain ⟨⟨a, b⟩, ⟨ha, hb, hlt⟩, rfl⟩ := hd
    simp [Finset.mem_Icc]
    constructor
    · omega
    · have ha_le : a ≤ N := by
        have := Finset.mem_range.mp (hA ha)
        omega
      omega
  -- Step 4: Card of image = card of pairs (by injectivity)
  have h_card_img : (Finset.image diff_map pairs).card = pairs.card := by
    exact Finset.card_image_of_injOn h_inj
  -- Step 5: Card of image ≤ N (since image ⊆ {1,...,N})
  have h_card_le_N : (Finset.image diff_map pairs).card ≤ N := by
    calc (Finset.image diff_map pairs).card
        ≤ (Finset.Icc 1 N).card := Finset.card_le_card h_range
      _ = N := by simp
  -- Step 6: pairs.card = k*(k-1)/2
  -- We use the fact that |{(a,b) ∈ A×A : b < a}| = C(k, 2) = k*(k-1)/2
  have h_pairs_card : pairs.card = A.card * (A.card - 1) / 2 := by
    -- The pairs with b < a biject with 2-element subsets of A
    rw [show pairs = Finset.filter (fun p : ℕ × ℕ => p.2 < p.1) (A ×ˢ A) from rfl]
    have : (Finset.filter (fun p : ℕ × ℕ => p.2 < p.1) (A ×ˢ A)).card
        = (Finset.powersetCard 2 A).card := by
      refine Finset.card_bij (fun p _ => {p.1, p.2}) ?_ ?_ ?_
      · intro ⟨a, b⟩ hp
        simp [Finset.mem_filter, Finset.mem_product] at hp
        rw [Finset.mem_powersetCard]
        refine ⟨?_, ?_⟩
        · intro x hx
          simp [Finset.mem_insert, Finset.mem_singleton] at hx
          cases hx with
          | inl h => exact h ▸ hp.1
          | inr h => exact h ▸ hp.2.1
        · rw [Finset.card_pair]
          exact Nat.ne_of_gt hp.2.2
      · intro ⟨a₁, b₁⟩ h₁ ⟨a₂, b₂⟩ h₂ heq
        simp [Finset.mem_filter, Finset.mem_product] at h₁ h₂
        have h₁lt := h₁.2.2
        have h₂lt := h₂.2.2
        -- From pair equality {a₁, b₁} = {a₂, b₂} with strict orderings, deduce component equality
        suffices a₁ = a₂ ∧ b₁ = b₂ from Prod.ext this.1 this.2
        have ha : a₁ ∈ ({a₂, b₂} : Finset ℕ) := heq ▸ mem_insert_self a₁ _
        simp only [mem_insert, mem_singleton] at ha
        rcases ha with rfl | rfl
        · -- a₁ = a₂: deduce b₁ = b₂
          refine ⟨rfl, ?_⟩
          have hb : b₁ ∈ ({a₂, b₂} : Finset ℕ) :=
            heq ▸ mem_insert.mpr (Or.inr (mem_singleton_self b₁))
          simp only [mem_insert, mem_singleton] at hb
          rcases hb with rfl | rfl
          · omega
          · rfl
        · -- a₁ = b₂: contradiction from orderings
          exfalso
          have ha₂ : a₂ ∈ ({b₂, b₁} : Finset ℕ) := heq.symm ▸ mem_insert_self a₂ _
          simp only [mem_insert, mem_singleton] at ha₂
          rcases ha₂ with rfl | rfl <;> omega
      · intro s hs
        rw [Finset.mem_powersetCard] at hs
        obtain ⟨hsub, hcard⟩ := hs
        rw [Finset.card_eq_two] at hcard
        obtain ⟨x, y, hxy, rfl⟩ := hcard
        cases Nat.lt_or_gt_of_ne hxy with
        | inl hlt =>
          exact ⟨(y, x), by simp [Finset.mem_filter, Finset.mem_product, hsub, hlt,
            Finset.pair_comm], by simp [Finset.pair_comm]⟩
        | inr hgt =>
          exact ⟨(x, y), by simp [Finset.mem_filter, Finset.mem_product, hsub, hgt], rfl⟩
    rw [this, Finset.card_powersetCard, Nat.choose_two_right]
  -- Step 7: Combine: k*(k-1)/2 ≤ N, so k*(k-1) ≤ 2*N
  have h_half_le : A.card * (A.card - 1) / 2 ≤ N := by
    rw [← h_pairs_card, ← h_card_img]
    exact h_card_le_N
  omega

/--
  Corollary: For a Sidon set A ⊆ {1,...,N}, |A|² ≤ 2N + |A|.
  This follows directly from |A|*(|A|-1) ≤ 2N by expanding.
  Since |A| ≤ 2N (trivially), this gives |A|² ≤ 4N, i.e., |A| ≤ 2√N.
-/
theorem sidon_card_sq_bound (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) :
    A.card ^ 2 ≤ 2 * N + A.card := by
  have h := sidon_difference_count A N hS hA
  have : A.card ^ 2 = A.card * (A.card - 1) + A.card := by
    cases A.card with
    | zero => simp
    | succ k => ring
  omega


end Erdos
