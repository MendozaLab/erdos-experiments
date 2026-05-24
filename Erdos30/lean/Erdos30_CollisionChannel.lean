/-
  Erdős Problem #30 — C2 Collision-Channel Framework

  This file extracts the finite counting spine behind the Sidon difference
  argument into reusable collision-channel lemmas. It deliberately avoids the
  unresolved Lindström, BFR, and Singer general-existence axioms.

  Status target: C2 automation/reusability, not progress on the open
  asymptotic Erdős #30 problem.
-/

import Mathlib
import Erdos30_Complete

open Finset

namespace Erdos.Sidon.CollisionChannel

open Erdos.Sidon

/-- A finite injective collision channel cannot have a larger source than
    codomain. This is the reusable pigeonhole spine behind the Sidon
    difference-counting proof. -/
theorem card_le_of_injective_channel
    {α β : Type*} [DecidableEq α] [DecidableEq β]
    (P : Finset α) (C : Finset β) (φ : α → β)
    (h_image : P.image φ ⊆ C)
    (h_inj : Set.InjOn φ ↑P) :
    P.card ≤ C.card := by
  calc
    P.card = (P.image φ).card := (Finset.card_image_of_injOn h_inj).symm
    _ ≤ C.card := Finset.card_le_card h_image

/-- Bounded-fiber version of the finite collision-channel count. This is the
    clean C2 extension point for B₂[g] / bounded-collision variants: if a
    source is counted as fibers over a codomain and every fiber has size at
    most `g`, then the source has size at most `g * C.card`. -/
theorem card_le_mul_of_fiber_card_bound
    {α β : Type*} [DecidableEq β]
    (P : Finset α) (C : Finset β) (fibers : β → Finset α) (g : ℕ)
    (h_count : P.card = C.sum (fun y => (fibers y).card))
    (h_fiber : ∀ y ∈ C, (fibers y).card ≤ g) :
    P.card ≤ g * C.card := by
  rw [h_count]
  calc
    C.sum (fun y => (fibers y).card) ≤ C.sum (fun _y => g) :=
      Finset.sum_le_sum h_fiber
    _ = C.card * g := by simp
    _ = g * C.card := by rw [Nat.mul_comm]

/-- The strict upper triangle of `A × A`: pairs `(a, b)` with `b < a`. -/
def strictUpper (A : Finset ℕ) : Finset (ℕ × ℕ) :=
  (A ×ˢ A).filter (fun p => p.2 < p.1)

/-- Differences on the strict upper triangle are injective for a Sidon set. -/
lemma diff_injOn_strictUpper (A : Finset ℕ) (hS : IsSidonSet A) :
    Set.InjOn (fun p : ℕ × ℕ => p.1 - p.2) ↑(strictUpper A) := by
  intro ⟨a₁, b₁⟩ h₁ ⟨a₂, b₂⟩ h₂ heq
  change (a₁, b₁) ∈ strictUpper A at h₁
  change (a₂, b₂) ∈ strictUpper A at h₂
  simp only [strictUpper, Finset.mem_filter, Finset.mem_product] at h₁ h₂
  obtain ⟨⟨ha₁, hb₁⟩, hlt₁⟩ := h₁
  obtain ⟨⟨ha₂, hb₂⟩, hlt₂⟩ := h₂
  have h := sidon_diff_injective A hS a₁ b₁ a₂ b₂ ha₁ hb₁ ha₂ hb₂ hlt₁ hlt₂ heq
  exact Prod.ext h.1 h.2

/-- Doubling identity for the strict upper triangle. Stated without division:
    `2 * |strictUpper A| = |A| * (|A| - 1)`. -/
lemma strictUpper_card_double (A : Finset ℕ) :
    2 * (strictUpper A).card = A.card * (A.card - 1) := by
  have h_off : A.offDiag.card = A.card * A.card - A.card := A.offDiag_card
  have h_swap : A.offDiag = strictUpper A ∪ (strictUpper A).image Prod.swap := by
    ext ⟨x, y⟩
    simp only [Finset.mem_offDiag, strictUpper, Finset.mem_union, Finset.mem_filter,
      Finset.mem_product, Finset.mem_image, Prod.swap, Prod.mk.injEq]
    constructor
    · rintro ⟨hx, hy, hne⟩
      rcases lt_or_gt_of_ne hne with h | h
      · right; exact ⟨(y, x), ⟨⟨hy, hx⟩, h⟩, rfl, rfl⟩
      · left; exact ⟨⟨hx, hy⟩, h⟩
    · rintro (⟨⟨hx, hy⟩, hlt⟩ | ⟨⟨a, b⟩, ⟨⟨ha, hb⟩, hlt⟩, hax, hby⟩)
      · exact ⟨hx, hy, ne_of_gt hlt⟩
      · subst hax; subst hby
        exact ⟨hb, ha, ne_of_lt hlt⟩
  have h_disj : Disjoint (strictUpper A) ((strictUpper A).image Prod.swap) := by
    rw [Finset.disjoint_left]
    intro ⟨x, y⟩ h1 h2
    simp only [strictUpper, Finset.mem_filter, Finset.mem_product] at h1
    simp only [Finset.mem_image, strictUpper, Finset.mem_filter, Finset.mem_product,
      Prod.swap, Prod.mk.injEq] at h2
    obtain ⟨_, hlt1⟩ := h1
    obtain ⟨⟨a, b⟩, ⟨_, hlt2⟩, hax, hby⟩ := h2
    subst hax; subst hby
    omega
  have h_img_card : ((strictUpper A).image Prod.swap).card = (strictUpper A).card := by
    apply Finset.card_image_of_injective
    intro ⟨a, b⟩ ⟨c, d⟩ h
    simp [Prod.swap, Prod.mk.injEq] at h
    ext <;> tauto
  have h_off_card :
      A.offDiag.card = (strictUpper A).card + ((strictUpper A).image Prod.swap).card := by
    rw [h_swap]; exact Finset.card_union_of_disjoint h_disj
  rw [h_img_card] at h_off_card
  rw [h_off] at h_off_card
  rw [Nat.mul_sub_one]
  omega

/-- The strict upper-triangle difference channel for a Sidon set in `[0, N]`
    lands in `Icc 1 N`. -/
lemma strictUpper_diff_image_subset_Icc (A : Finset ℕ) (N : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    (strictUpper A).image
        (fun p : ℕ × ℕ => p.1 - p.2) ⊆ Finset.Icc 1 N := by
  intro d hd
  simp only [Finset.mem_image] at hd
  obtain ⟨p, hp, rfl⟩ := hd
  change p ∈ strictUpper A at hp
  simp only [strictUpper, Finset.mem_filter, Finset.mem_product] at hp
  obtain ⟨⟨hp1, _hp2⟩, hplt⟩ := hp
  have hp1_le : p.1 ≤ N := by
    have := Finset.mem_range.mp (hA hp1)
    omega
  rw [Finset.mem_Icc]
  constructor <;> omega

/-- C2 specialization: the generic injective-channel lemma recovers the
    finite Sidon difference-channel bound on the strict upper triangle. -/
theorem strictUpper_card_le_of_collision_channel (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) :
    (strictUpper A).card ≤ N := by
  have h_channel := card_le_of_injective_channel
    (strictUpper A)
    (Finset.Icc 1 N)
    (fun p : ℕ × ℕ => p.1 - p.2)
    (strictUpper_diff_image_subset_Icc A N hA)
    (diff_injOn_strictUpper A hS)
  simpa [Nat.card_Icc] using h_channel

/-- C2 recovery theorem: the collision-channel abstraction proves the standard
    sharp Sidon difference-count bound for Erdős #30. -/
theorem sidon_collision_channel_bound (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) :
    A.card * (A.card - 1) ≤ 2 * N := by
  have h_upper := strictUpper_card_le_of_collision_channel A N hS hA
  have h_double := strictUpper_card_double A
  omega

end Erdos.Sidon.CollisionChannel
