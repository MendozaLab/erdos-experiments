/-
  Erdős Problem #30 — Lindström Upper Bound (1969)
  =================================================

  **Result:** h(N) ≤ √N + N^{1/4} + 1

  **Reference:** Lindström, B. (1969). "An inequality for B₂-sequences."
  J. Combinatorial Theory, Vol. 7, No. 1, pp. 134-138.

  **Correct proof technique (from Balogh-Füredi-Roy 2023, Section 2):**
  Let A = {a₁ < a₂ < ... < a_k} ⊆ [n] be a Sidon set.

  1. For i < j, call j - i the "order" of the difference a_j - a_i.
  2. Fix parameter ℓ (chosen later as ⌊n^{1/4}⌋).
  3. Count differences of orders 1 through ℓ:
     There are (k-1) + (k-2) + ... + (k-ℓ) = ℓ(k - (ℓ+1)/2) such differences.
  4. These differences are ALL DISTINCT positive integers (Sidon property).
  5. LOWER BOUND: sum ≥ 1 + 2 + ... + m > m²/2 where m = ℓ(k - (ℓ+1)/2).
  6. UPPER BOUND: sum of order-r differences ≤ r·n (telescoping), so
     total ≤ Σ_{r=1}^ℓ r·n = ℓ(ℓ+1)n/2.
  7. Combining: (1/2)ℓ²(k - (ℓ+1)/2)² < (1/2)ℓ(ℓ+1)n
  8. Simplify: ℓ(k - (ℓ+1)/2)² < (ℓ+1)n
  9. k < √(n(ℓ+1)/ℓ) + (ℓ+1)/2 ≈ √n + √n/(2ℓ) + ℓ/2 + 1/2
  10. With ℓ = ⌊n^{1/4}⌋: k < √n + n^{1/4} + 1.

  **This file also contains** a residue-class parametric bound (k² ≤ 2Nt + kt),
  which is an independent result useful for BFR but does NOT directly yield
  the Lindström bound.

  **Superseded by:** Balogh-Füredi-Roy (2023): h(N) ≤ √N + 0.998·N^{1/4}.
  See Erdos30_BFR.lean for that improvement.

  **Lean version**: leanprover/lean4:v4.24.0
  **Mathlib version**: f897ebcf72cd16f89ab4577d0c826cd14afaafc7
-/

import Mathlib
import Erdos30_Sidon_Defs
import Erdos30_Complete

open Finset Nat

namespace Erdos.Sidon

/-! ## Step 3–4: Pigeonhole on Residue Classes -/

/-- Residue class: elements of A congruent to r mod t. -/
def residueClass (A : Finset ℕ) (t r : ℕ) : Finset ℕ :=
  A.filter (fun a => a % t = r)

/-- Pigeonhole: some residue class mod t has ≥ ⌈k/t⌉ elements.
    Stated as: t * |A_r| ≥ |A| for some r. -/
theorem pigeonhole_residue (A : Finset ℕ) (t : ℕ) (ht : 0 < t) :
    ∃ r, r < t ∧ A.card ≤ t * (residueClass A t r).card := by
  by_contra h
  push_neg at h
  -- h : ∀ r, r < t → t * (residueClass A t r).card < A.card
  -- Partition identity: |A| = ∑_{r<t} |A_r|
  have hfib : A.card = ∑ r ∈ Finset.range t, (residueClass A t r).card := by
    simp only [residueClass]
    exact Finset.card_eq_sum_card_fiberwise
      (by intro a _; exact Finset.mem_range.mpr (Nat.mod_lt a ht))
  have hlt : ∀ r ∈ Finset.range t, t * (residueClass A t r).card < A.card :=
    fun r hr => h r (Finset.mem_range.mp hr)
  -- Sum: t * |A| = ∑ t*|A_r| < ∑ |A| = t * |A|. Contradiction.
  have : t * A.card < t * A.card := calc
    t * A.card
        = t * ∑ r ∈ Finset.range t, (residueClass A t r).card := by rw [hfib]
      _ = ∑ r ∈ Finset.range t, t * (residueClass A t r).card := by
          rw [Finset.mul_sum]
      _ < ∑ r ∈ Finset.range t, A.card := by
          apply Finset.sum_lt_sum
          · intro r hr; exact le_of_lt (hlt r hr)
          · exact ⟨0, Finset.mem_range.mpr (by omega),
              hlt 0 (Finset.mem_range.mpr (by omega))⟩
      _ = t * A.card := by rw [Finset.sum_const, Finset.card_range, smul_eq_mul]
  omega

/-! ## Step 5: Residue Class Preserves Sidon Property -/

/-- Scaling a residue class by 1/t preserves the Sidon property.
    If A is Sidon and B = {(a-r)/t : a ∈ A, a ≡ r (mod t)}, then B is Sidon. -/
theorem residue_class_scaled_sidon (A : Finset ℕ) (t r : ℕ) (ht : 0 < t)
    (hS : IsSidonSet A) :
    IsSidonSet ((residueClass A t r).image (fun a => a / t)) := by
  intro a₁ ha₁ b₁ hb₁ a₂ ha₂ b₂ hb₂ hab₁ hab₂ heq
  simp only [Finset.mem_image, residueClass, Finset.mem_filter] at ha₁ hb₁ ha₂ hb₂
  obtain ⟨x₁, ⟨hx₁A, hx₁r⟩, rfl⟩ := ha₁
  obtain ⟨y₁, ⟨hy₁A, hy₁r⟩, rfl⟩ := hb₁
  obtain ⟨x₂, ⟨hx₂A, hx₂r⟩, rfl⟩ := ha₂
  obtain ⟨y₂, ⟨hy₂A, hy₂r⟩, rfl⟩ := hb₂
  -- Key: x_i = t * (x_i / t) + x_i % t = t * (x_i / t) + r
  have hx₁d := Nat.div_add_mod x₁ t
  have hy₁d := Nat.div_add_mod y₁ t
  have hx₂d := Nat.div_add_mod x₂ t
  have hy₂d := Nat.div_add_mod y₂ t
  -- Multiply heq by t: t*(x₁/t) + t*(y₁/t) = t*(x₂/t) + t*(y₂/t)
  have sum_t : t * (x₁ / t) + t * (y₁ / t) = t * (x₂ / t) + t * (y₂ / t) := by
    have := congr_arg (t * ·) heq; simp only [mul_add] at this; exact this
  -- Therefore x₁ + y₁ = x₂ + y₂
  have sum_eq : x₁ + y₁ = x₂ + y₂ := by omega
  -- Ordering: x₁/t ≤ y₁/t implies x₁ ≤ y₁ (same remainder)
  have ord₁ : x₁ ≤ y₁ := by
    have := mul_le_mul_of_nonneg_left hab₁ (Nat.zero_le t); omega
  have ord₂ : x₂ ≤ y₂ := by
    have := mul_le_mul_of_nonneg_left hab₂ (Nat.zero_le t); omega
  -- Apply Sidon property of A
  have hs := hS x₁ hx₁A y₁ hy₁A x₂ hx₂A y₂ hy₂A ord₁ ord₂ sum_eq
  exact ⟨by rw [hs.1], by rw [hs.2]⟩

/-! ## Step 5b: Elementary Sidon Bound — CLOSED 2026-05-02

  Previously declared as axiom; now proved by direct bridge to
  `Erdos30_Complete.sidon_difference_count`, which has a full Lean 4 proof
  using only difference injectivity + pigeonhole on Finset.image.
  Mathlib v4.27.0 API drift (simp lemma changes, ▸ motive issues,
  Finset.pair_comm decidable inference) was patched in `Erdos30_Complete.lean`
  on the same date so both files compile against the same Mathlib pin.
-/

/-- **Elementary Sidon bound**: For Sidon A ⊆ {0,...,M}, |A|*(|A|-1) ≤ 2*M.
    Proved in `Erdos30_Complete.sidon_difference_count`. -/
theorem sidon_elem_bound (A : Finset ℕ) (M : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (M + 1)) : A.card * (A.card - 1) ≤ 2 * M :=
  sidon_difference_count A M hS hA

/-! ## Step 6: Elementary Bound on Scaled Set -/

/-- The scaled residue class B lives in {0,...,⌊N/t⌋}. -/
theorem scaled_range (A : Finset ℕ) (N t r : ℕ) (ht : 0 < t)
    (hA : A ⊆ Finset.range (N + 1)) :
    (residueClass A t r).image (fun a => a / t) ⊆ Finset.range (N / t + 1) := by
  intro b hb
  simp [Finset.mem_image, residueClass, Finset.mem_filter] at hb
  obtain ⟨a, ⟨ha, _⟩, rfl⟩ := hb
  simp [Finset.mem_range]
  have : a ≤ N := by
    have := Finset.mem_range.mp (hA ha); omega
  exact Nat.div_le_div_right this

/-! ## Step 7–8: The Key Inequality

  From steps 4–6:
  - |B| ≥ k/t (pigeonhole)
  - B is Sidon in {0,...,N/t} (residue scaling)
  - |B|(|B|-1) ≤ 2·(N/t) (elementary bound)

  Combining: (k/t)(k/t - 1) ≤ 2N/t
  Multiply by t²: k(k-t) ≤ 2Nt ... but we need to be careful with integer division.

  Working with integers: let m = |A_r| ≥ ⌈k/t⌉, so m·t ≥ k.
  m(m-1) ≤ 2·⌊N/t⌋.
  m ≤ t/2 + √(2N/t + t/4)  (complete the square).
  k ≤ m·t ≤ ...

  Actually, cleaner to state the parametric inequality directly.
-/

/-- **Parametric Lindström inequality.**
  For Sidon A ⊆ {0,...,N} and any t > 0:
    k² ≤ 2Nt + kt
  where k = |A|.

  This is the key step — everything else is optimization over t. -/
theorem lindstrom_parametric (A : Finset ℕ) (N t : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (ht : 0 < t) :
    A.card * A.card ≤ 2 * N * t + A.card * t := by
  -- Step 1: Pigeonhole gives residue class with k ≤ t * m elements
  obtain ⟨r, hr, hpig⟩ := pigeonhole_residue A t ht
  set Ar := residueClass A t r with hAr_def
  set m := Ar.card
  set B := Ar.image (fun a => a / t)
  -- Step 2: Division is injective on residue class (all have same remainder)
  have h_inj : Set.InjOn (fun a => a / t) (↑Ar) := by
    intro a ha b hb hab
    -- Extract filter condition: elements of Ar have remainder r mod t
    have ha_r : a % t = r := by
      have := Finset.mem_coe.mp ha
      rw [hAr_def, residueClass, Finset.mem_filter] at this
      exact this.2
    have hb_r : b % t = r := by
      have := Finset.mem_coe.mp hb
      rw [hAr_def, residueClass, Finset.mem_filter] at this
      exact this.2
    -- a = t*(a/t) + r, b = t*(b/t) + r, a/t = b/t → a = b
    have hab' : a / t = b / t := hab
    calc a = t * (a / t) + a % t := (Nat.div_add_mod a t).symm
      _ = t * (b / t) + r := by rw [hab', ha_r]
      _ = t * (b / t) + b % t := by rw [hb_r]
      _ = b := Nat.div_add_mod b t
  have hBcard : B.card = m := Finset.card_image_of_injOn h_inj
  -- Step 3: B is Sidon in {0,...,N/t}
  have hBS : IsSidonSet B := residue_class_scaled_sidon A t r ht hS
  have hBrange : B ⊆ Finset.range (N / t + 1) := scaled_range A N t r ht hA
  -- Step 4: Elementary bound → m*(m-1) ≤ 2*(N/t)
  have hBbound : m * (m - 1) ≤ 2 * (N / t) := by
    have := sidon_elem_bound B (N / t) hBS hBrange; rwa [hBcard] at this
  -- Step 5: Case split on k ≤ t (trivial) vs k > t
  by_cases hkt : A.card ≤ t
  · -- k ≤ t: k² ≤ kt ≤ 2Nt + kt
    calc A.card * A.card ≤ A.card * t := Nat.mul_le_mul_left _ hkt
       _ ≤ 2 * N * t + A.card * t := Nat.le_add_left _ _
  · -- k > t: decompose k² = k*(k-t) + k*t, bound k*(k-t) ≤ 2*N*t
    push_neg at hkt
    have hge : t ≤ A.card := Nat.le_of_lt hkt
    -- k² = k*(k-t) + k*t
    have h_split : A.card * A.card = A.card * (A.card - t) + A.card * t := by
      rw [← Nat.mul_add, Nat.sub_add_cancel hge]
    rw [h_split]
    suffices h : A.card * (A.card - t) ≤ 2 * N * t by linarith
    -- Nat arithmetic: t*m - t = t*(m-1)
    have htm_eq : t * m - t = t * (m - 1) := by
      cases m with
      | zero => simp
      | succ n => simp [Nat.mul_succ]
    -- Chain: k*(k-t) ≤ (t*m)*(t*m-t) = t²*m*(m-1) ≤ 2*(t*(N/t))*t ≤ 2*N*t
    have h4 : t * (N / t) ≤ N := Nat.mul_div_le N t
    calc A.card * (A.card - t)
        ≤ (t * m) * (t * m - t) :=
          Nat.mul_le_mul hpig (Nat.sub_le_sub_right hpig t)
      _ = (t * m) * (t * (m - 1)) := by rw [htm_eq]
      _ = t * t * (m * (m - 1)) := by ring
      _ ≤ t * t * (2 * (N / t)) :=
          Nat.mul_le_mul_left _ hBbound
      _ = 2 * (t * (N / t)) * t := by ring
      _ ≤ 2 * N * t :=
          Nat.mul_le_mul_right t (Nat.mul_le_mul_left 2 h4)

/-! ## The Actual Lindström Proof (order-of-differences argument)

  The parametric inequality k² ≤ 2Nt + kt (above) is useful for BFR but
  does NOT directly yield the Lindström bound (it gives O(N^{5/8})).

  The correct proof (BFR Section 2) uses a completely different technique:
  ordering elements and counting differences by "order" (index gap j - i).

  Key idea: For Sidon A = {a₁ < ... < a_k} ⊆ [N], ALL pairwise differences
  are distinct. The differences of orders 1..ℓ are m = ℓ(k-(ℓ+1)/2)
  distinct positive integers. Their sum is ≥ m²/2 (minimum sum of m distinct
  positive integers) and ≤ ℓ(ℓ+1)N/2 (telescoping: order-r diffs sum to ≤ rN).
  This gives ℓ(k-(ℓ+1)/2)² < (ℓ+1)N, and with ℓ = ⌊N^{1/4}⌋:
    k < √N + N^{1/4} + 1.
-/

/-- For a Sidon set, all pairwise positive differences are distinct.
    This is equivalent to the Sidon (B₂) property. -/
theorem sidon_distinct_differences (A : Finset ℕ) (hS : IsSidonSet A) :
    ∀ a₁ ∈ A, ∀ b₁ ∈ A, ∀ a₂ ∈ A, ∀ b₂ ∈ A,
      a₁ < b₁ → a₂ < b₂ → b₁ - a₁ = b₂ - a₂ → (a₁ = a₂ ∧ b₁ = b₂) := by
  intro a₁ ha₁ b₁ hb₁ a₂ ha₂ b₂ hb₂ h₁ h₂ heq
  have hsum : a₁ + b₂ = a₂ + b₁ := by omega
  clear heq  -- Remove Nat subtraction from context (confuses omega)
  -- 4-way case split on orderings needed for Sidon application
  by_cases hab : a₁ ≤ b₂
  · by_cases hab' : a₂ ≤ b₁
    · -- a₁ ≤ b₂, a₂ ≤ b₁: direct Sidon on (a₁,b₂) vs (a₂,b₁)
      have hs := hS a₁ ha₁ b₂ hb₂ a₂ ha₂ b₁ hb₁ hab hab' hsum
      exact ⟨hs.1, hs.2.symm⟩
    · -- a₁ ≤ b₂, b₁ < a₂: Sidon on (a₁,b₂) vs (b₁,a₂) → a₁=b₁, contradiction
      push_neg at hab'
      exfalso
      have hs := hS a₁ ha₁ b₂ hb₂ b₁ hb₁ a₂ ha₂ hab (by omega) (by omega)
      omega
  · by_cases hab' : a₂ ≤ b₁
    · -- b₂ < a₁, a₂ ≤ b₁: Sidon on (b₂,a₁) vs (a₂,b₁) → a₁=b₁, contradiction
      push_neg at hab
      exfalso
      have hs := hS b₂ hb₂ a₁ ha₁ a₂ ha₂ b₁ hb₁ (by omega) hab' (by omega)
      omega
    · -- b₂ < a₁, b₁ < a₂: chain a₁ < b₁ < a₂ < b₂ < a₁, contradiction
      push_neg at hab hab'
      omega

/-! ## Intermediate lemmas for the order-of-differences proof -/

/-- The set of positive pairwise differences from A. -/
def posDiffs (A : Finset ℕ) : Finset ℕ :=
  ((A ×ˢ A).filter (fun p => p.1 < p.2)).image (fun p => p.2 - p.1)

/-- For a Sidon set, the map (a,b) ↦ b-a is injective on ordered pairs,
    so |posDiffs A| = k(k-1)/2. -/
theorem card_posDiffs_sidon (A : Finset ℕ) (hS : IsSidonSet A) :
    (posDiffs A).card = A.card * (A.card - 1) / 2 := by
  unfold posDiffs
  -- Step 1: Injectivity from Sidon distinct differences
  have h_inj : Set.InjOn (fun p : ℕ × ℕ => p.2 - p.1)
      ↑((A ×ˢ A).filter (fun p => p.1 < p.2)) := by
    intro ⟨a₁, b₁⟩ h₁ ⟨a₂, b₂⟩ h₂ heq
    simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_product] at h₁ h₂
    have := sidon_distinct_differences A hS a₁ h₁.1.1 b₁ h₁.1.2 a₂ h₂.1.1 b₂ h₂.1.2 h₁.2 h₂.2 heq
    exact Prod.ext this.1 this.2
  rw [Finset.card_image_of_injOn h_inj]
  -- Step 2: filter(< on A×A) = filter(< on offDiag)
  have h_filter_eq : (A ×ˢ A).filter (fun p : ℕ × ℕ => p.1 < p.2) =
      A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2) := by
    ext ⟨a, b⟩
    simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_offDiag]
    constructor
    · intro ⟨⟨ha, hb⟩, hab⟩; exact ⟨⟨ha, hb, Nat.ne_of_lt hab⟩, hab⟩
    · intro ⟨⟨ha, hb, _⟩, hab⟩; exact ⟨⟨ha, hb⟩, hab⟩
  rw [h_filter_eq]
  -- Step 3: Partition offDiag into filter(<) ∪ filter(>)
  have h_union : A.offDiag =
      A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2) ∪
      A.offDiag.filter (fun p : ℕ × ℕ => p.2 < p.1) := by
    ext ⟨a, b⟩
    simp only [Finset.mem_offDiag, Finset.mem_union, Finset.mem_filter]
    constructor
    · intro ⟨ha, hb, hab⟩
      rcases lt_or_gt_of_ne hab with h | h
      · left; exact ⟨⟨ha, hb, hab⟩, h⟩
      · right; exact ⟨⟨ha, hb, hab⟩, h⟩
    · rintro (⟨h, _⟩ | ⟨h, _⟩) <;> exact h
  have h_disj : Disjoint
      (A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2))
      (A.offDiag.filter (fun p : ℕ × ℕ => p.2 < p.1)) :=
    Finset.disjoint_filter.mpr (fun ⟨a, b⟩ _ h1 h2 => absurd h1 (not_lt.mpr (le_of_lt h2)))
  -- Step 4: Swap bijection |filter(>)| = |filter(<)|
  have h_swap : (A.offDiag.filter (fun p : ℕ × ℕ => p.2 < p.1)).card =
      (A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2)).card :=
    Finset.card_bij' (fun p _ => (p.2, p.1)) (fun p _ => (p.2, p.1))
      (fun ⟨a, b⟩ h => by
        simp only [Finset.mem_filter, Finset.mem_offDiag] at h ⊢
        exact ⟨⟨h.1.2.1, h.1.1, Ne.symm h.1.2.2⟩, h.2⟩)
      (fun ⟨a, b⟩ h => by
        simp only [Finset.mem_filter, Finset.mem_offDiag] at h ⊢
        exact ⟨⟨h.1.2.1, h.1.1, Ne.symm h.1.2.2⟩, h.2⟩)
      (fun _ _ => rfl) (fun _ _ => rfl)
  -- Step 5: Combine cardinalities
  have h_card : A.offDiag.card =
      (A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2)).card +
      (A.offDiag.filter (fun p : ℕ × ℕ => p.2 < p.1)).card := by
    rw [← Finset.card_union_of_disjoint h_disj, ← h_union]
  rw [h_swap] at h_card
  -- h_card : |offDiag| = 2 * |filter(<)|
  -- offDiag_card : |offDiag| = k*k - k
  have h_mul_sub : A.card * (A.card - 1) = A.card * A.card - A.card := by
    cases A.card with
    | zero => simp
    | succ n =>
      simp only [Nat.succ_sub_one]
      rw [show (n + 1) * (n + 1) = (n + 1) * n + (n + 1) from by ring, Nat.add_sub_cancel]
  have h_offDiag_eq : A.offDiag.card = A.card * (A.card - 1) :=
    A.offDiag_card.trans h_mul_sub.symm
  -- 2 * |filter(<)| = k*(k-1)
  have h_2c : 2 * (A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2)).card =
      A.card * (A.card - 1) := by linarith
  -- Final: |filter(<)| = k*(k-1)/2
  have : A.card * (A.card - 1) / 2 =
      (A.offDiag.filter (fun p : ℕ × ℕ => p.1 < p.2)).card := by
    rw [← h_2c]; exact Nat.mul_div_cancel_left _ (by omega)
  exact this.symm

/-- All positive differences of A ⊆ {0,...,N} lie in {1,...,N}. -/
theorem posDiffs_subset_Icc (A : Finset ℕ) (N : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    ∀ d ∈ posDiffs A, d ≤ N := by
  intro d hd
  simp only [posDiffs, Finset.mem_image, Finset.mem_filter, Finset.mem_product] at hd
  obtain ⟨⟨a, b⟩, ⟨⟨ha, hb⟩, hab⟩, rfl⟩ := hd
  have : b ≤ N := by have := Finset.mem_range.mp (hA hb); omega
  omega

/-! ### Helper: Sum of distinct positive naturals -/

/-- Auxiliary: For Finset S ⊆ ℕ with |S| = m and all elements ≥ c,
    m² + 2mc ≤ 2·sum(S) + m.
    Equivalently: 2·sum(S) ≥ m(m-1) + 2mc (but stated without ℕ subtraction).
    Used with c = 1 for the Lindström lower bound. -/
private lemma sum_distinct_ge_aux : ∀ (m : ℕ) (S : Finset ℕ) (c : ℕ),
    S.card = m → (∀ x ∈ S, c ≤ x) →
    m * m + 2 * m * c ≤ 2 * S.sum id + m := by
  intro m
  induction m with
  | zero => intro S c hm _; simp [show S = ∅ from Finset.card_eq_zero.mp hm]
  | succ n ih =>
    intro S c hm hc
    have hne : S.Nonempty := Finset.card_pos.mp (by omega)
    set x := S.min' hne
    have hx_mem : x ∈ S := Finset.min'_mem S hne
    have hx_ge : c ≤ x := hc x hx_mem
    set S' := S.erase x
    have hS'_card : S'.card = n := by
      rw [Finset.card_erase_of_mem hx_mem, hm]; omega
    have hS'_bound : ∀ y ∈ S', c + 1 ≤ y := by
      intro y hy
      have hy_mem : y ∈ S := Finset.mem_of_mem_erase hy
      have hy_ne : y ≠ x := Finset.ne_of_mem_erase hy
      have := Finset.min'_le S y hy_mem
      omega
    have h_ih := ih S' (c + 1) hS'_card hS'_bound
    -- h_ih : n * n + 2 * n * (c + 1) ≤ 2 * S'.sum id + n
    have h_split : S.sum id = (S.erase x).sum id + x :=
      (Finset.sum_erase_add S id hx_mem).symm
    -- Expand nonlinear terms for linarith
    have e1 : (n + 1) * (n + 1) = n * n + 2 * n + 1 := by ring
    have e2 : 2 * (n + 1) * c = 2 * n * c + 2 * c := by ring
    have e3 : 2 * n * (c + 1) = 2 * n * c + 2 * n := by ring
    rw [e3] at h_ih
    rw [e1, e2, h_split]
    -- Goal: n*n + 2*n + 1 + (2*n*c + 2*c) ≤ 2*(S'.sum id + x) + (n+1)
    -- From h_ih: n*n + 2*n*c + 2*n ≤ 2*S'.sum id + n
    -- From hx_ge: c ≤ x
    linarith

/-- Sum of m distinct positive naturals is at least m(m+1)/2.
    Stated as m(m+1) ≤ 2·sum to avoid ℕ division. -/
lemma sum_distinct_pos_ge (S : Finset ℕ) (h_pos : ∀ x ∈ S, 0 < x) :
    S.card * (S.card + 1) ≤ 2 * S.sum id := by
  have h := sum_distinct_ge_aux S.card S 1 rfl (fun x hx => h_pos x hx)
  -- h : S.card * S.card + 2 * S.card * 1 ≤ 2 * S.sum id + S.card
  have e : S.card * (S.card + 1) = S.card * S.card + S.card := by ring
  linarith

/-! ### Combinatorial core: counting and summing order-bounded differences -/

/-! #### Sorted-enumeration helpers (private; localized to avoid stale OrderedElements) -/

/-- The sorted enumeration of A as a list. -/
private noncomputable def oList (A : Finset ℕ) : List ℕ := A.sort (· ≤ ·)

private lemma length_oList (A : Finset ℕ) : (oList A).length = A.card := by
  unfold oList; rw [Finset.length_sort]

private lemma sortedLT_oList (A : Finset ℕ) : StrictMono (oList A).get :=
  Finset.sortedLT_sort A

/-- The i-th element in the sorted enumeration; default 0 outside bounds. -/
private noncomputable def oGet (A : Finset ℕ) (i : ℕ) : ℕ := (oList A).getD i 0

private lemma oGet_eq_get {A : Finset ℕ} {i : ℕ} (h : i < A.card) :
    oGet A i = (oList A).get ⟨i, by rw [length_oList]; exact h⟩ := by
  unfold oGet
  have hl : i < (oList A).length := by rw [length_oList]; exact h
  exact List.getD_eq_getElem _ 0 hl

private lemma oGet_mem {A : Finset ℕ} {i : ℕ} (h : i < A.card) : oGet A i ∈ A := by
  rw [oGet_eq_get h]
  unfold oList
  apply (Finset.mem_sort (· ≤ ·)).mp
  exact List.get_mem _ _

private lemma oGet_lt {A : Finset ℕ} {i j : ℕ} (hi : i < A.card) (hj : j < A.card)
    (hij : i < j) : oGet A i < oGet A j := by
  rw [oGet_eq_get hi, oGet_eq_get hj]
  exact sortedLT_oList A hij

private lemma oGet_le {A : Finset ℕ} {i j : ℕ} (hi : i < A.card) (hj : j < A.card)
    (hij : i ≤ j) : oGet A i ≤ oGet A j := by
  rcases lt_or_eq_of_le hij with h | h
  · exact le_of_lt (oGet_lt hi hj h)
  · subst h; rfl

private lemma oGet_le_N {A : Finset ℕ} {N : ℕ} (hA : A ⊆ Finset.range (N + 1))
    {i : ℕ} (hi : i < A.card) : oGet A i ≤ N := by
  have h_mem : oGet A i ∈ A := oGet_mem hi
  have := Finset.mem_range.mp (hA h_mem)
  omega

/-! #### orderPairs and orderDiffs -/

/-- The set of index pairs (i, j) with i < j ≤ i + ℓ, both indices < A.card. -/
private noncomputable def orderPairs (A : Finset ℕ) (ℓ : ℕ) : Finset (ℕ × ℕ) :=
  ((Finset.range A.card) ×ˢ (Finset.range A.card)).filter
    (fun p => p.1 < p.2 ∧ p.2 ≤ p.1 + ℓ)

private lemma mem_orderPairs_iff (A : Finset ℕ) (ℓ : ℕ) (p : ℕ × ℕ) :
    p ∈ orderPairs A ℓ ↔ p.1 < A.card ∧ p.2 < A.card ∧ p.1 < p.2 ∧ p.2 ≤ p.1 + ℓ := by
  unfold orderPairs
  simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_range]
  tauto

/-- The set of differences a_j - a_i for (i, j) ∈ orderPairs. -/
private noncomputable def orderDiffs (A : Finset ℕ) (ℓ : ℕ) : Finset ℕ :=
  (orderPairs A ℓ).image (fun p => oGet A p.2 - oGet A p.1)

/-! #### Cardinality of orderPairs and orderDiffs -/

/-- Auxiliary: integer arithmetic for cardinality formula.
    2 * Σ_{r=0}^{ℓ-1} (k - (r+1)) = ℓ * (2k - ℓ - 1) when ℓ ≤ k. -/
private lemma sum_int_offset (ℓ : ℕ) (k : ℤ) :
    (∑ r ∈ Finset.range ℓ, (k - r - 1)) = ℓ * k - (∑ r ∈ Finset.range ℓ, (r : ℤ)) - ℓ := by
  induction ℓ with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, Finset.sum_range_succ, ih]
    push_cast
    ring

private lemma cardSum_doubled (k ℓ : ℕ) (hℓ : ℓ ≤ k) :
    2 * (∑ r ∈ Finset.range ℓ, (k - (r + 1))) = ℓ * (2 * k - ℓ - 1) := by
  rcases Nat.eq_zero_or_pos ℓ with hℓ0 | hℓpos
  · subst hℓ0; simp
  have h2kℓ : ℓ + 1 ≤ 2 * k := by omega
  zify [show 2 * k ≥ ℓ + 1 from by omega, hℓ]
  have h1 : ∀ r ∈ Finset.range ℓ, ((k - (r + 1) : ℕ) : ℤ) = (k : ℤ) - r - 1 := by
    intro r hr
    rw [Finset.mem_range] at hr
    omega
  rw [Finset.sum_congr rfl h1]
  rw [show ((2 * k - ℓ - 1 : ℕ) : ℤ) = 2 * (k : ℤ) - ℓ - 1 from by omega]
  rw [sum_int_offset ℓ (k : ℤ)]
  have hS : (∑ r ∈ Finset.range ℓ, (r : ℤ)) * 2 = (ℓ : ℤ) * (ℓ - 1) := by
    have := Finset.sum_range_id_mul_two ℓ
    zify at this
    have h_lift : ((ℓ - 1 : ℕ) : ℤ) = (ℓ : ℤ) - 1 := by omega
    rw [h_lift] at this
    exact this
  linarith

private lemma orderPairs_card_eq (A : Finset ℕ) (ℓ : ℕ) :
    (orderPairs A ℓ).card = ∑ r ∈ Finset.range ℓ, (A.card - (r + 1)) := by
  rw [show (∑ r ∈ Finset.range ℓ, (A.card - (r + 1))) =
      ((Finset.range ℓ).sigma (fun r => Finset.range (A.card - (r + 1)))).card from ?_]
  · apply Finset.card_bij (fun (p : ℕ × ℕ) _ => Sigma.mk (p.2 - p.1 - 1) p.1)
    · intro ⟨i, j⟩ hp
      simp only [orderPairs, Finset.mem_filter, Finset.mem_product, Finset.mem_range] at hp
      simp only [Finset.mem_sigma, Finset.mem_range]
      exact ⟨by omega, by omega⟩
    · intro ⟨i1, j1⟩ hp1 ⟨i2, j2⟩ hp2 heq
      simp only [orderPairs, Finset.mem_filter, Finset.mem_product, Finset.mem_range] at hp1 hp2
      simp only [Sigma.mk.injEq] at heq
      obtain ⟨_, hieq⟩ := heq
      subst hieq
      have h1 : i1 < j1 := hp1.2.1
      have h2 : i1 < j2 := hp2.2.1
      have hj : j1 = j2 := by omega
      simp [hj]
    · intro ⟨r, i⟩ hri
      simp only [Finset.mem_sigma, Finset.mem_range] at hri
      obtain ⟨hr, hi_bound⟩ := hri
      refine ⟨(i, i + r + 1), ?_, ?_⟩
      · simp only [orderPairs, Finset.mem_filter, Finset.mem_product, Finset.mem_range]
        exact ⟨⟨by omega, by omega⟩, by omega, by omega⟩
      · ext
        · simp; omega
        · simp
  · rw [Finset.card_sigma]
    congr 1
    ext r
    rw [Finset.card_range]

private lemma orderPairs_card_doubled (A : Finset ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ A.card) :
    (orderPairs A ℓ).card * 2 = ℓ * (2 * A.card - ℓ - 1) := by
  rw [orderPairs_card_eq, mul_comm]
  exact cardSum_doubled A.card ℓ hℓ

/-- Sidon injectivity: (i, j) → a_j - a_i is injective on orderPairs. -/
private lemma orderDiffs_injOn (A : Finset ℕ) (ℓ : ℕ) (hS : IsSidonSet A) :
    Set.InjOn (fun p : ℕ × ℕ => oGet A p.2 - oGet A p.1) (orderPairs A ℓ) := by
  intro ⟨i1, j1⟩ hp1 ⟨i2, j2⟩ hp2 heq
  rw [Finset.mem_coe, mem_orderPairs_iff] at hp1 hp2
  obtain ⟨hi1, hj1, hij1, _⟩ := hp1
  obtain ⟨hi2, hj2, hij2, _⟩ := hp2
  -- oGet A i_k ∈ A, and oGet A i_k < oGet A j_k
  have hi1_lt : oGet A i1 < oGet A j1 := oGet_lt hi1 hj1 hij1
  have hi2_lt : oGet A i2 < oGet A j2 := oGet_lt hi2 hj2 hij2
  have ha_i1 : oGet A i1 ∈ A := oGet_mem hi1
  have ha_j1 : oGet A j1 ∈ A := oGet_mem hj1
  have ha_i2 : oGet A i2 ∈ A := oGet_mem hi2
  have ha_j2 : oGet A j2 ∈ A := oGet_mem hj2
  have h_diff_eq := sidon_distinct_differences A hS
    (oGet A i1) ha_i1 (oGet A j1) ha_j1
    (oGet A i2) ha_i2 (oGet A j2) ha_j2 hi1_lt hi2_lt heq
  -- h_diff_eq : oGet A i1 = oGet A i2 ∧ oGet A j1 = oGet A j2
  -- Need: i1 = i2 and j1 = j2.
  -- oGet is strict mono, so injective on indices < A.card.
  ext
  · -- (i1, j1).1 = (i2, j2).1, i.e., i1 = i2
    have hii : oGet A i1 = oGet A i2 := h_diff_eq.1
    by_contra hne
    rcases Nat.lt_or_ge i1 i2 with h12 | h21
    · exact absurd hii (Nat.ne_of_lt (oGet_lt hi1 hi2 h12))
    · have h21' : i2 ≤ i1 := h21
      rcases lt_or_eq_of_le h21' with h | h
      · exact absurd hii (Nat.ne_of_gt (oGet_lt hi2 hi1 h))
      · exact hne h.symm
  · -- j1 = j2
    have hjj : oGet A j1 = oGet A j2 := h_diff_eq.2
    by_contra hne
    rcases Nat.lt_or_ge j1 j2 with h12 | h21
    · exact absurd hjj (Nat.ne_of_lt (oGet_lt hj1 hj2 h12))
    · have h21' : j2 ≤ j1 := h21
      rcases lt_or_eq_of_le h21' with h | h
      · exact absurd hjj (Nat.ne_of_gt (oGet_lt hj2 hj1 h))
      · exact hne h.symm

private lemma orderDiffs_card_eq (A : Finset ℕ) (ℓ : ℕ) (hS : IsSidonSet A) :
    (orderDiffs A ℓ).card = (orderPairs A ℓ).card := by
  unfold orderDiffs
  exact Finset.card_image_of_injOn (orderDiffs_injOn A ℓ hS)

private lemma orderDiffs_pos (A : Finset ℕ) (ℓ : ℕ) :
    ∀ d ∈ orderDiffs A ℓ, 0 < d := by
  intro d hd
  unfold orderDiffs at hd
  simp only [Finset.mem_image] at hd
  obtain ⟨p, hp, rfl⟩ := hd
  rw [mem_orderPairs_iff] at hp
  obtain ⟨hp1, hp2, hij, _⟩ := hp
  have h_lt : oGet A p.1 < oGet A p.2 := oGet_lt hp1 hp2 hij
  omega

/-! #### Sum bound: 2 * Σ orderDiffs ≤ ℓ * (ℓ + 1) * N

    Strategy: decompose orderPairs by row r := j - i - 1 ∈ [0, ℓ-1]; for fixed r,
    the row sum Σ_{i=0}^{k-r-2} (a_{i+r+1} - a_i) telescopes to ≤ (r+1) * N. -/

/-- Row r of orderPairs: pairs (i, i + r + 1) with i + r + 1 < A.card. -/
private noncomputable def rowPairs (A : Finset ℕ) (r : ℕ) : Finset (ℕ × ℕ) :=
  (Finset.range (A.card - (r + 1))).image (fun i => (i, i + r + 1))

private lemma mem_rowPairs (A : Finset ℕ) (r : ℕ) (p : ℕ × ℕ) :
    p ∈ rowPairs A r ↔ p.1 + r + 1 < A.card ∧ p.2 = p.1 + r + 1 := by
  unfold rowPairs
  simp only [Finset.mem_image, Finset.mem_range]
  constructor
  · rintro ⟨i, hi, rfl⟩
    refine ⟨by omega, rfl⟩
  · rintro ⟨h1, h2⟩
    refine ⟨p.1, by omega, ?_⟩
    ext
    · simp
    · simp [h2]

private lemma orderPairs_eq_biUnion_rowPairs (A : Finset ℕ) (ℓ : ℕ) :
    orderPairs A ℓ = (Finset.range ℓ).biUnion (fun r => rowPairs A r) := by
  ext ⟨i, j⟩
  simp only [Finset.mem_biUnion, Finset.mem_range]
  rw [mem_orderPairs_iff]
  constructor
  · rintro ⟨hi, hj, hij, hjℓ⟩
    refine ⟨j - i - 1, by omega, ?_⟩
    rw [mem_rowPairs]
    refine ⟨by omega, by omega⟩
  · rintro ⟨r, hr, hrow⟩
    rw [mem_rowPairs] at hrow
    obtain ⟨h1, h2⟩ := hrow
    -- p.1 + r + 1 < A.card, p.2 = p.1 + r + 1
    -- Goal: (i, j).1 < A.card ∧ (i, j).2 < A.card ∧ (i, j).1 < (i, j).2 ∧ ...
    refine ⟨by omega, by omega, by omega, by omega⟩

private lemma rowPairs_disjoint (A : Finset ℕ) :
    ∀ {r1 r2 : ℕ}, r1 ≠ r2 → Disjoint (rowPairs A r1) (rowPairs A r2) := by
  intro r1 r2 hne
  rw [Finset.disjoint_left]
  intro ⟨i, j⟩ h1 h2
  rw [mem_rowPairs] at h1 h2
  obtain ⟨_, h1eq⟩ := h1
  obtain ⟨_, h2eq⟩ := h2
  apply hne
  omega

private lemma rowPairs_pairwiseDisjoint (A : Finset ℕ) (ℓ : ℕ) :
    Set.PairwiseDisjoint (Finset.range ℓ : Set ℕ) (fun r => rowPairs A r) := by
  intro r1 _ r2 _ hne
  exact rowPairs_disjoint A hne

/-- Sum of orderPairs decomposes as a sum over rows. -/
private lemma sum_orderPairs_eq_sum_rowPairs (A : Finset ℕ) (ℓ : ℕ) (f : ℕ × ℕ → ℤ) :
    (∑ p ∈ orderPairs A ℓ, f p) = ∑ r ∈ Finset.range ℓ, ∑ p ∈ rowPairs A r, f p := by
  rw [orderPairs_eq_biUnion_rowPairs A ℓ]
  exact Finset.sum_biUnion (rowPairs_pairwiseDisjoint A ℓ)

/-- Row r of rowPairs sum telescopes (in ℤ) to a_{k-1} - a_? for the r=0 case,
    or a difference of partial sums in general. We just need an upper bound. -/
private lemma rowPairs_sum_bound (A : Finset ℕ) (N : ℕ) (r : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    (∑ p ∈ rowPairs A r, ((oGet A p.2 : ℤ) - oGet A p.1)) ≤ (r + 1) * N := by
  -- rowPairs A r = {(i, i+r+1) : i + r + 1 < A.card}
  -- Sum: Σ_i (a_{i+r+1} - a_i) where i ranges over [0, A.card - r - 2]
  unfold rowPairs
  rw [Finset.sum_image (by
    intro i _ j _ heq
    simp only [Prod.mk.injEq] at heq
    exact heq.1)]
  -- Goal: ∑ i ∈ range (A.card - (r+1)), (oGet A (i+r+1) - oGet A i) ≤ (r+1) * N
  by_cases h_empty : A.card ≤ r + 1
  · -- range is empty, sum is 0, RHS ≥ 0.
    rw [show A.card - (r + 1) = 0 from by omega]
    simp
  push_neg at h_empty
  -- A.card > r + 1, so A.card - (r+1) > 0
  -- Telescoping: a_{i+r+1} - a_i = Σ_{s=0}^{r} (a_{i+s+1} - a_{i+s})
  have h_telescope : ∀ i, i + r + 1 < A.card →
      ((oGet A (i + r + 1) : ℤ) - oGet A i) =
      ∑ s ∈ Finset.range (r + 1), ((oGet A (i + s + 1) : ℤ) - oGet A (i + s)) := by
    intro i _
    have h_eq := Finset.sum_range_sub (fun s => (oGet A (i + s) : ℤ)) (r + 1)
    -- h_eq : ∑ s in range (r+1), (oGet A (i+(s+1)) - oGet A (i+s)) = oGet A (i+(r+1)) - oGet A (i+0)
    have h_eq' : ∑ s ∈ Finset.range (r + 1), ((oGet A (i + s + 1) : ℤ) - oGet A (i + s)) =
                 (oGet A (i + (r + 1)) : ℤ) - oGet A (i + 0) := by
      rw [← h_eq]
      apply Finset.sum_congr rfl
      intro s _
      have : i + s + 1 = i + (s + 1) := by ring
      rw [this]
    rw [h_eq']
    have hr1 : i + (r + 1) = i + r + 1 := by ring
    have hi0 : i + 0 = i := by ring
    rw [hr1, hi0]
  rw [Finset.sum_congr rfl (fun i hi => h_telescope i (by
    rw [Finset.mem_range] at hi; omega))]
  rw [Finset.sum_comm]
  -- Goal: ∑ s in range (r+1), ∑ i in range (A.card - (r+1)), (a_{i+s+1} - a_{i+s}) ≤ (r+1) * N
  -- Inner sum (over i): telescopes to a_{(A.card - r - 1) + s} - a_s ≤ N (in ℤ)
  have h_inner : ∀ s ∈ Finset.range (r + 1),
      (∑ i ∈ Finset.range (A.card - (r + 1)), ((oGet A (i + s + 1) : ℤ) - oGet A (i + s))) ≤ N := by
    intro s hs
    rw [Finset.mem_range] at hs
    -- Σ_{i=0}^{m-1} (a_{i+s+1} - a_{i+s}) = a_{m+s} - a_s (telescoping with shift s)
    let m := A.card - (r + 1)
    have h_telescope_inner : ∀ M, (∑ i ∈ Finset.range M, ((oGet A (i + s + 1) : ℤ) - oGet A (i + s))) =
        oGet A (M + s) - oGet A s := by
      intro M
      induction M with
      | zero => simp
      | succ M' ih =>
        rw [Finset.sum_range_succ, ih]
        have : (M' + s + 1 : ℕ) = (M' + 1 + s : ℕ) := by ring
        rw [this]
        ring
    rw [h_telescope_inner]
    -- m + s = A.card - (r+1) + s. Since s ≤ r, m + s ≤ A.card - 1.
    have h_idx : m + s < A.card := by
      simp only [m]
      omega
    have h_aN : oGet A (m + s) ≤ N := oGet_le_N hA h_idx
    have h_a0 : 0 ≤ (oGet A s : ℤ) := Int.natCast_nonneg _
    -- m + s = A.card - (r+1) + s = A.card - 1 - (r - s); always ≤ A.card - 1
    -- oGet A (m + s) ≤ N
    have : (oGet A (m + s) : ℤ) ≤ N := by exact_mod_cast h_aN
    linarith
  -- Now apply outer sum bound
  calc (∑ s ∈ Finset.range (r + 1), ∑ i ∈ Finset.range (A.card - (r + 1)),
          ((oGet A (i + s + 1) : ℤ) - oGet A (i + s)))
      ≤ ∑ s ∈ Finset.range (r + 1), (N : ℤ) := by
        apply Finset.sum_le_sum h_inner
    _ = (r + 1) * N := by
        rw [Finset.sum_const, Finset.card_range]
        rw [nsmul_eq_mul]
        push_cast; ring

private lemma orderPairs_sum_bound (A : Finset ℕ) (N : ℕ) (ℓ : ℕ)
    (hA : A ⊆ Finset.range (N + 1)) :
    (∑ p ∈ orderPairs A ℓ, ((oGet A p.2 : ℤ) - oGet A p.1)) ≤ (ℓ * (ℓ + 1) / 2) * N := by
  rw [sum_orderPairs_eq_sum_rowPairs A ℓ (fun p => (oGet A p.2 : ℤ) - oGet A p.1)]
  -- Σ_{r=0}^{ℓ-1} rowPairs.sum ≤ Σ_r (r+1) * N = N * ℓ(ℓ+1)/2
  have h_each : ∀ r ∈ Finset.range ℓ,
      (∑ p ∈ rowPairs A r, ((oGet A p.2 : ℤ) - oGet A p.1)) ≤ (r + 1) * N := fun r _ =>
    rowPairs_sum_bound A N r hA
  calc (∑ r ∈ Finset.range ℓ, ∑ p ∈ rowPairs A r, ((oGet A p.2 : ℤ) - oGet A p.1))
      ≤ ∑ r ∈ Finset.range ℓ, ((r + 1) * N : ℤ) := Finset.sum_le_sum h_each
    _ = (ℓ * (ℓ + 1) / 2) * N := by
        -- Σ_{r=0}^{ℓ-1} (r+1) = ℓ(ℓ+1)/2
        have h_div : (2 : ℤ) ∣ (ℓ : ℤ) * (ℓ + 1) := by
          rcases Nat.even_or_odd ℓ with ⟨m, hm⟩ | ⟨m, hm⟩
          · exact ⟨m * (ℓ + 1), by push_cast [hm]; ring⟩
          · refine ⟨ℓ * (m+1), ?_⟩
            have : (ℓ : ℤ) + 1 = 2 * (m + 1) := by push_cast [hm]; ring
            rw [this]; ring
        have h_sum : ∑ r ∈ Finset.range ℓ, ((r + 1 : ℕ) : ℤ) = (ℓ : ℤ) * (ℓ + 1) / 2 := by
          have hS : (∑ r ∈ Finset.range ℓ, (r : ℤ)) * 2 = (ℓ : ℤ) * (ℓ - 1) := by
            have := Finset.sum_range_id_mul_two ℓ
            zify at this
            rcases Nat.eq_zero_or_pos ℓ with hℓ0 | hℓpos
            · subst hℓ0; simp at this ⊢
            · have h_lift : ((ℓ - 1 : ℕ) : ℤ) = (ℓ : ℤ) - 1 := by omega
              rw [h_lift] at this
              exact this
          have h_simp : ∑ r ∈ Finset.range ℓ, ((r + 1 : ℕ) : ℤ) =
              (∑ r ∈ Finset.range ℓ, (r : ℤ)) + ℓ := by
            rw [show (fun r : ℕ => ((r + 1 : ℕ) : ℤ)) = (fun r : ℕ => (r : ℤ) + 1) from
              funext (fun r => by push_cast; ring)]
            rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
            ring
          rw [h_simp]
          obtain ⟨q, hq⟩ := h_div
          have h_div2 : (2 : ℤ) ∣ (ℓ : ℤ) * ((ℓ : ℤ) - 1) := by
            rcases Nat.even_or_odd ℓ with ⟨m, hm⟩ | ⟨m, hm⟩
            · exact ⟨m * ((ℓ : ℤ) - 1), by push_cast [hm]; ring⟩
            · refine ⟨ℓ * m, ?_⟩
              have : (ℓ : ℤ) - 1 = 2 * m := by push_cast [hm]; ring
              rw [this]; ring
          obtain ⟨q', hq'⟩ := h_div2
          rw [hq, Int.mul_ediv_cancel_left _ (by norm_num : (2 : ℤ) ≠ 0)]
          rw [hq'] at hS
          have hS' : (∑ r ∈ Finset.range ℓ, (r : ℤ)) = q' := by linarith
          linarith [hS', hq, hq']
        rw [show (∑ r ∈ Finset.range ℓ, ((r + 1) * N : ℤ)) =
            (∑ r ∈ Finset.range ℓ, ((r + 1 : ℕ) : ℤ)) * N from by
          rw [Finset.sum_mul]
          apply Finset.sum_congr rfl
          intro r _; push_cast; ring]
        rw [h_sum]

/-- The sum of orderDiffs equals the sum over orderPairs of (oGet j - oGet i), in ℤ. -/
private lemma sum_orderDiffs_eq_sum_orderPairs (A : Finset ℕ) (ℓ : ℕ) (hS : IsSidonSet A) :
    (((orderDiffs A ℓ).sum id : ℕ) : ℤ) =
    ∑ p ∈ orderPairs A ℓ, ((oGet A p.2 : ℤ) - oGet A p.1) := by
  -- Step 1: sum over image = sum over orderPairs (using Sidon injection)
  have h_img : (orderDiffs A ℓ).sum id =
      (orderPairs A ℓ).sum (fun p => oGet A p.2 - oGet A p.1) := by
    unfold orderDiffs
    rw [Finset.sum_image (fun p1 hp1 p2 hp2 heq =>
      orderDiffs_injOn A ℓ hS hp1 hp2 heq)]
    rfl
  rw [h_img]
  -- Step 2: cast each ℕ-difference to ℤ-difference (uses oGet_lt to ensure non-negativity)
  push_cast
  apply Finset.sum_congr rfl
  intro p hp
  rw [mem_orderPairs_iff] at hp
  obtain ⟨hi, hj, hij, _⟩ := hp
  have h_lt : oGet A p.1 < oGet A p.2 := oGet_lt hi hj hij
  omega

/-- **Order-bounded difference counting for Sidon sets.**

    For Sidon A = {a₀ < ... < a_{k-1}} ⊆ {0,...,N}, the differences of
    orders 1 through ℓ (i.e., a_{i+r} - a_i for 1 ≤ r ≤ ℓ) form a set D of
    m distinct positive naturals where 2m = ℓ(2k-ℓ-1) and sum(D) ≤ ℓ(ℓ+1)N/2.

    **Proof:** Construct D := orderDiffs A ℓ (image of (i, j) ↦ a_j - a_i for
    pairs i < j ≤ i + ℓ). Sidon distinct-differences gives injectivity. The
    sum bound uses telescoping on each row r ∈ {0, …, ℓ-1}: the row sum
    Σ_i (a_{i+r+1} - a_i) telescopes to ≤ (r+1)·N via single-step differences,
    yielding total ≤ N · ℓ(ℓ+1)/2.

    **Reference:** Lindström (1969); BFR (2023) Section 2.
    Stated in doubled form (2m, 2·sum) to avoid ℕ division. -/
theorem order_diff_counting (A : Finset ℕ) (N : ℕ) (ℓ : ℕ)
    (hS : IsSidonSet A) (hA : A ⊆ Finset.range (N + 1))
    (hℓ : ℓ < A.card) (_hℓ_pos : 0 < ℓ) :
    ∃ D : Finset ℕ,
      D.card * 2 = ℓ * (2 * A.card - ℓ - 1) ∧
      (∀ d ∈ D, 0 < d) ∧
      2 * D.sum id ≤ ℓ * (ℓ + 1) * N := by
  refine ⟨orderDiffs A ℓ, ?_, orderDiffs_pos A ℓ, ?_⟩
  · rw [orderDiffs_card_eq A ℓ hS]
    exact orderPairs_card_doubled A ℓ (le_of_lt hℓ)
  · -- 2 * D.sum id ≤ ℓ * (ℓ + 1) * N
    -- Let S = (orderDiffs A ℓ).sum id : ℕ. Show 2 * S ≤ ℓ * (ℓ + 1) * N.
    set S : ℕ := (orderDiffs A ℓ).sum id with hS_def
    have h_sum_int : ((S : ℕ) : ℤ) =
        ∑ p ∈ orderPairs A ℓ, ((oGet A p.2 : ℤ) - oGet A p.1) := by
      rw [hS_def]
      exact sum_orderDiffs_eq_sum_orderPairs A ℓ hS
    have h_bound_int : ((S : ℕ) : ℤ) ≤ (ℓ * (ℓ + 1) / 2) * N := by
      rw [h_sum_int]
      exact orderPairs_sum_bound A N ℓ hA
    -- 2 * sum ≤ ℓ * (ℓ + 1) * N
    have h_div : (2 : ℤ) ∣ (ℓ : ℤ) * (ℓ + 1) := by
      rcases Nat.even_or_odd ℓ with ⟨m, hm⟩ | ⟨m, hm⟩
      · exact ⟨m * (ℓ + 1), by push_cast [hm]; ring⟩
      · refine ⟨ℓ * (m+1), ?_⟩
        have : (ℓ : ℤ) + 1 = 2 * (m + 1) := by push_cast [hm]; ring
        rw [this]; ring
    obtain ⟨q, hq⟩ := h_div
    have h_eq : (ℓ * (ℓ + 1) / 2 : ℤ) = q := by
      rw [hq, Int.mul_ediv_cancel_left _ (by norm_num : (2 : ℤ) ≠ 0)]
    rw [h_eq] at h_bound_int
    -- h_bound_int : (S : ℤ) ≤ q * N. From hq: 2q = ℓ(ℓ+1).
    have h_int : ((2 * S : ℕ) : ℤ) ≤ ((ℓ * (ℓ + 1) * N : ℕ) : ℤ) := by
      push_cast
      have h2q : (2 : ℤ) * q = (ℓ : ℤ) * (ℓ + 1) := by linarith
      calc (2 : ℤ) * (S : ℤ) ≤ 2 * (q * N) := by linarith
        _ = 2 * q * N := by ring
        _ = (ℓ : ℤ) * (ℓ + 1) * N := by rw [h2q]
    exact_mod_cast h_int

/-! ### Lindström core quadratic inequality -/

/-- **Core Lindström quadratic inequality (order-of-differences).**

    For Sidon A ⊆ {0,...,N} and parameter ℓ with ℓ < k:
      ℓ · (2k - ℓ - 1)² ≤ 4(ℓ+1) · N.

    **Statement note:** Multiplied through by 4 to avoid ℕ division in (ℓ+1)/2.
    The classical real-arithmetic form ℓ(k-(ℓ+1)/2)² ≤ (ℓ+1)N is equivalent
    via 2m = ℓ(2k-ℓ-1) where m = ∑_{r=1}^ℓ (k-r).

    **Proof:** From `order_diff_counting` + `sum_distinct_pos_ge` + algebra.
    Chain: m(m+1) ≤ 2·sum(D) ≤ ℓ(ℓ+1)N, then (2m)² ≤ 4m(m+1) ≤ 4ℓ(ℓ+1)N,
    and 2m = ℓ(2k-ℓ-1), so ℓ²(2k-ℓ-1)² ≤ 4ℓ(ℓ+1)N. Cancel ℓ.

    **Reference:** Lindström (1969), reproduced in BFR (2023) Section 2. -/
theorem lindstrom_quadratic (A : Finset ℕ) (N : ℕ) (ℓ : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (hℓ : ℓ < A.card) (hℓ_pos : 0 < ℓ) :
    ℓ * (2 * A.card - ℓ - 1) ^ 2 ≤ 4 * (ℓ + 1) * N := by
  obtain ⟨D, hcount, hpos, hsum⟩ := order_diff_counting A N ℓ hS hA hℓ hℓ_pos
  have h_lb := sum_distinct_pos_ge D hpos
  -- h_lb : D.card * (D.card + 1) ≤ 2 * D.sum id
  -- hcount : D.card * 2 = ℓ * (2 * A.card - ℓ - 1)
  -- hsum : 2 * D.sum id ≤ ℓ * (ℓ + 1) * N
  -- Chain: D.card * (D.card + 1) ≤ ℓ * (ℓ + 1) * N
  have hmm : D.card * (D.card + 1) ≤ ℓ * (ℓ + 1) * N := by linarith
  -- D.card² ≤ D.card*(D.card+1) ≤ ℓ*(ℓ+1)*N
  have h_sq : D.card * D.card ≤ ℓ * (ℓ + 1) * N := by nlinarith
  -- Multiply suffices by ℓ, then cancel:
  -- ℓ*(ℓ*q²) = (ℓ*q)² = (2*D.card)² = 4*D.card² ≤ 4*ℓ*(ℓ+1)*N = ℓ*(4*(ℓ+1)*N)
  suffices hsuff : ℓ * (ℓ * (2 * A.card - ℓ - 1) ^ 2) ≤ ℓ * (4 * (ℓ + 1) * N) by
    exact Nat.le_of_mul_le_mul_left hsuff hℓ_pos
  -- Rewrite ^2 to * for ring reasoning
  set q := 2 * A.card - ℓ - 1 with hq_def
  -- hcount : D.card * 2 = ℓ * q
  -- Goal: ℓ * (ℓ * q ^ 2) ≤ ℓ * (4 * (ℓ + 1) * N)
  -- LHS = (ℓ*q)*(ℓ*q) = (D.card*2)*(D.card*2) = 4*D.card*D.card
  -- RHS = 4*ℓ*(ℓ+1)*N ≥ 4*D.card*D.card (from h_sq)
  -- First show LHS = (D.card*2)*(D.card*2):
  have h_lhs : ℓ * (ℓ * q ^ 2) = D.card * 2 * (D.card * 2) := by
    have : q ^ 2 = q * q := sq q
    rw [this]
    nlinarith [hcount]
  rw [h_lhs]
  -- Goal: D.card * 2 * (D.card * 2) ≤ ℓ * (4 * (ℓ + 1) * N)
  nlinarith [h_sq]

/-- **Weak Lindström bound: k ≤ √(2N) + 1.**
    From lindstrom_quadratic with ℓ = 1. This is provable and gives a clean
    ℕ statement, though weaker than the full √N + N^{1/4} + 1. -/
theorem lindstrom_bound_weak (A : Finset ℕ) (N : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (hN : 0 < N)
    (hk : 1 < A.card) :
    A.card ≤ Nat.sqrt (2 * N) + 1 := by
  have hq := lindstrom_quadratic A N 1 hS hA hk (by omega)
  simp only [one_mul] at hq
  -- hq : (2 * A.card - 1 - 1) ^ 2 ≤ 4 * 2 * N
  set k := A.card
  -- Step 1: (k-1)*(k-1) ≤ 2*N from (2k-2)^2 ≤ 8N
  have h2k : 2 * k - 1 - 1 = 2 * (k - 1) := by omega
  rw [h2k] at hq
  -- hq : (2*(k-1))^2 ≤ 8*N
  have h_sq : (k - 1) * (k - 1) ≤ 2 * N := by nlinarith [hq, sq_nonneg (k - 1)]
  -- Step 2: (k-1) ≤ Nat.sqrt(2*N) via Nat.le_sqrt
  have h_le : k - 1 ≤ Nat.sqrt (2 * N) := Nat.le_sqrt.mpr h_sq
  omega

/-- **Lindström bound: k ≤ √N + ⁴√N + 1.**
    From lindstrom_quadratic with ℓ = Nat.sqrt (Nat.sqrt N).

    **⚠ Statement requires real-number intermediate step.**
    The quadratic `ℓ(2k-ℓ-1)² ≤ 4(ℓ+1)N` yields (over ℝ):
      k ≤ √(N(ℓ+1)/ℓ) + (ℓ+1)/2
    Converting to ℕ with ℓ = ⌊⁴√N⌋ requires bounding:
      Nat.sqrt(N + N/ℓ) ≤ Nat.sqrt N + corrections
    which needs either Taylor-expansion-style ℕ bounds on Nat.sqrt
    or a direct proof that Nat.sqrt rounding doesn't accumulate.

    `lindstrom_bound_weak` (k ≤ √(2N)+1) is the strongest version
    provable from `lindstrom_quadratic` via pure ℕ arithmetic.

    **Reference:** Lindström (1969), J. Combinatorial Theory 7(1). -/
axiom lindstrom_bound (A : Finset ℕ) (N : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (hN : 0 < N) :
    A.card ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 1

/-! ## Summary

  **Proven (zero sorries):**
  1. pigeonhole_residue: partition + pigeonhole
  2. residue_class_scaled_sidon: Sidon preservation under modular scaling
  3. scaled_range: range containment of scaled residue class
  4. lindstrom_parametric: k² ≤ 2Nt + kt (residue class parametric bound)
  5. sidon_distinct_differences: all positive differences are distinct
  6. card_posDiffs_sidon: |posDiffs A| = k(k-1)/2
  7. posDiffs_subset_Icc: differences ≤ N
  8. sum_distinct_pos_ge: sum of m distinct positive naturals ≥ m(m+1)/2
  9. lindstrom_quadratic: ℓ(2k-ℓ-1)² ≤ 4(ℓ+1)N (from axiom + algebra)
  10. lindstrom_bound_weak: k ≤ √(2N) + 1 (from 9 with ℓ=1 + Nat.le_sqrt)

  **Axioms (3 total, all with references):**
  - sidon_elem_bound: k(k-1) ≤ 2M (proved in Erdos30_Complete, pending Mathlib API port)
  - order_diff_counting: sorted enumeration + telescoping (requires orderEmbOfFin)
  - lindstrom_bound: k ≤ ⌊√N⌋ + ⌊⁴√N⌋ + 1 (requires ℝ→ℕ Nat.sqrt bounding)

  **Proof dependency DAG:**
    pigeonhole_residue ──┐
    residue_class_scaled_sidon ──┤
    scaled_range ──┤──→ lindstrom_parametric (for BFR)
    sidon_elem_bound ──┘

    order_diff_counting ──┐
    sum_distinct_pos_ge ──┤──→ lindstrom_quadratic ──→ lindstrom_bound_weak ✓
                          │                        └──→ lindstrom_bound (axiom: ℝ→ℕ step)

  **Note on proof architecture:**
  The residue class machinery (pigeonhole, scaling, parametric) is needed
  for the BFR 0.998 improvement (Erdos30_BFR.lean), not for the basic
  Lindström bound. The Lindström bound uses a different technique:
  counting differences by order in the sorted sequence.

  **Statement correction (prior session):** The original ℕ formulation
  ℓ·(k-(ℓ+1)/2)² ≤ (ℓ+1)·N is false for even ℓ with tight parameters
  (e.g. ℓ=4, k=10, N=47). Fixed to ℓ·(2k-ℓ-1)² ≤ 4·(ℓ+1)·N which
  avoids ℕ division entirely by multiplying through by 4.
-/

end Erdos.Sidon

