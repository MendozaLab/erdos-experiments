/-
  Erdős #30 — Lindström Quadratic Inequality (Assembly v2)
  =========================================================

  Combines all proven B2 components into the Lindström quadratic:
    ℓ * (2k - ℓ - 1)² ≤ 4 * (ℓ + 1) * N

  Proven theorems (inlined, compiled):
  - sorted_enum_exists: Finset.sort construction (PROVEN — was axiom)
  - sidon_diff_injective: 4-way case split
  - order_diff_inj: uses sidon_diff_injective
  - order_diff_card_formula: induction + division avoidance
  - per_order_bound: nested telescoping
  - total_bound: aggregation over per_order_bound
  - order_diff_lower_bound: distinct differences + positive sum (NEW — was axiom)
  - lindstrom_quadratic: algebraic chain from axioms to bound

  Axiom inventory (1 remaining):
  - lindstrom_bound_extraction: ℝ→ℕ sqrt extraction

  Score: 9 proven theorems, 1 axiom (down from 2).

  2026-04-12 — Lean 4 v4.24.0, Mathlib f897ebcf
-/

import Mathlib.Data.Finset.Sort
import Mathlib.Data.Finset.Prod
import Mathlib.Data.Finset.Image
import Mathlib.Data.Finset.Interval
import Mathlib.Order.Interval.Finset.Nat
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Push

open Finset Nat

namespace Erdos.Lindstrom

/-! ═══════════════════════════════════════════════════════════
    DEFINITIONS
    ═══════════════════════════════════════════════════════════ -/

abbrev IsSidonSet (A : Finset ℕ) : Prop :=
  ∀ a₁ ∈ A, ∀ b₁ ∈ A, ∀ a₂ ∈ A, ∀ b₂ ∈ A,
    a₁ ≤ b₁ → a₂ ≤ b₂ → a₁ + b₁ = a₂ + b₂ → (a₁ = a₂ ∧ b₁ = b₂)

/-! ═══════════════════════════════════════════════════════════
    SECTION A: DIFFERENCE INJECTIVITY (proven)
    ═══════════════════════════════════════════════════════════ -/

theorem sidon_diff_injective (A : Finset ℕ) (hS : IsSidonSet A)
    (a₁ b₁ a₂ b₂ : ℕ)
    (ha₁ : a₁ ∈ A) (hb₁ : b₁ ∈ A) (ha₂ : a₂ ∈ A) (hb₂ : b₂ ∈ A)
    (hlt₁ : b₁ < a₁) (hlt₂ : b₂ < a₂)
    (heq : a₁ - b₁ = a₂ - b₂) :
    a₁ = a₂ ∧ b₁ = b₂ := by
  have h_sum : a₁ + b₂ = a₂ + b₁ := by omega
  by_cases h1 : b₂ ≤ a₁
  · by_cases h2 : b₁ ≤ a₂
    · obtain ⟨heq1, heq2⟩ := hS b₂ hb₂ a₁ ha₁ b₁ hb₁ a₂ ha₂ h1 h2 (by omega)
      exact ⟨heq2, heq1.symm⟩
    · push_neg at h2
      obtain ⟨heq1, heq2⟩ := hS b₂ hb₂ a₁ ha₁ a₂ ha₂ b₁ hb₁ h1 (le_of_lt h2) (by omega)
      omega
  · push_neg at h1
    by_cases h2 : b₁ ≤ a₂
    · obtain ⟨heq1, heq2⟩ := hS a₁ ha₁ b₂ hb₂ b₁ hb₁ a₂ ha₂ (le_of_lt h1) h2 (by omega)
      omega
    · push_neg at h2
      obtain ⟨heq1, heq2⟩ := hS a₁ ha₁ b₂ hb₂ a₂ ha₂ b₁ hb₁ (le_of_lt h1) (le_of_lt h2) (by omega)
      exact ⟨heq1, heq2.symm⟩

theorem order_diff_inj (A : Finset ℕ) (hS : IsSidonSet A) (k ℓ : ℕ)
    (r r' i i' : ℕ)
    (hr : 1 ≤ r ∧ r ≤ ℓ) (hr' : 1 ≤ r' ∧ r' ≤ ℓ)
    (hi : i + r < k) (hi' : i' + r' < k)
    (a : ℕ → ℕ) (ha_mem : ∀ j, j < k → a j ∈ A)
    (ha_strict_mono : ∀ j₁ j₂ : ℕ, j₁ < j₂ → j₂ < k → a j₁ < a j₂)
    (heq_diff : a (i + r) - a i = a (i' + r') - a i') :
    r = r' ∧ i = i' := by
  have h_i_lt_k : i < k := by omega
  have h_i'_lt_k : i' < k := by omega
  obtain ⟨h_eq_a, h_eq_b⟩ := sidon_diff_injective A hS
    (a (i + r)) (a i) (a (i' + r')) (a i')
    (ha_mem (i + r) hi) (ha_mem i h_i_lt_k)
    (ha_mem (i' + r') hi') (ha_mem i' h_i'_lt_k)
    (ha_strict_mono i (i + r) (by omega) hi)
    (ha_strict_mono i' (i' + r') (by omega) hi')
    heq_diff
  have h_i_eq : i = i' := by
    by_contra h_ne
    rcases Nat.lt_or_gt_of_ne h_ne with h_lt | h_gt
    · exact absurd h_eq_b (_root_.ne_of_lt (ha_strict_mono i i' h_lt h_i'_lt_k))
    · exact absurd h_eq_b.symm (_root_.ne_of_lt (ha_strict_mono i' i h_gt h_i_lt_k))
  have h_ir_eq : i + r = i' + r' := by
    by_contra h_ne
    rcases Nat.lt_or_gt_of_ne h_ne with h_lt | h_gt
    · exact absurd h_eq_a (_root_.ne_of_lt (ha_strict_mono (i + r) (i' + r') h_lt hi'))
    · exact absurd h_eq_a.symm (_root_.ne_of_lt (ha_strict_mono (i' + r') (i + r) h_gt hi))
  exact ⟨by omega, h_i_eq⟩

/-! ═══════════════════════════════════════════════════════════
    SECTION B: CARDINALITY FORMULA (proven)
    ═══════════════════════════════════════════════════════════ -/

lemma doubled_sum_aux (k ℓ : ℕ) (hℓ : ℓ ≤ k) :
    2 * (∑ r ∈ Finset.range ℓ, (k - r - 1)) + ℓ * (ℓ + 1) = 2 * ℓ * k := by
  induction ℓ with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ]
    have h_n_le_k : n ≤ k := by omega
    have ih' := ih h_n_le_k
    have h1 : (n + 1) * (n + 1 + 1) = n * (n + 1) + 2 * (n + 1) := by ring
    have h2 : 2 * (n + 1) * k = 2 * n * k + 2 * k := by ring
    omega

theorem order_diff_card_formula (k ℓ : ℕ) (hℓ : ℓ < k) :
    ∃ m : ℕ, m * 2 = ℓ * (2 * k - ℓ - 1) := by
  use ∑ r ∈ Finset.range ℓ, (k - r - 1)
  have h_aux := doubled_sum_aux k ℓ (le_of_lt hℓ)
  have h_cancel : (2 * k - ℓ - 1) + (ℓ + 1) = 2 * k := by omega
  have h_rhs : ℓ * (2 * k - ℓ - 1) + ℓ * (ℓ + 1) = 2 * ℓ * k := by
    rw [← mul_add, h_cancel]; ring
  linarith

/-! ═══════════════════════════════════════════════════════════
    SECTION C: TELESCOPING BOUNDS (proven — 0 axioms)
    Copied exactly from compiled B2_telescope_full.lean
    ═══════════════════════════════════════════════════════════ -/

lemma mono_step_to_le (a : ℕ → ℕ) (j n : ℕ)
    (h_step : ∀ m, m < n → a (j + m) ≤ a (j + m + 1)) :
    a j ≤ a (j + n) := by
  induction n with
  | zero => simp
  | succ n' ih =>
    have h_prev : ∀ m, m < n' → a (j + m) ≤ a (j + m + 1) :=
      fun m hm => h_step m (by omega)
    have h1 := ih h_prev
    have h2 := h_step n' (by omega)
    show a j ≤ a (j + n' + 1)
    linarith

lemma consecutive_telescope (a : ℕ → ℕ) (j n : ℕ)
    (h_step : ∀ m, m < n → a (j + m) ≤ a (j + m + 1)) :
    (∑ m ∈ range n, (a (j + m + 1) - a (j + m))) = a (j + n) - a j := by
  induction n with
  | zero => simp
  | succ n' ih =>
    rw [sum_range_succ]
    have h_prev : ∀ m, m < n' → a (j + m) ≤ a (j + m + 1) :=
      fun m hm => h_step m (by omega)
    rw [ih h_prev]
    have h_j_le : a j ≤ a (j + n') := mono_step_to_le a j n' h_prev
    have h_n_le : a (j + n') ≤ a (j + n' + 1) := h_step n' (by omega)
    have h_trans : a j ≤ a (j + n' + 1) := le_trans h_j_le h_n_le
    have eq1 := Nat.sub_add_cancel h_j_le
    have eq2 := Nat.sub_add_cancel h_n_le
    have eq3 := Nat.sub_add_cancel h_trans
    show a (j + n') - a j + (a (j + n' + 1) - a (j + n')) = a (j + n' + 1) - a j
    linarith

lemma difference_decomposition (a : ℕ → ℕ) (i r : ℕ)
    (h_step : ∀ m, m < r → a (i + m) ≤ a (i + m + 1)) :
    a (i + r) - a i = ∑ j ∈ range r, (a (i + j + 1) - a (i + j)) := by
  exact (consecutive_telescope a i r h_step).symm

lemma inner_sum_telescopes (a : ℕ → ℕ) (k r j : ℕ)
    (hj : j < r) (hrk : r < k)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    ∑ i ∈ range (k - r), (a (i + j + 1) - a (i + j))
    = a (k - r + j) - a j := by
  have h_reindex : ∀ i ∈ range (k - r), a (i + j + 1) - a (i + j)
                 = a (j + i + 1) - a (j + i) := by
    intro i _
    have h1 : i + j + 1 = j + i + 1 := by omega
    have h2 : i + j = j + i := by omega
    rw [h1, h2]
  rw [sum_congr rfl h_reindex]
  have h_step : ∀ m, m < (k - r) → a (j + m) ≤ a (j + m + 1) := by
    intro m hm
    exact le_of_lt (ha_mono (j + m) (j + m + 1) (by omega) (by omega))
  rw [consecutive_telescope a j (k - r) h_step]
  have h_eq : j + (k - r) = k - r + j := by omega
  rw [h_eq]

lemma bounded_telescoped_sum (a : ℕ → ℕ) (k r j N : ℕ)
    (hj : j < r) (hrk : r < k)
    (ha_bound : ∀ x, x < k → a x ≤ N) :
    a (k - r + j) - a j ≤ N := by
  have h_bound_endpoint : k - r + j < k := by omega
  have h_a_end_bound := ha_bound (k - r + j) h_bound_endpoint
  omega

theorem per_order_bound (a : ℕ → ℕ) (k r N : ℕ)
    (hr : 1 ≤ r) (hrk : r < k)
    (ha_bound : ∀ j, j < k → a j ≤ N)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    ∑ i ∈ range (k - r), (a (i + r) - a i) ≤ r * N := by
  have h_decomp : ∀ i ∈ range (k - r),
      a (i + r) - a i = ∑ j ∈ range r, (a (i + j + 1) - a (i + j)) := by
    intro i hi
    simp only [mem_range] at hi
    have h_step : ∀ m, m < r → a (i + m) ≤ a (i + m + 1) := by
      intro m hm
      exact le_of_lt (ha_mono (i + m) (i + m + 1) (by omega) (by omega))
    exact difference_decomposition a i r h_step
  rw [sum_congr rfl h_decomp, Finset.sum_comm]
  have h_each : ∀ j ∈ range r,
      ∑ i ∈ range (k - r), (a (i + j + 1) - a (i + j)) ≤ N := by
    intro j hj
    simp only [mem_range] at hj
    rw [inner_sum_telescopes a k r j hj hrk ha_mono]
    exact bounded_telescoped_sum a k r j N hj hrk ha_bound
  calc ∑ j ∈ range r, ∑ i ∈ range (k - r), (a (i + j + 1) - a (i + j))
      ≤ ∑ j ∈ range r, N := sum_le_sum h_each
    _ = r * N := by rw [sum_const, card_range]; ring

lemma doubled_arith_series (ℓ : ℕ) :
    2 * ∑ r ∈ range ℓ, (r + 1) = ℓ * (ℓ + 1) := by
  induction ℓ with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ]
    have h1 : (n + 1) * (n + 1 + 1) = n * (n + 1) + 2 * (n + 1) := by ring
    omega

lemma sum_mul_factor (ℓ N : ℕ) :
    ∑ r ∈ range ℓ, ((r + 1) * N) = (∑ r ∈ range ℓ, (r + 1)) * N := by
  induction ℓ with
  | zero => simp
  | succ n ih =>
    rw [sum_range_succ, sum_range_succ, ih]
    ring

theorem total_bound (a : ℕ → ℕ) (k ℓ N : ℕ)
    (hℓ : 1 ≤ ℓ) (hℓk : ℓ < k)
    (ha_bound : ∀ j, j < k → a j ≤ N)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    2 * ∑ r ∈ range ℓ, ∑ i ∈ range (k - (r + 1)), (a (i + (r + 1)) - a i)
    ≤ ℓ * (ℓ + 1) * N := by
  have h_per : ∀ r ∈ range ℓ,
      ∑ i ∈ range (k - (r + 1)), (a (i + (r + 1)) - a i) ≤ (r + 1) * N := by
    intro r hr
    simp only [mem_range] at hr
    exact per_order_bound a k (r + 1) N (by omega) (by omega) ha_bound ha_mono
  have h_sum : ∑ r ∈ range ℓ, ∑ i ∈ range (k - (r + 1)), (a (i + (r + 1)) - a i)
              ≤ ∑ r ∈ range ℓ, ((r + 1) * N) := Finset.sum_le_sum h_per
  rw [sum_mul_factor] at h_sum
  have h_series := doubled_arith_series ℓ
  nlinarith

/-! ═══════════════════════════════════════════════════════════
    SECTION D: AXIOM 2 HELPER DEFINITIONS & LEMMAS (NEW)
    From axiom2_full.lean — distinct difference sum bound
    ═══════════════════════════════════════════════════════════ -/

theorem distinct_positive_sum_bound : ∀ (n : ℕ) (S : Finset ℕ),
    S.card = n → (∀ x ∈ S, 0 < x) → n * (n + 1) ≤ 2 * ∑ x ∈ S, x := by
  intro n
  induction n with
  | zero => intro S hcard _; simp at hcard; simp [hcard]
  | succ m ih =>
    intro S hcard h_pos
    have hne : S.Nonempty := Finset.card_pos.mp (by omega)
    set a := S.min' hne with ha_def
    have ha_mem : a ∈ S := Finset.min'_mem S hne
    have ha_pos : 0 < a := h_pos a ha_mem
    have ha_min : ∀ x ∈ S, a ≤ x := fun x hx => Finset.min'_le S x hx
    set S' := S.erase a with hS'_def
    have hS'_card : S'.card = m := by
      rw [hS'_def, Finset.card_erase_of_mem ha_mem, hcard]; omega
    have hS'_ge2 : ∀ x ∈ S', 2 ≤ x := by
      intro x hx
      have h1 : x ∈ S := Finset.mem_of_mem_erase hx
      have h2 : x ≠ a := Finset.ne_of_mem_erase hx
      have := ha_min x h1; omega
    have h_sum_split : ∑ x ∈ S, x = a + ∑ x ∈ S', x := by
      rw [hS'_def, ← Finset.add_sum_erase S (fun x => x) ha_mem]
    have h_inj : ∀ x ∈ S', ∀ y ∈ S', (x - 1 : ℕ) = y - 1 → x = y := by
      intro x hx y hy hxy; have := hS'_ge2 x hx; have := hS'_ge2 y hy; omega
    set T := S'.image (· - 1 : ℕ → ℕ) with hT_def
    have hT_card : T.card = m := by
      rw [hT_def, Finset.card_image_of_injOn]; · exact hS'_card
      · intro x hx y hy hxy
        exact h_inj x (Finset.mem_coe.mp hx) y (Finset.mem_coe.mp hy) hxy
    have hT_pos : ∀ y ∈ T, 0 < y := by
      intro y hy; rw [hT_def, Finset.mem_image] at hy
      obtain ⟨x, hx, rfl⟩ := hy; have := hS'_ge2 x hx; omega
    have ih_T := ih T hT_card hT_pos
    have h_sum_T : ∑ y ∈ T, y = ∑ x ∈ S', (x - 1) := by
      rw [hT_def]; exact Finset.sum_image h_inj
    have h_shift_sum : (∑ x ∈ S', (x - 1)) + S'.card = ∑ x ∈ S', x := by
      have h_congr : ∀ x ∈ S', (x - 1) + 1 = x := by
        intro x hx; have := hS'_ge2 x hx; omega
      have h_split : ∑ x ∈ S', ((x - 1) + 1) = ∑ x ∈ S', (x - 1) + ∑ x ∈ S', 1 :=
        Finset.sum_add_distrib
      have h_eq : ∑ x ∈ S', ((x - 1) + 1) = ∑ x ∈ S', x :=
        Finset.sum_congr rfl h_congr
      have h_ones : ∑ _x ∈ S', (1 : ℕ) = S'.card := by simp
      linarith
    rw [h_sum_T] at ih_T
    have h_S'_sum_ge : m * (m + 1) + 2 * m ≤ 2 * ∑ x ∈ S', x := by
      have h_eq : (∑ x ∈ S', (x - 1)) + m = ∑ x ∈ S', x := by
        rw [← hS'_card]; exact h_shift_sum
      omega
    calc (m + 1) * (m + 1 + 1)
        = m * (m + 1) + 2 * m + 2 := by ring
      _ ≤ 2 * ∑ x ∈ S', x + 2 := by linarith [h_S'_sum_ge]
      _ ≤ 2 * ∑ x ∈ S', x + 2 * a := by linarith [ha_pos]
      _ = 2 * (a + ∑ x ∈ S', x) := by ring
      _ = 2 * ∑ x ∈ S, x := by rw [h_sum_split]

/-- The index Finset for order-differences. -/
def diffPairSet (k ℓ : ℕ) : Finset (Σ _ : ℕ, ℕ) :=
  (range ℓ).sigma (fun r => range (k - (r + 1)))

/-- The difference function. -/
def diffFn (a : ℕ → ℕ) : (Σ _ : ℕ, ℕ) → ℕ :=
  fun p => a (p.2 + (p.1 + 1)) - a p.2

/-- Double sum = sum over sigma. -/
lemma double_sum_eq_sigma_sum (a : ℕ → ℕ) (k ℓ : ℕ) :
    ∑ r ∈ range ℓ, ∑ i ∈ range (k - (r + 1)), (a (i + (r + 1)) - a i)
    = ∑ p ∈ diffPairSet k ℓ, diffFn a p := by
  unfold diffPairSet diffFn
  rw [Finset.sum_sigma]

/-- Helper: 2 * ∑_{r<ℓ} (k-r-1) + ℓ*(ℓ+1) = 2*ℓ*k. -/
lemma doubled_sum_aux_card (k ℓ : ℕ) (hℓk : ℓ ≤ k) :
    2 * (∑ r ∈ Finset.range ℓ, (k - r - 1)) + ℓ * (ℓ + 1) = 2 * ℓ * k := by
  induction ℓ with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ]
    have h_n_le_k : n ≤ k := by omega
    have ih' := ih h_n_le_k
    have h1 : (n + 1) * (n + 1 + 1) = n * (n + 1) + 2 * (n + 1) := by ring
    have h2 : 2 * (n + 1) * k = 2 * n * k + 2 * k := by ring
    omega

/-- Cardinality of index set = m (when m*2 = ℓ*(2k-ℓ-1)). -/
lemma diffPairSet_card_eq (k ℓ m : ℕ) (hm : m * 2 = ℓ * (2 * k - ℓ - 1)) (hℓk : ℓ < k) :
    (diffPairSet k ℓ).card = m := by
  unfold diffPairSet
  rw [Finset.card_sigma]
  simp only [Finset.card_range]
  -- Goal: ∑ r ∈ range ℓ, (k - (r + 1)) = m
  -- Step 1: rewrite k-(r+1) = k-r-1
  have h_congr : ∀ r ∈ range ℓ, k - (r + 1) = k - r - 1 := by
    intro r hr; simp only [Finset.mem_range] at hr; omega
  rw [Finset.sum_congr rfl h_congr]
  -- Step 2: from doubled_sum_aux_card, 2*∑(k-r-1) + ℓ*(ℓ+1) = 2*ℓ*k
  have h_doubled := doubled_sum_aux_card k ℓ (le_of_lt hℓk)
  -- Factor: ℓ*(2k-ℓ-1) + ℓ*(ℓ+1) = ℓ*((2k-ℓ-1)+(ℓ+1)) = ℓ*2k = 2ℓk
  have h_sub_ok : ℓ + 1 ≤ 2 * k := by omega
  have h_add_back : (2 * k - ℓ - 1) + (ℓ + 1) = 2 * k := by omega
  have h_factor : ℓ * (2 * k - ℓ - 1) + ℓ * (ℓ + 1) = 2 * ℓ * k := by
    rw [← mul_add, h_add_back]; ring
  -- From h_doubled and h_factor: 2*∑ + ℓ*(ℓ+1) = ℓ*(2k-ℓ-1) + ℓ*(ℓ+1)
  -- So 2*∑ = ℓ*(2k-ℓ-1) = m*2, hence ∑ = m
  omega

/-- The diff function is injective on the pair set. -/
lemma diffFn_injOn (a : ℕ → ℕ) (k ℓ : ℕ) (hℓk : ℓ < k)
    (A : Finset ℕ) (hS : IsSidonSet A)
    (ha_mem : ∀ j, j < k → a j ∈ A)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    ∀ p ∈ diffPairSet k ℓ, ∀ q ∈ diffPairSet k ℓ,
      diffFn a p = diffFn a q → p = q := by
  intro p hp q hq heq
  unfold diffPairSet at hp hq
  simp at hp hq
  unfold diffFn at heq
  obtain ⟨r, i⟩ := p
  obtain ⟨r', i'⟩ := q
  simp at hp hq heq ⊢
  have h_inj := order_diff_inj A hS k ℓ
    (r + 1) (r' + 1) i i'
    ⟨by omega, by omega⟩ ⟨by omega, by omega⟩
    (by omega) (by omega)
    a ha_mem ha_mono heq
  exact ⟨by omega, h_inj.2⟩

/-- All differences are positive. -/
lemma diffFn_pos (a : ℕ → ℕ) (k ℓ : ℕ) (hℓk : ℓ < k)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    ∀ p ∈ diffPairSet k ℓ, 0 < diffFn a p := by
  intro p hp
  unfold diffPairSet at hp; simp at hp
  unfold diffFn
  have : a p.2 < a (p.2 + (p.1 + 1)) :=
    ha_mono p.2 (p.2 + (p.1 + 1)) (by omega) (by omega)
  omega

/-- Sum over image = sum over domain (injectivity). -/
lemma sum_image_eq (a : ℕ → ℕ) (k ℓ : ℕ) (hℓk : ℓ < k)
    (A : Finset ℕ) (hS : IsSidonSet A)
    (ha_mem : ∀ j, j < k → a j ∈ A)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂) :
    ∑ x ∈ (diffPairSet k ℓ).image (diffFn a), x
    = ∑ p ∈ diffPairSet k ℓ, diffFn a p := by
  rw [Finset.sum_image]
  exact diffFn_injOn a k ℓ hℓk A hS ha_mem ha_mono

/-! ═══════════════════════════════════════════════════════════
    SECTION E: ASSEMBLY AXIOMS (1 remaining)
    ═══════════════════════════════════════════════════════════ -/

/-- **Proven** (was axiom): A finite set of naturals has a sorted enumeration.
    Construction: Finset.sort + List.get with default. -/
theorem sorted_enum_exists (A : Finset ℕ) (_hS : IsSidonSet A)
    (N : ℕ) (hA : A ⊆ Finset.range (N + 1)) :
    ∃ (a : ℕ → ℕ),
      (∀ j, j < A.card → a j ∈ A) ∧
      (∀ j, j < A.card → a j ≤ N) ∧
      (∀ j₁ j₂, j₁ < j₂ → j₂ < A.card → a j₁ < a j₂) := by
  let s := A.sort (· ≤ ·)
  have hs_len : s.length = A.card := Finset.length_sort (· ≤ ·)
  let a : ℕ → ℕ := fun j => if h : j < s.length then s.get ⟨j, h⟩ else 0
  use a
  refine ⟨?_, ?_, ?_⟩
  · -- Property 1: membership
    intro j hj
    have hj' : j < s.length := by rw [hs_len]; exact hj
    simp only [a, dif_pos hj']
    have : s.get ⟨j, hj'⟩ ∈ s := List.get_mem s ⟨j, hj'⟩
    rwa [Finset.mem_sort] at this
  · -- Property 2: bounded by N
    intro j hj
    have hj' : j < s.length := by rw [hs_len]; exact hj
    simp only [a, dif_pos hj']
    have hmem : s.get ⟨j, hj'⟩ ∈ s := List.get_mem s ⟨j, hj'⟩
    rw [Finset.mem_sort] at hmem
    have := hA hmem
    rw [Finset.mem_range] at this
    omega
  · -- Property 3: strictly increasing
    intro j₁ j₂ hlt hj₂
    have hj₂' : j₂ < s.length := by rw [hs_len]; exact hj₂
    have hj₁' : j₁ < s.length := by omega
    simp only [a, dif_pos hj₁', dif_pos hj₂']
    have h_sorted := Finset.sort_sorted (· ≤ ·) A
    have h_nodup := Finset.sort_nodup (· ≤ ·) A
    have h_le : s.get ⟨j₁, hj₁'⟩ ≤ s.get ⟨j₂, hj₂'⟩ :=
      List.Sorted.get_mono h_sorted (Fin.mk_le_mk.mpr (le_of_lt hlt))
    have h_ne : s.get ⟨j₁, hj₁'⟩ ≠ s.get ⟨j₂, hj₂'⟩ := by
      intro heq
      have := (List.Nodup.get_inj_iff h_nodup).mp heq
      simp [Fin.ext_iff] at this
      omega
    omega

/-- **Proven** (was axiom): The order-differences lower bound.
    m = ℓ(2k-ℓ-1)/2 order-differences are distinct (by order_diff_inj) and
    positive. By minimum-sum of m distinct positives: m(m+1) ≤ 2 * (double sum).
    Proof: construct Finset D as image of diff function, prove injectivity from
    order_diff_inj, apply distinct_positive_sum_bound. -/
theorem order_diff_lower_bound (a : ℕ → ℕ) (k ℓ m : ℕ)
    (A : Finset ℕ) (hS : IsSidonSet A)
    (ha_mem : ∀ j, j < k → a j ∈ A)
    (ha_mono : ∀ j₁ j₂, j₁ < j₂ → j₂ < k → a j₁ < a j₂)
    (hm : m * 2 = ℓ * (2 * k - ℓ - 1))
    (hℓ : ℓ < k) :
    m * (m + 1) ≤ 2 * ∑ r ∈ range ℓ, ∑ i ∈ range (k - (r + 1)), (a (i + (r + 1)) - a i) := by
  -- The image Finset D = {a(i+r) - a(i) | (r,i) ∈ pairs}
  set D := (diffPairSet k ℓ).image (diffFn a) with hD_def
  -- D has card m (injectivity + pair set card)
  have hD_card : D.card = m := by
    rw [hD_def, Finset.card_image_of_injOn]
    · exact diffPairSet_card_eq k ℓ m hm hℓ
    · intro p hp q hq heq
      exact diffFn_injOn a k ℓ hℓ A hS ha_mem ha_mono p
        (Finset.mem_coe.mp hp) q (Finset.mem_coe.mp hq) heq
  -- All elements of D are positive
  have hD_pos : ∀ x ∈ D, 0 < x := by
    intro x hx
    rw [hD_def, Finset.mem_image] at hx
    obtain ⟨p, hp, rfl⟩ := hx
    exact diffFn_pos a k ℓ hℓ ha_mono p hp
  -- Apply core lemma
  have h_bound := distinct_positive_sum_bound m D hD_card hD_pos
  -- Connect sums
  have h_sum_eq : ∑ x ∈ D, x = ∑ p ∈ diffPairSet k ℓ, diffFn a p :=
    sum_image_eq a k ℓ hℓ A hS ha_mem ha_mono
  have h_sigma_eq := double_sum_eq_sigma_sum a k ℓ
  linarith

/-- Axiom: √N + ⁴√N + 1 extraction from the quadratic.
    Citation: Standard ℝ→ℕ floor argument with Nat.sqrt.
    Closable but requires Real imports (blocked in sandbox). -/
axiom lindstrom_bound_extraction (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A) (hA : A ⊆ Finset.range (N + 1)) (hN : 0 < N)
    (h_quad : ∀ ℓ : ℕ, 0 < ℓ → ℓ < A.card →
      ℓ * (2 * A.card - ℓ - 1) * (2 * A.card - ℓ - 1) ≤ 4 * (ℓ + 1) * N) :
    A.card ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 1

/-! ═══════════════════════════════════════════════════════════
    SECTION F: MAIN THEOREMS (proven from 1 axiom)
    ═══════════════════════════════════════════════════════════ -/

/-- **Lindström Quadratic Inequality** (proven modulo 1 axiom):
    For a Sidon set A ⊆ [0,N] with ℓ < |A|:
    ℓ * (2|A| - ℓ - 1)² ≤ 4(ℓ+1)*N

    Chain: order_diff_lower_bound gives m(m+1) ≤ 2*Σ.
    total_bound gives 2*Σ ≤ ℓ(ℓ+1)*N.
    So m(m+1) ≤ ℓ(ℓ+1)*N, hence m² ≤ ℓ(ℓ+1)*N.
    Since (2m)² = ℓ²(2k-ℓ-1)² and ℓ*(2k-ℓ-1)² = 4m²/ℓ,
    we get ℓ(2k-ℓ-1)² ≤ 4(ℓ+1)*N by cancelling ℓ from ℓ²(2k-ℓ-1)² ≤ 4ℓ(ℓ+1)*N. -/
theorem lindstrom_quadratic (A : Finset ℕ) (N ℓ : ℕ)
    (hS : IsSidonSet A) (hA : A ⊆ Finset.range (N + 1))
    (hℓ_pos : 0 < ℓ) (hℓ : ℓ < A.card) :
    ℓ * (2 * A.card - ℓ - 1) * (2 * A.card - ℓ - 1) ≤ 4 * (ℓ + 1) * N := by
  set k := A.card with hk_def
  obtain ⟨a, ha_mem, ha_bound, ha_mono⟩ := sorted_enum_exists A hS N hA
  -- Upper bound on double sum from total_bound
  have h_tb := total_bound a k ℓ N (by omega) hℓ ha_bound ha_mono
  -- Cardinality of order-differences
  obtain ⟨m, hm_eq⟩ := order_diff_card_formula k ℓ hℓ
  -- Lower bound from distinct differences
  have h_lb := order_diff_lower_bound a k ℓ m A hS ha_mem ha_mono hm_eq hℓ
  -- Chain: m(m+1) ≤ 2*Σ ≤ ℓ(ℓ+1)*N
  have h_key : m * (m + 1) ≤ ℓ * (ℓ + 1) * N := by linarith
  -- Algebraic deduction: ℓ(2k-ℓ-1)² ≤ 4(ℓ+1)*N
  -- Strategy: multiply goal by ℓ, prove ℓ²D² ≤ 4ℓ(ℓ+1)N via calc, cancel ℓ.
  suffices h_ℓ_mul : ℓ * (ℓ * (2 * k - ℓ - 1) * (2 * k - ℓ - 1))
                   ≤ ℓ * (4 * (ℓ + 1) * N) from
    Nat.le_of_mul_le_mul_left h_ℓ_mul hℓ_pos
  calc ℓ * (ℓ * (2 * k - ℓ - 1) * (2 * k - ℓ - 1))
      = ℓ * (2 * k - ℓ - 1) * (ℓ * (2 * k - ℓ - 1)) := by ring
    _ = m * 2 * (m * 2) := by rw [← hm_eq]
    _ ≤ 4 * m * (m + 1) := by nlinarith
    _ ≤ 4 * (ℓ * (ℓ + 1) * N) := by nlinarith [h_key]
    _ = ℓ * (4 * (ℓ + 1) * N) := by ring

/-- **Lindström Bound** (proven modulo 1 axiom):
    |A| ≤ √N + ⁴√N + 1 for any Sidon set A ⊆ [0,N]. -/
theorem lindstrom_bound (A : Finset ℕ) (N : ℕ)
    (hS : IsSidonSet A) (hA : A ⊆ Finset.range (N + 1)) (hN : 0 < N) :
    A.card ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 1 := by
  exact lindstrom_bound_extraction A N hS hA hN (fun ℓ hℓ_pos hℓ =>
    lindstrom_quadratic A N ℓ hS hA hℓ_pos hℓ)

end Erdos.Lindstrom
