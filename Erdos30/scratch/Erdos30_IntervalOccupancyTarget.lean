import Mathlib
import Erdos30_Sidon_Defs
import Erdos30_OrderedElements

open Finset Nat

namespace Erdos.Sidon

/-!
  Erdős #30 — Dense finite Sidon ordered-element interface
  =======================================================

  This file is a **research target**, not part of the imported core package.
  It records the current honest formal state of the dense-finite rigidity
  program after checking the literature scale.

  Local, axiom-free combinatorics previously lived here and has been lifted to
  `Erdos30_OrderedElements.lean`. What remains in this scratch file is exactly
  the content that depends on an external-theorem interface: the
  Balasubramanian–Dutta axiom and the theorems derived from it.

  Important:
  - The external ordered-element theorem is not proved locally.
  - The earlier sharp prefix-discrepancy axiom was removed because it was
    stronger than the published Balasubramanian–Dutta scale.
  - The file is kept in `scratch/` until a real local proof route is identified
    for the external interface.

  **Style note.** The axiom and its derived consequences are stated in ℝ-valued
  form (`Real.sqrt`, `Real.rpow`) to match the published paper. The local
  combinatorics imported from `Erdos30_OrderedElements.lean` are ℕ-valued; casts
  are performed at the point of use.
-/

/-!
### Perplexity Pre-Submission Gate — 2026-04-20

**Query:** Verify Balasubramanian–Dutta, *The m-th Element of a Sidon Set*
(`arXiv:2409.01986`, J. Number Theory Vol. 279, DOI
`10.1016/j.jnt.2025.07.007`) Theorem 3 statement, error exponents,
deficiency dependence, and hypotheses.

**Result: PASS — external theorem interface matches Theorem 3 up to O-constant
absorption.**

- Published form: `a_m = m · n^{1/2} + O(n^{7/8}) + O(L^{1/2} · n^{3/4})`,
  where `L = max(0, n^{1/2} - |A|)` and `m` ranges over `{1, ..., |A|}`.
- No macroscopic-index restriction; bound is uniform in `m`.
- Hypotheses: Sidon-ness and the deficiency relation.
- No known counterexamples, gaps, or sharper follow-ups found by the gate as of
  2026-04-20.

**Absorbed discrepancies:**
- Paper: `A ⊆ {1, ..., n}`; this file: `A ⊆ Finset.range (n+1)` = `{0, ..., n}`.
  Shift by `1` in the ambient range is absorbed into the O-constant.
- Paper: `L = max(0, √n - |A|)` with real square root; this file uses the
  integer-floor deficiency relation `A.card + L = Nat.sqrt n`.
  The floor difference is `< 1`, again absorbed into the O-constant.

**Future-proof note.** One proof route in the paper assumes `L ≤ n^(21/80)`,
but the main theorem statement does not. Any attempt to replace the external
interface with a local proof must handle both deficiency regimes.

**Lean-side pitfalls flagged by gate:**
- `n^(7/8)` and `n^(3/4)` are expressed here with `Real.rpow`.
- O-notation is non-constructive; this file packages the hidden constant as an
  existential `C`.
-/

/-- External theorem interface for Balasubramanian–Dutta,
*The m-th Element of a Sidon Set*.

Verified against the open arXiv version `arXiv:2409.01986` and the Journal of
Number Theory publication (Vol. 279, February 2026, DOI
`10.1016/j.jnt.2025.07.007`).

This is not a local proof in the present file. It records the literature-scale
ordered-element estimate that is currently justified for dense Sidon sets.
-/
axiom dense_sidon_ordered_element_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (_hDense : DenseSidonAtScale A n L)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) - (((i.1 + 1 : ℕ) : ℝ) * Real.sqrt (n : ℝ))|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)

/-- Literature consequence at the actual cutpoints `t = a_i`.

Combining the external Balasubramanian–Dutta theorem with the exact local
identity `|(A ∩ [0, a_i])| = i+1` gives a genuine prefix-count statement at the
ordered-element cutpoints. This is weaker than a uniform prefix discrepancy
theorem, but it is an honest consequence of the verified literature interface.
-/
theorem dense_sidon_prefix_cutpoint_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i : Fin A.card),
        |(orderedElement A i : ℝ) -
            ((intervalSlice A 0 (orderedElement A i)).card : ℝ) * Real.sqrt (n : ℝ)|
          ≤ C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i
  have hmain := hC hDense i
  have hcard :
      (((i.1 + 1 : ℕ) : ℝ)) =
        ((intervalSlice A 0 (orderedElement A i)).card : ℝ) := by
    norm_num [ordered_prefix_card_target A i]
  simpa [hcard, mul_comm, mul_left_comm, mul_assoc] using hmain

/-- Consecutive-gap corollary of the external ordered-element theorem.

Two applications of Balasubramanian–Dutta at indices `i` and `i+1`, stitched
by the triangle inequality, give that consecutive ordered elements differ by
`√n` up to twice the Balasubramanian–Dutta error. This is the stepping stone
to a nearby-prefix theorem for general `t`: once the gap is controlled, the
prefix count cannot jump more than `√n + error` across any single step, so
moving `t` off a cutpoint perturbs `|t − s(t)·√n|` by at most one gap.

No new axioms used; this is pure triangle inequality on
`dense_sidon_ordered_element_external`.
-/
theorem dense_sidon_consecutive_gap_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i : Fin A.card) (hi : i.1 + 1 < A.card),
        |((orderedElement A ⟨i.1 + 1, hi⟩ : ℝ) - (orderedElement A i : ℝ))
            - Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i hi
  set j : Fin A.card := ⟨i.1 + 1, hi⟩ with hj
  have hi_est := hC hDense i
  have hj_est := hC hDense j
  set a_i : ℝ := (orderedElement A i : ℝ)
  set a_j : ℝ := (orderedElement A j : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hi_bd : |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ E := hi_est
  have hj_bd : |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| ≤ E := hj_est
  have hjval : ((j.1 + 1 : ℕ) : ℝ) = ((i.1 + 1 : ℕ) : ℝ) + 1 := by
    push_cast [hj]
    ring
  have hj_bd' : |a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s| ≤ E := by
    simpa [hjval] using hj_bd
  have hkey :
      (a_j - a_i) - s =
        (a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s) -
          (a_i - ((i.1 + 1 : ℕ) : ℝ) * s) := by ring
  calc
    |(a_j - a_i) - s|
        = |(a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s) -
            (a_i - ((i.1 + 1 : ℕ) : ℝ) * s)| := by rw [hkey]
    _ ≤ |a_j - (((i.1 + 1 : ℕ) : ℝ) + 1) * s| +
          |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| := abs_sub _ _
    _ ≤ E + E := add_le_add hj_bd' hi_bd
    _ = 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
          2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
        show (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) +
              (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) = _
        ring

/-- Index-difference bound: multi-step generalization of the consecutive gap.

For any two ordered-element indices `i ≤ j`, the displacement `a_j − a_i` is
within `2·E` of `(j − i)·√n`, where `E` is the Balasubramanian–Dutta error.
Two applications of the axiom at `i` and `j` plus triangle inequality.

This subsumes `dense_sidon_consecutive_gap_external` (`j = i + 1` case) and
is the natural tool for bounding `t` that sits between two ordered elements:
once `t` is bracketed by `a_i ≤ t ≤ a_j`, interval control follows from this
lemma plus monotonicity of the prefix count.

The gap between two cutpoints `(j − i) · √n` grows linearly with index
separation, while the error stays at B–D scale — so widely separated cutpoints
give a sharper relative bound than consecutive ones.

No new axioms used.
-/
theorem dense_sidon_index_difference_external :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ {A : Finset ℕ} {n L : ℕ} (hDense : DenseSidonAtScale A n L)
        (i j : Fin A.card) (hij : i.1 ≤ j.1),
        |((orderedElement A j : ℝ) - (orderedElement A i : ℝ))
            - ((j.1 - i.1 : ℕ) : ℝ) * Real.sqrt (n : ℝ)|
          ≤ 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
  rcases dense_sidon_ordered_element_external with ⟨C, hCnonneg, hC⟩
  refine ⟨C, hCnonneg, ?_⟩
  intro A n L hDense i j hij
  have hi_bd := hC hDense i
  have hj_bd := hC hDense j
  set a_i : ℝ := (orderedElement A i : ℝ)
  set a_j : ℝ := (orderedElement A j : ℝ)
  set s : ℝ := Real.sqrt (n : ℝ)
  set E : ℝ := C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
      C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)
  have hi_abs : |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| ≤ E := hi_bd
  have hj_abs : |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| ≤ E := hj_bd
  have hdiff :
      ((j.1 - i.1 : ℕ) : ℝ) = ((j.1 + 1 : ℕ) : ℝ) - ((i.1 + 1 : ℕ) : ℝ) := by
    rw [Nat.cast_sub hij]
    push_cast
    ring
  have hkey :
      (a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s =
        (a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
          (a_i - ((i.1 + 1 : ℕ) : ℝ) * s) := by
    rw [hdiff]; ring
  calc
    |(a_j - a_i) - ((j.1 - i.1 : ℕ) : ℝ) * s|
        = |(a_j - ((j.1 + 1 : ℕ) : ℝ) * s) -
            (a_i - ((i.1 + 1 : ℕ) : ℝ) * s)| := by rw [hkey]
    _ ≤ |a_j - ((j.1 + 1 : ℕ) : ℝ) * s| +
          |a_i - ((i.1 + 1 : ℕ) : ℝ) * s| := abs_sub _ _
    _ ≤ E + E := add_le_add hj_abs hi_abs
    _ = 2 * C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
          2 * C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ) := by
        show (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) +
              (C * Real.rpow (n : ℝ) (7 / 8 : ℝ) +
            C * Real.sqrt (L : ℝ) * Real.rpow (n : ℝ) (3 / 4 : ℝ)) = _
        ring

end Erdos.Sidon
