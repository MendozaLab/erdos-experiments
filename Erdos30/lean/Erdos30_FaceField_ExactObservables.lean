import Mathlib

/-!
# Exact Observable Helpers for Erdos #30 Face/Field Certificates

This file defines the arithmetic observables used by the finite face/field
certificates without importing floating-point packet scores.

The density-adjusted mass observable is fully computable as an integer:

`|2 * sum(A) - n * (|A| + 1)|`.

This is twice the packet's density-adjusted mass deviation, so it gives the
same ordering while avoiding rationals and floats. The prefix definitions here
are an exact scaffold for the next climb; minimizer theorems over prefix
residuals require algebraic comparison over square roots and are intentionally
not asserted here.
-/

namespace Erdos30FaceFieldExactObservables

def witnessSum (A : Finset Nat) : Nat :=
  A.sum (fun a => a)

def prefixCount (A : Finset Nat) (t : Nat) : Nat :=
  (A.filter (fun a => a ≤ t)).card

def densityAdjustedMassCenterTwice (n : Nat) (A : Finset Nat) : Nat :=
  n * (A.card + 1)

def densityAdjustedMassTwice (n : Nat) (A : Finset Nat) : Nat :=
  let lhs := 2 * witnessSum A
  let rhs := densityAdjustedMassCenterTwice n A
  if lhs ≤ rhs then rhs - lhs else lhs - rhs

noncomputable def prefixDeviationAt (n : Nat) (A : Finset Nat) (t : Nat) : ℝ :=
  |(t : ℝ) - (prefixCount A t : ℝ) * Real.sqrt (n : ℝ)|

noncomputable def prefixDrift (n : Nat) (A : Finset Nat) : ℝ :=
  max |(A.card : ℝ) - Real.sqrt (n : ℝ)| 1 * Real.sqrt (n : ℝ)

noncomputable def prefixResidualAt (n : Nat) (A : Finset Nat) (t : Nat) : ℝ :=
  max 0 (prefixDeviationAt n A t - prefixDrift n A)

theorem prefixResidualAt_nonneg (n : Nat) (A : Finset Nat) (t : Nat) :
    0 ≤ prefixResidualAt n A t := by
  unfold prefixResidualAt
  exact le_max_left 0 _

def fullPrefixResidualZero (n : Nat) (A : Finset Nat) : Prop :=
  ∀ t : Nat, t ≤ n -> prefixResidualAt n A t = 0

noncomputable def fullPrefixResidualMax (n : Nat) (A : Finset Nat) : ℝ :=
  (Finset.range (n + 1)).sup' (by exact ⟨0, by simp⟩)
    (fun t => prefixResidualAt n A t)

theorem fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    {n : Nat} {A : Finset Nat}
    (h : fullPrefixResidualZero n A) :
    fullPrefixResidualMax n A = 0 := by
  unfold fullPrefixResidualMax
  apply le_antisymm
  · apply Finset.sup'_le
    intro t ht
    have htle : t ≤ n := by
      simp at ht
      omega
    rw [h t htle]
  · have hmem : 0 ∈ Finset.range (n + 1) := by
      simp
    have hle := Finset.le_sup' (fun t => prefixResidualAt n A t) hmem
    rw [h 0 (Nat.zero_le n)] at hle
    exact hle

theorem fullPrefixResidualMax_pos_of_pos_at
    {n : Nat} {A : Finset Nat} {t : Nat}
    (ht : t ≤ n)
    (hpos : 0 < prefixResidualAt n A t) :
    0 < fullPrefixResidualMax n A := by
  unfold fullPrefixResidualMax
  have hmem : t ∈ Finset.range (n + 1) := by
    simp
    omega
  have hle := Finset.le_sup' (fun t => prefixResidualAt n A t) hmem
  exact lt_of_lt_of_le hpos hle

theorem fullPrefixResidualMax_nonneg (n : Nat) (A : Finset Nat) :
    0 ≤ fullPrefixResidualMax n A := by
  unfold fullPrefixResidualMax
  have hmem : 0 ∈ Finset.range (n + 1) := by
    simp
  have hle := Finset.le_sup' (fun t => prefixResidualAt n A t) hmem
  exact le_trans (prefixResidualAt_nonneg n A 0) hle

theorem fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    {n : Nat} {A : Finset Nat} {t : Nat} {v : ℝ}
    (ht : t ≤ n)
    (hupper : ∀ u : Nat, u ≤ n -> prefixResidualAt n A u ≤ v)
    (hat : prefixResidualAt n A t = v) :
    fullPrefixResidualMax n A = v := by
  unfold fullPrefixResidualMax
  apply le_antisymm
  · apply Finset.sup'_le
    intro u hu
    have hule : u ≤ n := by
      simp at hu
      omega
    exact hupper u hule
  · have hmem : t ∈ Finset.range (n + 1) := by
      simp
      omega
    have hle := Finset.le_sup' (fun u => prefixResidualAt n A u) hmem
    rwa [hat] at hle

noncomputable def prefixResidualProbeCard10 (n t prefix_count : Nat) : ℝ :=
  max 0 (|(t : ℝ) - (prefix_count : ℝ) * Real.sqrt (n : ℝ)| -
    (10 * Real.sqrt (n : ℝ) - n))

noncomputable def prefixProbeMassJointKeyCard10
    (n t prefix_count mass_twice : Nat) : ℝ :=
  prefixResidualProbeCard10 n t prefix_count * Real.sqrt (n : ℝ) +
    (mass_twice : ℝ) / 2

noncomputable def densityAdjustedMassDev (n : Nat) (A : Finset Nat) : ℝ :=
  (densityAdjustedMassTwice n A : ℝ) / 2

noncomputable def jointKeyAt (n : Nat) (A : Finset Nat) (t : Nat) : ℝ :=
  prefixResidualAt n A t * Real.sqrt (n : ℝ) + densityAdjustedMassDev n A

end Erdos30FaceFieldExactObservables
