import Erdos30_FaceField
import Erdos30_Sidon_Defs

/-!
# Erdos #30 Sidon Face/Field Window Certificate

This file instantiates the abstract face/field lemma over the Atheneum Gauntlet
Sidon window `n = 20..30`.

Source packet:
`ATH-GAUNTLET-SIDON-20-30-2026-04-30_RESULTS.json`.

For each row, the packet-exported prefix-field and mass-field witnesses are
recorded as literal finite sets. Lean checks that each witness is a Sidon set,
has the packet cardinality, lies in `[0,n]`, differs from the other witness, and
therefore gives a two-point exposed witness family under the abstract
face/field lemma.

This is not a theorem about the full extremal face and not a proof of Erdos #30.
It is an exact finite certificate over the packet-exported exposed witnesses.
-/

open Finset Nat

namespace Erdos30FaceFieldWindowCertificate

def pairFamily (x y : Finset Nat) : Finset (Finset Nat) :=
  ({x, y} : Finset (Finset Nat))

def leftObservable (x A : Finset Nat) : Nat :=
  if A = x then 0 else 1

def rightObservable (y A : Finset Nat) : Nat :=
  if A = y then 0 else 1

theorem pair_field_split {x y : Finset Nat} (hxy : x ≠ y) :
    Erdos.Collider.FieldSplit (pairFamily x y)
      (leftObservable x) (rightObservable y) := by
  refine ⟨x, y, ?_, ?_, hxy⟩
  · constructor
    · simp [pairFamily]
    · intro z hz
      simp [leftObservable]
  · constructor
    · simp [pairFamily]
    · intro z hz
      simp [rightObservable]

/-! ## n = 20 -/

def n20PrefixWitness : Finset Nat :=
  ([0, 3, 7, 12, 18, 20] : List Nat).toFinset

def n20MassWitness : Finset Nat :=
  ([0, 3, 13, 15, 19, 20] : List Nat).toFinset

theorem n20_prefix_card :
    n20PrefixWitness.card = 6 := by
  native_decide

theorem n20_mass_card :
    n20MassWitness.card = 6 := by
  native_decide

theorem n20_prefix_in_range :
    n20PrefixWitness ⊆ Finset.range 21 := by
  native_decide

theorem n20_mass_in_range :
    n20MassWitness ⊆ Finset.range 21 := by
  native_decide

theorem n20_prefix_sidon :
    Erdos.Sidon.IsSidonSet n20PrefixWitness := by
  native_decide

theorem n20_mass_sidon :
    Erdos.Sidon.IsSidonSet n20MassWitness := by
  native_decide

theorem n20_prefix_mass_distinct :
    n20PrefixWitness ≠ n20MassWitness := by
  native_decide

theorem n20_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n20PrefixWitness n20MassWitness)
      (leftObservable n20PrefixWitness)
      (rightObservable n20MassWitness) :=
  pair_field_split n20_prefix_mass_distinct

theorem n20_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n20PrefixWitness n20MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n20_prefix_mass_field_split

/-! ## n = 21 -/

def n21PrefixWitness : Finset Nat :=
  ([0, 3, 8, 15, 17, 21] : List Nat).toFinset

def n21MassWitness : Finset Nat :=
  ([0, 3, 14, 16, 20, 21] : List Nat).toFinset

theorem n21_prefix_card :
    n21PrefixWitness.card = 6 := by
  native_decide

theorem n21_mass_card :
    n21MassWitness.card = 6 := by
  native_decide

theorem n21_prefix_in_range :
    n21PrefixWitness ⊆ Finset.range 22 := by
  native_decide

theorem n21_mass_in_range :
    n21MassWitness ⊆ Finset.range 22 := by
  native_decide

theorem n21_prefix_sidon :
    Erdos.Sidon.IsSidonSet n21PrefixWitness := by
  native_decide

theorem n21_mass_sidon :
    Erdos.Sidon.IsSidonSet n21MassWitness := by
  native_decide

theorem n21_prefix_mass_distinct :
    n21PrefixWitness ≠ n21MassWitness := by
  native_decide

theorem n21_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n21PrefixWitness n21MassWitness)
      (leftObservable n21PrefixWitness)
      (rightObservable n21MassWitness) :=
  pair_field_split n21_prefix_mass_distinct

theorem n21_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n21PrefixWitness n21MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n21_prefix_mass_field_split

/-! ## n = 22 -/

def n22PrefixWitness : Finset Nat :=
  ([0, 4, 9, 19, 20, 22] : List Nat).toFinset

def n22MassWitness : Finset Nat :=
  ([0, 4, 14, 16, 21, 22] : List Nat).toFinset

theorem n22_prefix_card :
    n22PrefixWitness.card = 6 := by
  native_decide

theorem n22_mass_card :
    n22MassWitness.card = 6 := by
  native_decide

theorem n22_prefix_in_range :
    n22PrefixWitness ⊆ Finset.range 23 := by
  native_decide

theorem n22_mass_in_range :
    n22MassWitness ⊆ Finset.range 23 := by
  native_decide

theorem n22_prefix_sidon :
    Erdos.Sidon.IsSidonSet n22PrefixWitness := by
  native_decide

theorem n22_mass_sidon :
    Erdos.Sidon.IsSidonSet n22MassWitness := by
  native_decide

theorem n22_prefix_mass_distinct :
    n22PrefixWitness ≠ n22MassWitness := by
  native_decide

theorem n22_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n22PrefixWitness n22MassWitness)
      (leftObservable n22PrefixWitness)
      (rightObservable n22MassWitness) :=
  pair_field_split n22_prefix_mass_distinct

theorem n22_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n22PrefixWitness n22MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n22_prefix_mass_field_split

/-! ## n = 23 -/

def n23PrefixWitness : Finset Nat :=
  ([0, 4, 9, 15, 22, 23] : List Nat).toFinset

def n23MassWitness : Finset Nat :=
  ([0, 2, 14, 19, 22, 23] : List Nat).toFinset

theorem n23_prefix_card :
    n23PrefixWitness.card = 6 := by
  native_decide

theorem n23_mass_card :
    n23MassWitness.card = 6 := by
  native_decide

theorem n23_prefix_in_range :
    n23PrefixWitness ⊆ Finset.range 24 := by
  native_decide

theorem n23_mass_in_range :
    n23MassWitness ⊆ Finset.range 24 := by
  native_decide

theorem n23_prefix_sidon :
    Erdos.Sidon.IsSidonSet n23PrefixWitness := by
  native_decide

theorem n23_mass_sidon :
    Erdos.Sidon.IsSidonSet n23MassWitness := by
  native_decide

theorem n23_prefix_mass_distinct :
    n23PrefixWitness ≠ n23MassWitness := by
  native_decide

theorem n23_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n23PrefixWitness n23MassWitness)
      (leftObservable n23PrefixWitness)
      (rightObservable n23MassWitness) :=
  pair_field_split n23_prefix_mass_distinct

theorem n23_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n23PrefixWitness n23MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n23_prefix_mass_field_split

/-! ## n = 24 -/

def n24PrefixWitness : Finset Nat :=
  ([0, 5, 11, 15, 23, 24] : List Nat).toFinset

def n24MassWitness : Finset Nat :=
  ([0, 2, 15, 20, 23, 24] : List Nat).toFinset

theorem n24_prefix_card :
    n24PrefixWitness.card = 6 := by
  native_decide

theorem n24_mass_card :
    n24MassWitness.card = 6 := by
  native_decide

theorem n24_prefix_in_range :
    n24PrefixWitness ⊆ Finset.range 25 := by
  native_decide

theorem n24_mass_in_range :
    n24MassWitness ⊆ Finset.range 25 := by
  native_decide

theorem n24_prefix_sidon :
    Erdos.Sidon.IsSidonSet n24PrefixWitness := by
  native_decide

theorem n24_mass_sidon :
    Erdos.Sidon.IsSidonSet n24MassWitness := by
  native_decide

theorem n24_prefix_mass_distinct :
    n24PrefixWitness ≠ n24MassWitness := by
  native_decide

theorem n24_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n24PrefixWitness n24MassWitness)
      (leftObservable n24PrefixWitness)
      (rightObservable n24MassWitness) :=
  pair_field_split n24_prefix_mass_distinct

theorem n24_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n24PrefixWitness n24MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n24_prefix_mass_field_split

/-! ## n = 25 -/

def n25PrefixWitness : Finset Nat :=
  ([0, 1, 7, 11, 20, 23, 25] : List Nat).toFinset

def n25MassWitness : Finset Nat :=
  ([0, 4, 9, 15, 22, 23, 25] : List Nat).toFinset

theorem n25_prefix_card :
    n25PrefixWitness.card = 7 := by
  native_decide

theorem n25_mass_card :
    n25MassWitness.card = 7 := by
  native_decide

theorem n25_prefix_in_range :
    n25PrefixWitness ⊆ Finset.range 26 := by
  native_decide

theorem n25_mass_in_range :
    n25MassWitness ⊆ Finset.range 26 := by
  native_decide

theorem n25_prefix_sidon :
    Erdos.Sidon.IsSidonSet n25PrefixWitness := by
  native_decide

theorem n25_mass_sidon :
    Erdos.Sidon.IsSidonSet n25MassWitness := by
  native_decide

theorem n25_prefix_mass_distinct :
    n25PrefixWitness ≠ n25MassWitness := by
  native_decide

theorem n25_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n25PrefixWitness n25MassWitness)
      (leftObservable n25PrefixWitness)
      (rightObservable n25MassWitness) :=
  pair_field_split n25_prefix_mass_distinct

theorem n25_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n25PrefixWitness n25MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n25_prefix_mass_field_split

/-! ## n = 26 -/

def n26PrefixWitness : Finset Nat :=
  ([0, 1, 6, 14, 17, 24, 26] : List Nat).toFinset

def n26MassWitness : Finset Nat :=
  ([0, 2, 12, 18, 21, 25, 26] : List Nat).toFinset

theorem n26_prefix_card :
    n26PrefixWitness.card = 7 := by
  native_decide

theorem n26_mass_card :
    n26MassWitness.card = 7 := by
  native_decide

theorem n26_prefix_in_range :
    n26PrefixWitness ⊆ Finset.range 27 := by
  native_decide

theorem n26_mass_in_range :
    n26MassWitness ⊆ Finset.range 27 := by
  native_decide

theorem n26_prefix_sidon :
    Erdos.Sidon.IsSidonSet n26PrefixWitness := by
  native_decide

theorem n26_mass_sidon :
    Erdos.Sidon.IsSidonSet n26MassWitness := by
  native_decide

theorem n26_prefix_mass_distinct :
    n26PrefixWitness ≠ n26MassWitness := by
  native_decide

theorem n26_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n26PrefixWitness n26MassWitness)
      (leftObservable n26PrefixWitness)
      (rightObservable n26MassWitness) :=
  pair_field_split n26_prefix_mass_distinct

theorem n26_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n26PrefixWitness n26MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n26_prefix_mass_field_split

/-! ## n = 27 -/

def n27PrefixWitness : Finset Nat :=
  ([0, 2, 9, 15, 23, 26, 27] : List Nat).toFinset

def n27MassWitness : Finset Nat :=
  ([0, 2, 12, 18, 23, 26, 27] : List Nat).toFinset

theorem n27_prefix_card :
    n27PrefixWitness.card = 7 := by
  native_decide

theorem n27_mass_card :
    n27MassWitness.card = 7 := by
  native_decide

theorem n27_prefix_in_range :
    n27PrefixWitness ⊆ Finset.range 28 := by
  native_decide

theorem n27_mass_in_range :
    n27MassWitness ⊆ Finset.range 28 := by
  native_decide

theorem n27_prefix_sidon :
    Erdos.Sidon.IsSidonSet n27PrefixWitness := by
  native_decide

theorem n27_mass_sidon :
    Erdos.Sidon.IsSidonSet n27MassWitness := by
  native_decide

theorem n27_prefix_mass_distinct :
    n27PrefixWitness ≠ n27MassWitness := by
  native_decide

theorem n27_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n27PrefixWitness n27MassWitness)
      (leftObservable n27PrefixWitness)
      (rightObservable n27MassWitness) :=
  pair_field_split n27_prefix_mass_distinct

theorem n27_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n27PrefixWitness n27MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n27_prefix_mass_field_split

/-! ## n = 28 -/

def n28PrefixWitness : Finset Nat :=
  ([0, 2, 7, 16, 24, 27, 28] : List Nat).toFinset

def n28MassWitness : Finset Nat :=
  ([0, 4, 15, 18, 20, 27, 28] : List Nat).toFinset

theorem n28_prefix_card :
    n28PrefixWitness.card = 7 := by
  native_decide

theorem n28_mass_card :
    n28MassWitness.card = 7 := by
  native_decide

theorem n28_prefix_in_range :
    n28PrefixWitness ⊆ Finset.range 29 := by
  native_decide

theorem n28_mass_in_range :
    n28MassWitness ⊆ Finset.range 29 := by
  native_decide

theorem n28_prefix_sidon :
    Erdos.Sidon.IsSidonSet n28PrefixWitness := by
  native_decide

theorem n28_mass_sidon :
    Erdos.Sidon.IsSidonSet n28MassWitness := by
  native_decide

theorem n28_prefix_mass_distinct :
    n28PrefixWitness ≠ n28MassWitness := by
  native_decide

theorem n28_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n28PrefixWitness n28MassWitness)
      (leftObservable n28PrefixWitness)
      (rightObservable n28MassWitness) :=
  pair_field_split n28_prefix_mass_distinct

theorem n28_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n28PrefixWitness n28MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n28_prefix_mass_field_split

/-! ## n = 29 -/

def n29PrefixWitness : Finset Nat :=
  ([0, 3, 8, 15, 19, 28, 29] : List Nat).toFinset

def n29MassWitness : Finset Nat :=
  ([0, 3, 12, 22, 23, 27, 29] : List Nat).toFinset

theorem n29_prefix_card :
    n29PrefixWitness.card = 7 := by
  native_decide

theorem n29_mass_card :
    n29MassWitness.card = 7 := by
  native_decide

theorem n29_prefix_in_range :
    n29PrefixWitness ⊆ Finset.range 30 := by
  native_decide

theorem n29_mass_in_range :
    n29MassWitness ⊆ Finset.range 30 := by
  native_decide

theorem n29_prefix_sidon :
    Erdos.Sidon.IsSidonSet n29PrefixWitness := by
  native_decide

theorem n29_mass_sidon :
    Erdos.Sidon.IsSidonSet n29MassWitness := by
  native_decide

theorem n29_prefix_mass_distinct :
    n29PrefixWitness ≠ n29MassWitness := by
  native_decide

theorem n29_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n29PrefixWitness n29MassWitness)
      (leftObservable n29PrefixWitness)
      (rightObservable n29MassWitness) :=
  pair_field_split n29_prefix_mass_distinct

theorem n29_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n29PrefixWitness n29MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n29_prefix_mass_field_split

/-! ## n = 30 -/

def n30PrefixWitness : Finset Nat :=
  ([0, 3, 9, 16, 20, 28, 30] : List Nat).toFinset

def n30MassWitness : Finset Nat :=
  ([0, 1, 16, 21, 24, 28, 30] : List Nat).toFinset

theorem n30_prefix_card :
    n30PrefixWitness.card = 7 := by
  native_decide

theorem n30_mass_card :
    n30MassWitness.card = 7 := by
  native_decide

theorem n30_prefix_in_range :
    n30PrefixWitness ⊆ Finset.range 31 := by
  native_decide

theorem n30_mass_in_range :
    n30MassWitness ⊆ Finset.range 31 := by
  native_decide

theorem n30_prefix_sidon :
    Erdos.Sidon.IsSidonSet n30PrefixWitness := by
  native_decide

theorem n30_mass_sidon :
    Erdos.Sidon.IsSidonSet n30MassWitness := by
  native_decide

theorem n30_prefix_mass_distinct :
    n30PrefixWitness ≠ n30MassWitness := by
  native_decide

theorem n30_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n30PrefixWitness n30MassWitness)
      (leftObservable n30PrefixWitness)
      (rightObservable n30MassWitness) :=
  pair_field_split n30_prefix_mass_distinct

theorem n30_exposed_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n30PrefixWitness n30MassWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n30_prefix_mass_field_split

theorem sidon_window_20_30_all_prefix_mass_split :
    Erdos.Collider.FieldSplit
      (pairFamily n20PrefixWitness n20MassWitness)
      (leftObservable n20PrefixWitness)
      (rightObservable n20MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n21PrefixWitness n21MassWitness)
      (leftObservable n21PrefixWitness)
      (rightObservable n21MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n22PrefixWitness n22MassWitness)
      (leftObservable n22PrefixWitness)
      (rightObservable n22MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n23PrefixWitness n23MassWitness)
      (leftObservable n23PrefixWitness)
      (rightObservable n23MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n24PrefixWitness n24MassWitness)
      (leftObservable n24PrefixWitness)
      (rightObservable n24MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n25PrefixWitness n25MassWitness)
      (leftObservable n25PrefixWitness)
      (rightObservable n25MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n26PrefixWitness n26MassWitness)
      (leftObservable n26PrefixWitness)
      (rightObservable n26MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n27PrefixWitness n27MassWitness)
      (leftObservable n27PrefixWitness)
      (rightObservable n27MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n28PrefixWitness n28MassWitness)
      (leftObservable n28PrefixWitness)
      (rightObservable n28MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n29PrefixWitness n29MassWitness)
      (leftObservable n29PrefixWitness)
      (rightObservable n29MassWitness) ∧
    Erdos.Collider.FieldSplit
      (pairFamily n30PrefixWitness n30MassWitness)
      (leftObservable n30PrefixWitness)
      (rightObservable n30MassWitness) := by
  exact ⟨n20_prefix_mass_field_split, ⟨n21_prefix_mass_field_split, ⟨n22_prefix_mass_field_split, ⟨n23_prefix_mass_field_split, ⟨n24_prefix_mass_field_split, ⟨n25_prefix_mass_field_split, ⟨n26_prefix_mass_field_split, ⟨n27_prefix_mass_field_split, ⟨n28_prefix_mass_field_split, ⟨n29_prefix_mass_field_split, n30_prefix_mass_field_split⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem sidon_window_20_30_all_exposed_pairs_have_two_witnesses :
    2 ≤ (pairFamily n20PrefixWitness n20MassWitness).card ∧
    2 ≤ (pairFamily n21PrefixWitness n21MassWitness).card ∧
    2 ≤ (pairFamily n22PrefixWitness n22MassWitness).card ∧
    2 ≤ (pairFamily n23PrefixWitness n23MassWitness).card ∧
    2 ≤ (pairFamily n24PrefixWitness n24MassWitness).card ∧
    2 ≤ (pairFamily n25PrefixWitness n25MassWitness).card ∧
    2 ≤ (pairFamily n26PrefixWitness n26MassWitness).card ∧
    2 ≤ (pairFamily n27PrefixWitness n27MassWitness).card ∧
    2 ≤ (pairFamily n28PrefixWitness n28MassWitness).card ∧
    2 ≤ (pairFamily n29PrefixWitness n29MassWitness).card ∧
    2 ≤ (pairFamily n30PrefixWitness n30MassWitness).card := by
  exact ⟨n20_exposed_pair_has_at_least_two_witnesses, ⟨n21_exposed_pair_has_at_least_two_witnesses, ⟨n22_exposed_pair_has_at_least_two_witnesses, ⟨n23_exposed_pair_has_at_least_two_witnesses, ⟨n24_exposed_pair_has_at_least_two_witnesses, ⟨n25_exposed_pair_has_at_least_two_witnesses, ⟨n26_exposed_pair_has_at_least_two_witnesses, ⟨n27_exposed_pair_has_at_least_two_witnesses, ⟨n28_exposed_pair_has_at_least_two_witnesses, ⟨n29_exposed_pair_has_at_least_two_witnesses, n30_exposed_pair_has_at_least_two_witnesses⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

end Erdos30FaceFieldWindowCertificate
