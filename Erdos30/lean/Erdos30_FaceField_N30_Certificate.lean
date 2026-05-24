import Erdos30_FaceField
import Erdos30_Sidon_Defs

/-!
# Erdos #30 n=30 Face/Field Certificate

This file instantiates the abstract face/field lemma on one exact finite row
from the Atheneum Gauntlet Sidon packet:

`ATH-GAUNTLET-SIDON-20-30-2026-04-30_RESULTS.json`, row `n = 30`.

It is deliberately a finite certificate, not a theorem about the Sidon problem.
The full exact face at n=30 has 1618 ground states in the packet.  This file
only records three packet-exported field-exposed witnesses and proves that the
abstract exposed-face lemma applies to this certified witness family.
-/

open Finset Nat

namespace Erdos30FaceFieldN30Certificate

def prefixWitness : Finset Nat :=
  ([0, 3, 9, 16, 20, 28, 30] : List Nat).toFinset

def massWitness : Finset Nat :=
  ([0, 1, 16, 21, 24, 28, 30] : List Nat).toFinset

def jointWitness : Finset Nat :=
  ([0, 3, 17, 19, 25, 26, 30] : List Nat).toFinset

def exposedWitnessFamily : Finset (Finset Nat) :=
  ({prefixWitness, massWitness, jointWitness} : Finset (Finset Nat))

def n30PrefixObservable (A : Finset Nat) : Nat :=
  if A = prefixWitness ∨ A = jointWitness then 0 else 1

def n30MassObservable (A : Finset Nat) : Nat :=
  if A = massWitness ∨ A = jointWitness then 0 else 1

theorem prefixWitness_card :
    prefixWitness.card = 7 := by
  native_decide

theorem massWitness_card :
    massWitness.card = 7 := by
  native_decide

theorem jointWitness_card :
    jointWitness.card = 7 := by
  native_decide

theorem prefixWitness_in_range31 :
    prefixWitness ⊆ Finset.range 31 := by
  native_decide

theorem massWitness_in_range31 :
    massWitness ⊆ Finset.range 31 := by
  native_decide

theorem jointWitness_in_range31 :
    jointWitness ⊆ Finset.range 31 := by
  native_decide

theorem prefixWitness_sidon :
    Erdos.Sidon.IsSidonSet prefixWitness := by
  native_decide

theorem massWitness_sidon :
    Erdos.Sidon.IsSidonSet massWitness := by
  native_decide

theorem jointWitness_sidon :
    Erdos.Sidon.IsSidonSet jointWitness := by
  native_decide

theorem n30_prefix_mass_distinct :
    prefixWitness ≠ massWitness := by
  native_decide

theorem n30_prefix_joint_distinct :
    prefixWitness ≠ jointWitness := by
  native_decide

theorem n30_mass_joint_distinct :
    massWitness ≠ jointWitness := by
  native_decide

theorem n30_prefix_mass_field_split :
    Erdos.Collider.FieldSplit exposedWitnessFamily
      n30PrefixObservable n30MassObservable := by
  refine ⟨prefixWitness, massWitness, ?_, ?_, n30_prefix_mass_distinct⟩
  · constructor
    · simp [exposedWitnessFamily]
    · intro y hy
      simp only [exposedWitnessFamily, Finset.mem_insert, Finset.mem_singleton] at hy
      rcases hy with rfl | rfl | rfl
      · simp [n30PrefixObservable]
      · simp [n30PrefixObservable]
      · simp [n30PrefixObservable]
  · constructor
    · simp [exposedWitnessFamily]
    · intro y hy
      simp only [exposedWitnessFamily, Finset.mem_insert, Finset.mem_singleton] at hy
      rcases hy with rfl | rfl | rfl
      · simp [n30MassObservable, n30_prefix_mass_distinct]
      · simp [n30MassObservable]
      · simp [n30MassObservable]

theorem n30_exposed_family_has_at_least_two_witnesses :
    2 ≤ exposedWitnessFamily.card :=
  Erdos.Collider.fieldSplit_card_two_le n30_prefix_mass_field_split

theorem n30_exposed_family_not_singleton :
    ¬ exposedWitnessFamily.card ≤ 1 :=
  Erdos.Collider.fieldSplit_not_card_le_one n30_prefix_mass_field_split

end Erdos30FaceFieldN30Certificate
