import Erdos30_FaceField
import Erdos30_Sidon_Defs
import Erdos30_Singer57_Certificate

/-!
# Erdos #30 n = 57/58 Face/Field Branch Certificate

This file is a narrow finite certificate over packet-exported witnesses from
the 56--58 face handoff reports and the Singer57/N58 Lean certificate.

It does not prove Erdos #30, does not certify a full extremal enumeration, and
does not state an asymptotic theorem.  It only records exact exported-face facts:

* at `n = 57`, the mass and joint selected witnesses differ while sharing the
  same positive-difference skeleton;
* at `n = 58`, a second positive-difference skeleton appears through the prefix
  side, while the Pareto chain remains the translated `n = 57` endpoint chain.
-/

open Finset Nat

namespace Erdos30FaceField5758Certificate


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

/-! ## n = 57: field handoff inside one difference skeleton -/

def n57MassWitness : Finset Nat :=
  Erdos30Singer57Certificate.N57_3.toFinset

def n57JointWitness : Finset Nat :=
  Erdos30Singer57Certificate.N57_5.toFinset

def n57PrefixWinnerIndices : List Nat := [4, 5]
def n57MassWinnerIndices : List Nat := [3]
def n57JointWinnerIndices : List Nat := [5]
def n57ParetoWinnerIndices : List Nat := [3, 5]

theorem n57_mass_witness_card :
    n57MassWitness.card = 10 := by
  native_decide

theorem n57_joint_witness_card :
    n57JointWitness.card = 10 := by
  native_decide

theorem n57_mass_witness_in_range :
    n57MassWitness ⊆ Finset.range 58 := by
  native_decide

theorem n57_joint_witness_in_range :
    n57JointWitness ⊆ Finset.range 58 := by
  native_decide

theorem n57_mass_witness_sidon :
    Erdos.Sidon.IsSidonSet n57MassWitness := by
  native_decide

theorem n57_joint_witness_sidon :
    Erdos.Sidon.IsSidonSet n57JointWitness := by
  native_decide

theorem n57_mass_joint_distinct :
    n57MassWitness ≠ n57JointWitness := by
  native_decide

theorem n57_mass_joint_share_positive_difference_skeleton :
    Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N57_3) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N57_5) = true := by
  native_decide

theorem n57_joint_is_mass_plus_one :
    Erdos30Singer57Certificate.N57_5 = Erdos30Singer57Certificate.translate 1 Erdos30Singer57Certificate.N57_3 := by
  native_decide

theorem n57_field_winner_indices_match_handoff_packet :
    (n57PrefixWinnerIndices == [4, 5]) &&
    (n57MassWinnerIndices == [3]) &&
    (n57JointWinnerIndices == [5]) &&
    (n57ParetoWinnerIndices == [3, 5]) = true := by
  native_decide

theorem n57_mass_joint_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n57MassWitness n57JointWitness)
      (leftObservable n57MassWitness)
      (rightObservable n57JointWitness) :=
  pair_field_split n57_mass_joint_distinct

theorem n57_exported_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n57MassWitness n57JointWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n57_mass_joint_field_split

/-! ## n = 58: prefix-side branch, Pareto chain persistence -/

def n58PrefixBranchWitness : Finset Nat :=
  Erdos30Singer57Certificate.N58_2.toFinset

def n58MassParetoWitness : Finset Nat :=
  Erdos30Singer57Certificate.N58_7.toFinset

def n58TranslatedParetoWitness : Finset Nat :=
  Erdos30Singer57Certificate.N58_9.toFinset

def n58PrefixWinnerIndices : List Nat := [2, 8, 9]
def n58MassWinnerIndices : List Nat := [7]
def n58JointWinnerIndices : List Nat := [7]
def n58ParetoWinnerIndices : List Nat := [7, 9]

theorem n58_prefix_branch_card :
    n58PrefixBranchWitness.card = 10 := by
  native_decide

theorem n58_mass_pareto_card :
    n58MassParetoWitness.card = 10 := by
  native_decide

theorem n58_translated_pareto_card :
    n58TranslatedParetoWitness.card = 10 := by
  native_decide

theorem n58_prefix_branch_in_range :
    n58PrefixBranchWitness ⊆ Finset.range 59 := by
  native_decide

theorem n58_mass_pareto_in_range :
    n58MassParetoWitness ⊆ Finset.range 59 := by
  native_decide

theorem n58_translated_pareto_in_range :
    n58TranslatedParetoWitness ⊆ Finset.range 59 := by
  native_decide

theorem n58_prefix_branch_sidon :
    Erdos.Sidon.IsSidonSet n58PrefixBranchWitness := by
  native_decide

theorem n58_mass_pareto_sidon :
    Erdos.Sidon.IsSidonSet n58MassParetoWitness := by
  native_decide

theorem n58_translated_pareto_sidon :
    Erdos.Sidon.IsSidonSet n58TranslatedParetoWitness := by
  native_decide

theorem n58_prefix_branch_differs_from_mass_pareto :
    n58PrefixBranchWitness ≠ n58MassParetoWitness := by
  native_decide

theorem n58_prefix_branch_has_new_positive_difference_skeleton :
    Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_2) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_7) = false := by
  native_decide

theorem n58_prefix_branch_matches_companion_skeleton :
    Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_2) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_3) = true := by
  native_decide

theorem n58_pareto_candidates_share_translated_skeleton :
    Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_7) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_9) = true := by
  native_decide

theorem n58_pareto_candidate_is_translation :
    Erdos30Singer57Certificate.N58_9 = Erdos30Singer57Certificate.translate 1 Erdos30Singer57Certificate.N58_7 := by
  native_decide

theorem n58_pareto_chain_persists_from_n57 :
    (Erdos30Singer57Certificate.N58_7 == Erdos30Singer57Certificate.N57_5) &&
    (Erdos30Singer57Certificate.N58_9 == Erdos30Singer57Certificate.translate 1 Erdos30Singer57Certificate.N57_5) &&
    (Erdos30Singer57Certificate.N58_9 == Erdos30Singer57Certificate.translate 2 Erdos30Singer57Certificate.N57_3) = true := by
  native_decide

theorem n58_field_winner_indices_match_branch_packet :
    (n58PrefixWinnerIndices == [2, 8, 9]) &&
    (n58MassWinnerIndices == [7]) &&
    (n58JointWinnerIndices == [7]) &&
    (n58ParetoWinnerIndices == [7, 9]) &&
    Erdos30Singer57Certificate.containsNat 2 n58PrefixWinnerIndices &&
    (Erdos30Singer57Certificate.containsNat 2 n58ParetoWinnerIndices == false) = true := by
  native_decide

theorem n58_prefix_mass_field_split :
    Erdos.Collider.FieldSplit
      (pairFamily n58PrefixBranchWitness n58MassParetoWitness)
      (leftObservable n58PrefixBranchWitness)
      (rightObservable n58MassParetoWitness) :=
  pair_field_split n58_prefix_branch_differs_from_mass_pareto

theorem n58_exported_pair_has_at_least_two_witnesses :
    2 ≤ (pairFamily n58PrefixBranchWitness n58MassParetoWitness).card :=
  Erdos.Collider.fieldSplit_card_two_le n58_prefix_mass_field_split

/-! ## Bounded combined finite certificate -/

theorem finite_57_58_branch_field_response_certificate :
    (Erdos30Singer57Certificate.allListNat Erdos30Singer57Certificate.intervalSidon10 Erdos30Singer57Certificate.n57Witnesses) &&
    (Erdos30Singer57Certificate.allListNat Erdos30Singer57Certificate.intervalSidon10 Erdos30Singer57Certificate.n58Witnesses) &&
    (Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N57_3) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N57_5)) &&
    (Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_2) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_7) == false) &&
    (Erdos30Singer57Certificate.sameNatSet (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_7) (Erdos30Singer57Certificate.positiveDiffs Erdos30Singer57Certificate.N58_9)) &&
    (Erdos30Singer57Certificate.N57_5 == Erdos30Singer57Certificate.translate 1 Erdos30Singer57Certificate.N57_3) &&
    (Erdos30Singer57Certificate.N58_9 == Erdos30Singer57Certificate.translate 1 Erdos30Singer57Certificate.N58_7) &&
    (Erdos30Singer57Certificate.N58_7 == Erdos30Singer57Certificate.N57_5) &&
    (Erdos30Singer57Certificate.containsNat 2 n58PrefixWinnerIndices) &&
    (Erdos30Singer57Certificate.containsNat 2 n58ParetoWinnerIndices == false) = true := by
  native_decide

end Erdos30FaceField5758Certificate
