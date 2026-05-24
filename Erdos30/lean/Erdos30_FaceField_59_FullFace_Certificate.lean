import Erdos30_FaceField
import Erdos30_FaceField_ExactObservables
import Erdos30_Sidon_Defs
import Mathlib

/-!
# Erdos #30 n = 59 Scalar Full-Prefix Joint Microcertificate

This file instantiates the reusable scalar-joint minimizer transfer over the
complete exported `n = 59` ground face from:

`EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01_RESULTS.json`.

It is a finite exported-face certificate only. It is not a proof of Erdos #30,
not a proof that the export was exhaustive beyond the source packet's exact
count claim, and not an asymptotic theorem.
-/

open Finset Nat

namespace Erdos30FaceField59FullFaceCertificate

/-! ## n = 59: complete exported face -/

def n59W0 : Finset Nat :=
  ([0, 1, 6, 10, 23, 26, 34, 41, 53, 55] : List Nat).toFinset

def n59W1 : Finset Nat :=
  ([0, 2, 8, 13, 22, 23, 40, 52, 56, 59] : List Nat).toFinset

def n59W2 : Finset Nat :=
  ([0, 2, 14, 21, 29, 32, 45, 49, 54, 55] : List Nat).toFinset

def n59W3 : Finset Nat :=
  ([0, 2, 15, 21, 22, 32, 46, 50, 55, 58] : List Nat).toFinset

def n59W4 : Finset Nat :=
  ([0, 3, 4, 10, 22, 30, 43, 45, 54, 59] : List Nat).toFinset

def n59W5 : Finset Nat :=
  ([0, 3, 7, 19, 36, 37, 46, 51, 57, 59] : List Nat).toFinset

def n59W6 : Finset Nat :=
  ([0, 3, 8, 12, 26, 36, 37, 43, 56, 58] : List Nat).toFinset

def n59W7 : Finset Nat :=
  ([0, 5, 14, 16, 29, 37, 49, 55, 56, 59] : List Nat).toFinset

def n59W8 : Finset Nat :=
  ([1, 2, 7, 11, 24, 27, 35, 42, 54, 56] : List Nat).toFinset

def n59W9 : Finset Nat :=
  ([1, 3, 15, 22, 30, 33, 46, 50, 55, 56] : List Nat).toFinset

def n59W10 : Finset Nat :=
  ([1, 3, 16, 22, 23, 33, 47, 51, 56, 59] : List Nat).toFinset

def n59W11 : Finset Nat :=
  ([1, 4, 9, 13, 27, 37, 38, 44, 57, 59] : List Nat).toFinset

def n59W12 : Finset Nat :=
  ([2, 3, 8, 12, 25, 28, 36, 43, 55, 57] : List Nat).toFinset

def n59W13 : Finset Nat :=
  ([2, 4, 16, 23, 31, 34, 47, 51, 56, 57] : List Nat).toFinset

def n59W14 : Finset Nat :=
  ([3, 4, 9, 13, 26, 29, 37, 44, 56, 58] : List Nat).toFinset

def n59W15 : Finset Nat :=
  ([3, 5, 17, 24, 32, 35, 48, 52, 57, 58] : List Nat).toFinset

def n59W16 : Finset Nat :=
  ([4, 5, 10, 14, 27, 30, 38, 45, 57, 59] : List Nat).toFinset

def n59W17 : Finset Nat :=
  ([4, 6, 18, 25, 33, 36, 49, 53, 58, 59] : List Nat).toFinset

def n59Face : Finset (Finset Nat) :=
  ([n59W0, n59W1, n59W2, n59W3, n59W4, n59W5, n59W6, n59W7, n59W8, n59W9, n59W10, n59W11, n59W12, n59W13, n59W14, n59W15, n59W16, n59W17] : List (Finset Nat)).toFinset

def n59FaceIndices : Finset Nat :=
  Finset.range 18

def n59IndexList : List Nat :=
  List.range 18

def n59WitnessOfIndex : Nat -> Finset Nat

  | 0 => n59W0

  | 1 => n59W1

  | 2 => n59W2

  | 3 => n59W3

  | 4 => n59W4

  | 5 => n59W5

  | 6 => n59W6

  | 7 => n59W7

  | 8 => n59W8

  | 9 => n59W9

  | 10 => n59W10

  | 11 => n59W11

  | 12 => n59W12

  | 13 => n59W13

  | 14 => n59W14

  | 15 => n59W15

  | 16 => n59W16

  | 17 => n59W17

  | _ => ∅

theorem n59_indexed_face_matches_exported_face :
    (n59IndexList.map n59WitnessOfIndex).toFinset = n59Face := by
  native_decide

theorem n59_exported_face_card :
    n59Face.card = 18 := by
  native_decide

theorem n59_face_indices_card :
    n59FaceIndices.card = 18 := by
  native_decide

theorem n59_all_exported_witnesses_have_card_h :
    n59W0.card = 10 ∧
    n59W1.card = 10 ∧
    n59W2.card = 10 ∧
    n59W3.card = 10 ∧
    n59W4.card = 10 ∧
    n59W5.card = 10 ∧
    n59W6.card = 10 ∧
    n59W7.card = 10 ∧
    n59W8.card = 10 ∧
    n59W9.card = 10 ∧
    n59W10.card = 10 ∧
    n59W11.card = 10 ∧
    n59W12.card = 10 ∧
    n59W13.card = 10 ∧
    n59W14.card = 10 ∧
    n59W15.card = 10 ∧
    n59W16.card = 10 ∧
    n59W17.card = 10 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_all_exported_witnesses_in_range :
    n59W0 ⊆ Finset.range 60 ∧
    n59W1 ⊆ Finset.range 60 ∧
    n59W2 ⊆ Finset.range 60 ∧
    n59W3 ⊆ Finset.range 60 ∧
    n59W4 ⊆ Finset.range 60 ∧
    n59W5 ⊆ Finset.range 60 ∧
    n59W6 ⊆ Finset.range 60 ∧
    n59W7 ⊆ Finset.range 60 ∧
    n59W8 ⊆ Finset.range 60 ∧
    n59W9 ⊆ Finset.range 60 ∧
    n59W10 ⊆ Finset.range 60 ∧
    n59W11 ⊆ Finset.range 60 ∧
    n59W12 ⊆ Finset.range 60 ∧
    n59W13 ⊆ Finset.range 60 ∧
    n59W14 ⊆ Finset.range 60 ∧
    n59W15 ⊆ Finset.range 60 ∧
    n59W16 ⊆ Finset.range 60 ∧
    n59W17 ⊆ Finset.range 60 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_all_exported_witnesses_are_sidon :
    Erdos.Sidon.IsSidonSet n59W0 ∧
    Erdos.Sidon.IsSidonSet n59W1 ∧
    Erdos.Sidon.IsSidonSet n59W2 ∧
    Erdos.Sidon.IsSidonSet n59W3 ∧
    Erdos.Sidon.IsSidonSet n59W4 ∧
    Erdos.Sidon.IsSidonSet n59W5 ∧
    Erdos.Sidon.IsSidonSet n59W6 ∧
    Erdos.Sidon.IsSidonSet n59W7 ∧
    Erdos.Sidon.IsSidonSet n59W8 ∧
    Erdos.Sidon.IsSidonSet n59W9 ∧
    Erdos.Sidon.IsSidonSet n59W10 ∧
    Erdos.Sidon.IsSidonSet n59W11 ∧
    Erdos.Sidon.IsSidonSet n59W12 ∧
    Erdos.Sidon.IsSidonSet n59W13 ∧
    Erdos.Sidon.IsSidonSet n59W14 ∧
    Erdos.Sidon.IsSidonSet n59W15 ∧
    Erdos.Sidon.IsSidonSet n59W16 ∧
    Erdos.Sidon.IsSidonSet n59W17 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

def n59PrefixRank : Nat -> Nat
  | 0 => 5
  | 1 => 6
  | 2 => 5
  | 3 => 1
  | 4 => 3
  | 5 => 0
  | 6 => 1
  | 7 => 0
  | 8 => 4
  | 9 => 4
  | 10 => 0
  | 11 => 0
  | 12 => 2
  | 13 => 2
  | 14 => 1
  | 15 => 1
  | 16 => 0
  | 17 => 0
  | _ => 99

def n59MassRank : Nat -> Nat
  | 0 => 13
  | 1 => 9
  | 2 => 6
  | 3 => 6
  | 4 => 10
  | 5 => 3
  | 6 => 8
  | 7 => 1
  | 8 => 12
  | 9 => 4
  | 10 => 4
  | 11 => 7
  | 12 => 11
  | 13 => 0
  | 14 => 8
  | 15 => 2
  | 16 => 7
  | 17 => 5
  | _ => 99

def n59JointRank : Nat -> Nat
  | 0 => 15
  | 1 => 14
  | 2 => 10
  | 3 => 6
  | 4 => 12
  | 5 => 1
  | 6 => 9
  | 7 => 0
  | 8 => 13
  | 9 => 8
  | 10 => 2
  | 11 => 7
  | 12 => 11
  | 13 => 5
  | 14 => 9
  | 15 => 3
  | 16 => 7
  | 17 => 4
  | _ => 99

def n59ExactMassTwice (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.densityAdjustedMassTwice 59 (n59WitnessOfIndex i)

def n59MassWinnerIndicesByExactMass : List Nat :=
  n59IndexList.filter (fun i => n59ExactMassTwice i == 7)

def n59PacketPrefixLocation : Nat -> Nat
  | 0 => 55
  | 1 => 23
  | 2 => 55
  | 3 => 58
  | 4 => 10
  | 5 => 59
  | 6 => 58
  | 7 => 59
  | 8 => 56
  | 9 => 56
  | 10 => 59
  | 11 => 59
  | 12 => 57
  | 13 => 57
  | 14 => 58
  | 15 => 58
  | 16 => 59
  | 17 => 59
  | _ => 0

def n59PacketPrefixCount : Nat -> Nat
  | 0 => 10
  | 1 => 6
  | 2 => 10
  | 3 => 10
  | 4 => 4
  | 5 => 10
  | 6 => 10
  | 7 => 10
  | 8 => 10
  | 9 => 10
  | 10 => 10
  | 11 => 10
  | 12 => 10
  | 13 => 10
  | 14 => 10
  | 15 => 10
  | 16 => 10
  | 17 => 10
  | _ => 0

noncomputable def n59ExactPrefixProbe (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.prefixResidualProbeCard10 59 (n59PacketPrefixLocation i) (n59PacketPrefixCount i)

noncomputable def n59FullPrefixResidualMax (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.fullPrefixResidualMax 59 (n59WitnessOfIndex i)


noncomputable def n59ExactJointKey (i : Nat) : ℝ :=
  n59ExactPrefixProbe i * Real.sqrt (59 : ℝ) +
    (n59ExactMassTwice i : ℝ) / 2


noncomputable def n59ScalarFullPrefixJointKey (i : Nat) : ℝ :=
  n59FullPrefixResidualMax i * Real.sqrt (59 : ℝ) +
    (n59ExactMassTwice i : ℝ) / 2


def n59PrefixCountAtPacketLocation (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.prefixCount (n59WitnessOfIndex i) (n59PacketPrefixLocation i)

def n59PrefixWinnerIndicesByRank : List Nat :=
  n59IndexList.filter (fun i => n59PrefixRank i == 0)

def n59MassWinnerIndicesByRank : List Nat :=
  n59IndexList.filter (fun i => n59MassRank i == 0)

def n59JointWinnerIndicesByRank : List Nat :=
  n59IndexList.filter (fun i => n59JointRank i == 0)

def n59DominatesPrefixMass (i j : Nat) : Bool :=
  (n59PrefixRank i <= n59PrefixRank j) &&
  (n59MassRank i <= n59MassRank j) &&
  ((n59PrefixRank i < n59PrefixRank j) ||
   (n59MassRank i < n59MassRank j))

def n59ParetoIndicesByRanks : List Nat :=
  n59IndexList.filter (fun i => !(n59IndexList.any (fun j => n59DominatesPrefixMass j i)))

theorem n59_prefix_winner_match_packet :
    n59PrefixWinnerIndicesByRank = [5, 7, 10, 11, 16, 17] := by
  native_decide

theorem n59_mass_winner_match_packet :
    n59MassWinnerIndicesByRank = [13] := by
  native_decide

theorem n59_joint_winner_match_packet :
    n59JointWinnerIndicesByRank = [7] := by
  native_decide

theorem n59_pareto_minimal_match_packet :
    n59ParetoIndicesByRanks = [7, 13] := by
  native_decide

theorem n59_exact_mass_twice_table :
    n59IndexList.map n59ExactMassTwice = [151, 99, 47, 47, 109, 19, 91, 9, 131, 27, 27, 71, 111, 7, 91, 13, 71, 33] := by
  native_decide

theorem n59_exact_mass_winner_match_packet :
    n59MassWinnerIndicesByExactMass = [13] := by
  native_decide

theorem n59_exact_mass_winner_matches_rank_winner :
    n59MassWinnerIndicesByExactMass = n59MassWinnerIndicesByRank := by
  native_decide

theorem n59_packet_prefix_location_table :
    n59IndexList.map n59PacketPrefixLocation = [55, 23, 55, 58, 10, 59, 58, 59, 56, 56, 59, 59, 57, 57, 58, 58, 59, 59] := by
  native_decide

theorem n59_prefix_count_at_packet_location_table :
    n59IndexList.map n59PrefixCountAtPacketLocation = [10, 6, 10, 10, 4, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10] := by
  native_decide

lemma sqrt59_ge_prefix_lower : (23/3 : ℝ) ≤ Real.sqrt (59 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 23/3) (by norm_num : (0:ℝ) ≤ 59)]
  norm_num

lemma sqrt59_lt_prefix_upper : Real.sqrt (59 : ℝ) < (8 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 8)]
  norm_num

lemma sqrt59_lt_prefix_tight_upper : Real.sqrt (59 : ℝ) < (39/5 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 39/5)]
  norm_num

def n59W5PrefixCountTable : Nat -> Nat
  | 0 => 1
  | 1 => 1
  | 2 => 1
  | 3 => 2
  | 4 => 2
  | 5 => 2
  | 6 => 2
  | 7 => 3
  | 8 => 3
  | 9 => 3
  | 10 => 3
  | 11 => 3
  | 12 => 3
  | 13 => 3
  | 14 => 3
  | 15 => 3
  | 16 => 3
  | 17 => 3
  | 18 => 3
  | 19 => 4
  | 20 => 4
  | 21 => 4
  | 22 => 4
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 4
  | 30 => 4
  | 31 => 4
  | 32 => 4
  | 33 => 4
  | 34 => 4
  | 35 => 4
  | 36 => 5
  | 37 => 6
  | 38 => 6
  | 39 => 6
  | 40 => 6
  | 41 => 6
  | 42 => 6
  | 43 => 6
  | 44 => 6
  | 45 => 6
  | 46 => 7
  | 47 => 7
  | 48 => 7
  | 49 => 7
  | 50 => 7
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 8
  | 57 => 9
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W5_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W5 t) =
      [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

theorem n59W5_prefix_count_model_table :
    (List.range 60).map n59W5PrefixCountTable =
      [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

lemma n59W5_prefix_segment_0_2_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_3_6_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 6 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (6 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_3_6 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_7_18_bound :
    ∀ t : Nat, 7 ≤ t -> t ≤ 18 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (7 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (18 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_7_18 :
    ∀ t : Nat, 7 ≤ t -> t ≤ 18 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_19_35_bound :
    ∀ t : Nat, 19 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (19 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_19_35 :
    ∀ t : Nat, 19 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_36_36_bound :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_36_36 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_37_45_bound :
    ∀ t : Nat, 37 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_37_45 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_46_50_bound :
    ∀ t : Nat, 46 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (50 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_46_50 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_51_56_bound :
    ∀ t : Nat, 51 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (51 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_51_56 :
    ∀ t : Nat, 51 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_57_58_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_57_58 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W5_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W5_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W5 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W5_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 6 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 7 ≤ t -> t ≤ 18 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 19 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 37 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 46 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 51 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W5_prefix_segment_0_2_bound, ⟨n59W5_prefix_segment_3_6_bound, ⟨n59W5_prefix_segment_7_18_bound, ⟨n59W5_prefix_segment_19_35_bound, ⟨n59W5_prefix_segment_36_36_bound, ⟨n59W5_prefix_segment_37_45_bound, ⟨n59W5_prefix_segment_46_50_bound, ⟨n59W5_prefix_segment_51_56_bound, ⟨n59W5_prefix_segment_57_58_bound, n59W5_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W5_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W5 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W5.card = 10 := by
    native_decide
  have hcard : (n59W5.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W5_prefix_residual_segment_0_2_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_0_2_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_3_6_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_3_6_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_3_6 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_7_18_zero :
    ∀ t : Nat, 7 ≤ t -> t ≤ 18 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_7_18_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_7_18 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_19_35_zero :
    ∀ t : Nat, 19 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_19_35_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_19_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_36_36_zero :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_36_36_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_36_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_37_45_zero :
    ∀ t : Nat, 37 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_37_45_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_37_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_46_50_zero :
    ∀ t : Nat, 46 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_46_50_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_46_50 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_51_56_zero :
    ∀ t : Nat, 51 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_51_56_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_51_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_57_58_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_57_58_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_57_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W5_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W5_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W5_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W5 := by
  intro t ht
  interval_cases t
  · exact n59W5_prefix_residual_segment_0_2_zero 0 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_0_2_zero 1 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_0_2_zero 2 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_3_6_zero 3 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_3_6_zero 4 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_3_6_zero 5 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_3_6_zero 6 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 7 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 8 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 9 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 10 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 11 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 12 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 13 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 14 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 15 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 16 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 17 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_7_18_zero 18 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 19 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 20 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 21 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 22 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 23 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 24 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 25 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 26 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 27 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 28 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 29 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 30 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 31 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 32 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 33 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 34 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_19_35_zero 35 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_36_36_zero 36 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 37 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 38 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 39 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 40 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 41 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 42 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 43 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 44 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_37_45_zero 45 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_46_50_zero 46 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_46_50_zero 47 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_46_50_zero 48 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_46_50_zero 49 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_46_50_zero 50 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 51 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 52 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 53 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 54 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 55 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_51_56_zero 56 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_57_58_zero 57 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_57_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W5_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)

def n59W7PrefixCountTable : Nat -> Nat
  | 0 => 1
  | 1 => 1
  | 2 => 1
  | 3 => 1
  | 4 => 1
  | 5 => 2
  | 6 => 2
  | 7 => 2
  | 8 => 2
  | 9 => 2
  | 10 => 2
  | 11 => 2
  | 12 => 2
  | 13 => 2
  | 14 => 3
  | 15 => 3
  | 16 => 4
  | 17 => 4
  | 18 => 4
  | 19 => 4
  | 20 => 4
  | 21 => 4
  | 22 => 4
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 5
  | 30 => 5
  | 31 => 5
  | 32 => 5
  | 33 => 5
  | 34 => 5
  | 35 => 5
  | 36 => 5
  | 37 => 6
  | 38 => 6
  | 39 => 6
  | 40 => 6
  | 41 => 6
  | 42 => 6
  | 43 => 6
  | 44 => 6
  | 45 => 6
  | 46 => 6
  | 47 => 6
  | 48 => 6
  | 49 => 7
  | 50 => 7
  | 51 => 7
  | 52 => 7
  | 53 => 7
  | 54 => 7
  | 55 => 8
  | 56 => 9
  | 57 => 9
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W7_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W7 t) =
      [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 9, 9, 9, 10] := by
  native_decide

theorem n59W7_prefix_count_model_table :
    (List.range 60).map n59W7PrefixCountTable =
      [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 9, 9, 9, 10] := by
  native_decide

lemma n59W7_prefix_segment_0_4_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (4 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_0_4 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_5_13_bound :
    ∀ t : Nat, 5 ≤ t -> t ≤ 13 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (5 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_5_13 :
    ∀ t : Nat, 5 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_14_15_bound :
    ∀ t : Nat, 14 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (15 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_14_15 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_16_28_bound :
    ∀ t : Nat, 16 ≤ t -> t ≤ 28 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (16 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_16_28 :
    ∀ t : Nat, 16 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_29_36_bound :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_29_36 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_37_48_bound :
    ∀ t : Nat, 37 ≤ t -> t ≤ 48 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_37_48 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_49_54_bound :
    ∀ t : Nat, 49 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_49_54 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_55_55_bound :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_55_55 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_56_58_bound :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_56_58 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W7_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W7_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W7 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W7_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 5 ≤ t -> t ≤ 13 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 14 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 16 ≤ t -> t ≤ 28 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 37 ≤ t -> t ≤ 48 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 49 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W7_prefix_segment_0_4_bound, ⟨n59W7_prefix_segment_5_13_bound, ⟨n59W7_prefix_segment_14_15_bound, ⟨n59W7_prefix_segment_16_28_bound, ⟨n59W7_prefix_segment_29_36_bound, ⟨n59W7_prefix_segment_37_48_bound, ⟨n59W7_prefix_segment_49_54_bound, ⟨n59W7_prefix_segment_55_55_bound, ⟨n59W7_prefix_segment_56_58_bound, n59W7_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W7_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W7 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W7.card = 10 := by
    native_decide
  have hcard : (n59W7.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W7_prefix_residual_segment_0_4_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_0_4_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_0_4 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_5_13_zero :
    ∀ t : Nat, 5 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_5_13_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_5_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_14_15_zero :
    ∀ t : Nat, 14 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_14_15_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_14_15 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_16_28_zero :
    ∀ t : Nat, 16 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_16_28_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_16_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_29_36_zero :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_29_36_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_29_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_37_48_zero :
    ∀ t : Nat, 37 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_37_48_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_37_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_49_54_zero :
    ∀ t : Nat, 49 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_49_54_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_49_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_55_55_zero :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_55_55_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_55_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_56_58_zero :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_56_58_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_56_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W7_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W7_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W7_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W7 := by
  intro t ht
  interval_cases t
  · exact n59W7_prefix_residual_segment_0_4_zero 0 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_0_4_zero 1 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_0_4_zero 2 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_0_4_zero 3 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_0_4_zero 4 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 5 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 6 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 7 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 8 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 9 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 10 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 11 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 12 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_5_13_zero 13 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_14_15_zero 14 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_14_15_zero 15 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 16 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 17 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 18 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 19 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 20 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 21 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 22 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 23 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 24 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 25 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 26 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 27 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_16_28_zero 28 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 29 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 30 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 31 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 32 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 33 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 34 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 35 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_29_36_zero 36 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 37 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 38 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 39 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 40 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 41 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 42 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 43 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 44 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 45 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 46 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 47 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_37_48_zero 48 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 49 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 50 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 51 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 52 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 53 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_49_54_zero 54 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_55_55_zero 55 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_56_58_zero 56 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_56_58_zero 57 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_56_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W7_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)

def n59W10PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 1
  | 2 => 1
  | 3 => 2
  | 4 => 2
  | 5 => 2
  | 6 => 2
  | 7 => 2
  | 8 => 2
  | 9 => 2
  | 10 => 2
  | 11 => 2
  | 12 => 2
  | 13 => 2
  | 14 => 2
  | 15 => 2
  | 16 => 3
  | 17 => 3
  | 18 => 3
  | 19 => 3
  | 20 => 3
  | 21 => 3
  | 22 => 4
  | 23 => 5
  | 24 => 5
  | 25 => 5
  | 26 => 5
  | 27 => 5
  | 28 => 5
  | 29 => 5
  | 30 => 5
  | 31 => 5
  | 32 => 5
  | 33 => 6
  | 34 => 6
  | 35 => 6
  | 36 => 6
  | 37 => 6
  | 38 => 6
  | 39 => 6
  | 40 => 6
  | 41 => 6
  | 42 => 6
  | 43 => 6
  | 44 => 6
  | 45 => 6
  | 46 => 6
  | 47 => 7
  | 48 => 7
  | 49 => 7
  | 50 => 7
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 9
  | 57 => 9
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W10_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W10 t) =
      [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10] := by
  native_decide

theorem n59W10_prefix_count_model_table :
    (List.range 60).map n59W10PrefixCountTable =
      [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10] := by
  native_decide

lemma n59W10_prefix_segment_0_0_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_1_2_bound :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_1_2 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_3_15_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (15 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_3_15 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_16_21_bound :
    ∀ t : Nat, 16 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (16 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_16_21 :
    ∀ t : Nat, 16 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_22_22_bound :
    ∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_22_22 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_23_32_bound :
    ∀ t : Nat, 23 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_23_32 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_33_46_bound :
    ∀ t : Nat, 33 ≤ t -> t ≤ 46 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (46 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_33_46 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_47_50_bound :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (47 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (50 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_47_50 :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_51_55_bound :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (51 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_51_55 :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_56_58_bound :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_56_58 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W10_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W10_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W10 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W10_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 16 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 23 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 33 ≤ t -> t ≤ 46 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W10_prefix_segment_0_0_bound, ⟨n59W10_prefix_segment_1_2_bound, ⟨n59W10_prefix_segment_3_15_bound, ⟨n59W10_prefix_segment_16_21_bound, ⟨n59W10_prefix_segment_22_22_bound, ⟨n59W10_prefix_segment_23_32_bound, ⟨n59W10_prefix_segment_33_46_bound, ⟨n59W10_prefix_segment_47_50_bound, ⟨n59W10_prefix_segment_51_55_bound, ⟨n59W10_prefix_segment_56_58_bound, n59W10_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W10_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W10 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W10.card = 10 := by
    native_decide
  have hcard : (n59W10.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W10_prefix_residual_segment_0_0_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_0_0_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_1_2_zero :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_1_2_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_1_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_3_15_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_3_15_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_3_15 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_16_21_zero :
    ∀ t : Nat, 16 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_16_21_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_16_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_22_22_zero :
    ∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_22_22_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_22_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_23_32_zero :
    ∀ t : Nat, 23 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_23_32_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_23_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_33_46_zero :
    ∀ t : Nat, 33 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_33_46_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_33_46 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_47_50_zero :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_47_50_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_47_50 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_51_55_zero :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_51_55_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_51_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_56_58_zero :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_56_58_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_56_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W10_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W10_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W10_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W10 := by
  intro t ht
  interval_cases t
  · exact n59W10_prefix_residual_segment_0_0_zero 0 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_1_2_zero 1 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_1_2_zero 2 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 3 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 4 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 5 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 6 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 7 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 8 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 9 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 10 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 11 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 12 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 13 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 14 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_3_15_zero 15 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 16 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 17 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 18 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 19 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 20 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_16_21_zero 21 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_22_22_zero 22 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 23 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 24 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 25 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 26 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 27 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 28 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 29 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 30 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 31 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_23_32_zero 32 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 33 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 34 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 35 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 36 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 37 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 38 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 39 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 40 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 41 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 42 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 43 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 44 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 45 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_33_46_zero 46 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_47_50_zero 47 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_47_50_zero 48 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_47_50_zero 49 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_47_50_zero 50 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_51_55_zero 51 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_51_55_zero 52 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_51_55_zero 53 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_51_55_zero 54 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_51_55_zero 55 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_56_58_zero 56 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_56_58_zero 57 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_56_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W10_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)

def n59W11PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 1
  | 2 => 1
  | 3 => 1
  | 4 => 2
  | 5 => 2
  | 6 => 2
  | 7 => 2
  | 8 => 2
  | 9 => 3
  | 10 => 3
  | 11 => 3
  | 12 => 3
  | 13 => 4
  | 14 => 4
  | 15 => 4
  | 16 => 4
  | 17 => 4
  | 18 => 4
  | 19 => 4
  | 20 => 4
  | 21 => 4
  | 22 => 4
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 5
  | 28 => 5
  | 29 => 5
  | 30 => 5
  | 31 => 5
  | 32 => 5
  | 33 => 5
  | 34 => 5
  | 35 => 5
  | 36 => 5
  | 37 => 6
  | 38 => 7
  | 39 => 7
  | 40 => 7
  | 41 => 7
  | 42 => 7
  | 43 => 7
  | 44 => 8
  | 45 => 8
  | 46 => 8
  | 47 => 8
  | 48 => 8
  | 49 => 8
  | 50 => 8
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 8
  | 57 => 9
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W11_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W11 t) =
      [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

theorem n59W11_prefix_count_model_table :
    (List.range 60).map n59W11PrefixCountTable =
      [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

lemma n59W11_prefix_segment_0_0_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_1_3_bound :
    ∀ t : Nat, 1 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_1_3 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_4_8_bound :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (8 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_4_8 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_9_12_bound :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (9 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (12 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_9_12 :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_13_26_bound :
    ∀ t : Nat, 13 ≤ t -> t ≤ 26 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (13 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_13_26 :
    ∀ t : Nat, 13 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_27_36_bound :
    ∀ t : Nat, 27 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_27_36 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_37_37_bound :
    ∀ t : Nat, 37 ≤ t -> t ≤ 37 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (37 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_37_37 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 37 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_38_43_bound :
    ∀ t : Nat, 38 ≤ t -> t ≤ 43 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (38 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (43 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_38_43 :
    ∀ t : Nat, 38 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_44_56_bound :
    ∀ t : Nat, 44 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (44 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_44_56 :
    ∀ t : Nat, 44 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_57_58_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_57_58 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W11_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W11_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W11 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W11_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 1 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 13 ≤ t -> t ≤ 26 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 27 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 37 ≤ t -> t ≤ 37 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 38 ≤ t -> t ≤ 43 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 44 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W11_prefix_segment_0_0_bound, ⟨n59W11_prefix_segment_1_3_bound, ⟨n59W11_prefix_segment_4_8_bound, ⟨n59W11_prefix_segment_9_12_bound, ⟨n59W11_prefix_segment_13_26_bound, ⟨n59W11_prefix_segment_27_36_bound, ⟨n59W11_prefix_segment_37_37_bound, ⟨n59W11_prefix_segment_38_43_bound, ⟨n59W11_prefix_segment_44_56_bound, ⟨n59W11_prefix_segment_57_58_bound, n59W11_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W11_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W11 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W11.card = 10 := by
    native_decide
  have hcard : (n59W11.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W11_prefix_residual_segment_0_0_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_0_0_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_1_3_zero :
    ∀ t : Nat, 1 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_1_3_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_1_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_4_8_zero :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_4_8_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_4_8 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_9_12_zero :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_9_12_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_9_12 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_13_26_zero :
    ∀ t : Nat, 13 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_13_26_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_13_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_27_36_zero :
    ∀ t : Nat, 27 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_27_36_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_27_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_37_37_zero :
    ∀ t : Nat, 37 ≤ t -> t ≤ 37 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_37_37_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_37_37 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_38_43_zero :
    ∀ t : Nat, 38 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_38_43_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_38_43 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_44_56_zero :
    ∀ t : Nat, 44 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_44_56_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_44_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_57_58_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_57_58_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_57_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W11_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W11_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W11_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W11 := by
  intro t ht
  interval_cases t
  · exact n59W11_prefix_residual_segment_0_0_zero 0 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_1_3_zero 1 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_1_3_zero 2 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_1_3_zero 3 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_4_8_zero 4 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_4_8_zero 5 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_4_8_zero 6 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_4_8_zero 7 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_4_8_zero 8 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_9_12_zero 9 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_9_12_zero 10 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_9_12_zero 11 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_9_12_zero 12 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 13 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 14 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 15 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 16 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 17 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 18 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 19 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 20 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 21 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 22 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 23 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 24 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 25 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_13_26_zero 26 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 27 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 28 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 29 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 30 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 31 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 32 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 33 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 34 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 35 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_27_36_zero 36 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_37_37_zero 37 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 38 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 39 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 40 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 41 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 42 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_38_43_zero 43 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 44 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 45 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 46 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 47 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 48 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 49 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 50 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 51 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 52 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 53 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 54 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 55 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_44_56_zero 56 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_57_58_zero 57 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_57_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W11_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)

def n59W16PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 0
  | 3 => 0
  | 4 => 1
  | 5 => 2
  | 6 => 2
  | 7 => 2
  | 8 => 2
  | 9 => 2
  | 10 => 3
  | 11 => 3
  | 12 => 3
  | 13 => 3
  | 14 => 4
  | 15 => 4
  | 16 => 4
  | 17 => 4
  | 18 => 4
  | 19 => 4
  | 20 => 4
  | 21 => 4
  | 22 => 4
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 5
  | 28 => 5
  | 29 => 5
  | 30 => 6
  | 31 => 6
  | 32 => 6
  | 33 => 6
  | 34 => 6
  | 35 => 6
  | 36 => 6
  | 37 => 6
  | 38 => 7
  | 39 => 7
  | 40 => 7
  | 41 => 7
  | 42 => 7
  | 43 => 7
  | 44 => 7
  | 45 => 8
  | 46 => 8
  | 47 => 8
  | 48 => 8
  | 49 => 8
  | 50 => 8
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 8
  | 57 => 9
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W16_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W16 t) =
      [0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

theorem n59W16_prefix_count_model_table :
    (List.range 60).map n59W16PrefixCountTable =
      [0, 0, 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

lemma n59W16_prefix_segment_0_3_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_0_3 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_4_4_bound :
    ∀ t : Nat, 4 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (4 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_4_4 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_5_9_bound :
    ∀ t : Nat, 5 ≤ t -> t ≤ 9 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (5 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_5_9 :
    ∀ t : Nat, 5 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_10_13_bound :
    ∀ t : Nat, 10 ≤ t -> t ≤ 13 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_10_13 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_14_26_bound :
    ∀ t : Nat, 14 ≤ t -> t ≤ 26 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_14_26 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_27_29_bound :
    ∀ t : Nat, 27 ≤ t -> t ≤ 29 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_27_29 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_30_37_bound :
    ∀ t : Nat, 30 ≤ t -> t ≤ 37 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (37 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_30_37 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 37 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_38_44_bound :
    ∀ t : Nat, 38 ≤ t -> t ≤ 44 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (38 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_38_44 :
    ∀ t : Nat, 38 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_45_56_bound :
    ∀ t : Nat, 45 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_45_56 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_57_58_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_57_58 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W16_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W16_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W16 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W16_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 4 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 5 ≤ t -> t ≤ 9 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 10 ≤ t -> t ≤ 13 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 14 ≤ t -> t ≤ 26 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 27 ≤ t -> t ≤ 29 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 30 ≤ t -> t ≤ 37 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 38 ≤ t -> t ≤ 44 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 45 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W16_prefix_segment_0_3_bound, ⟨n59W16_prefix_segment_4_4_bound, ⟨n59W16_prefix_segment_5_9_bound, ⟨n59W16_prefix_segment_10_13_bound, ⟨n59W16_prefix_segment_14_26_bound, ⟨n59W16_prefix_segment_27_29_bound, ⟨n59W16_prefix_segment_30_37_bound, ⟨n59W16_prefix_segment_38_44_bound, ⟨n59W16_prefix_segment_45_56_bound, ⟨n59W16_prefix_segment_57_58_bound, n59W16_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W16_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W16 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W16.card = 10 := by
    native_decide
  have hcard : (n59W16.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W16_prefix_residual_segment_0_3_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_0_3_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_0_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_4_4_zero :
    ∀ t : Nat, 4 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_4_4_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_4_4 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_5_9_zero :
    ∀ t : Nat, 5 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_5_9_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_5_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_10_13_zero :
    ∀ t : Nat, 10 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_10_13_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_10_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_14_26_zero :
    ∀ t : Nat, 14 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_14_26_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_14_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_27_29_zero :
    ∀ t : Nat, 27 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_27_29_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_27_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_30_37_zero :
    ∀ t : Nat, 30 ≤ t -> t ≤ 37 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_30_37_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_30_37 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_38_44_zero :
    ∀ t : Nat, 38 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_38_44_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_38_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_45_56_zero :
    ∀ t : Nat, 45 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_45_56_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_45_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_57_58_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_57_58_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_57_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W16_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W16_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W16_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W16 := by
  intro t ht
  interval_cases t
  · exact n59W16_prefix_residual_segment_0_3_zero 0 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_0_3_zero 1 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_0_3_zero 2 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_0_3_zero 3 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_4_4_zero 4 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_5_9_zero 5 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_5_9_zero 6 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_5_9_zero 7 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_5_9_zero 8 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_5_9_zero 9 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_10_13_zero 10 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_10_13_zero 11 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_10_13_zero 12 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_10_13_zero 13 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 14 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 15 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 16 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 17 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 18 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 19 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 20 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 21 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 22 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 23 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 24 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 25 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_14_26_zero 26 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_27_29_zero 27 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_27_29_zero 28 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_27_29_zero 29 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 30 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 31 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 32 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 33 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 34 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 35 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 36 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_30_37_zero 37 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 38 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 39 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 40 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 41 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 42 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 43 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_38_44_zero 44 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 45 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 46 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 47 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 48 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 49 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 50 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 51 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 52 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 53 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 54 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 55 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_45_56_zero 56 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_57_58_zero 57 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_57_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W16_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)

def n59W17PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 0
  | 3 => 0
  | 4 => 1
  | 5 => 1
  | 6 => 2
  | 7 => 2
  | 8 => 2
  | 9 => 2
  | 10 => 2
  | 11 => 2
  | 12 => 2
  | 13 => 2
  | 14 => 2
  | 15 => 2
  | 16 => 2
  | 17 => 2
  | 18 => 3
  | 19 => 3
  | 20 => 3
  | 21 => 3
  | 22 => 3
  | 23 => 3
  | 24 => 3
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 4
  | 30 => 4
  | 31 => 4
  | 32 => 4
  | 33 => 5
  | 34 => 5
  | 35 => 5
  | 36 => 6
  | 37 => 6
  | 38 => 6
  | 39 => 6
  | 40 => 6
  | 41 => 6
  | 42 => 6
  | 43 => 6
  | 44 => 6
  | 45 => 6
  | 46 => 6
  | 47 => 6
  | 48 => 6
  | 49 => 7
  | 50 => 7
  | 51 => 7
  | 52 => 7
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 8
  | 57 => 8
  | 58 => 9
  | 59 => 10
  | _ => 10
theorem n59W17_prefix_count_table_matches_witness :
    (List.range 60).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n59W17 t) =
      [0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

theorem n59W17_prefix_count_model_table :
    (List.range 60).map n59W17PrefixCountTable =
      [0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

lemma n59W17_prefix_segment_0_3_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_0_3 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_4_5_bound :
    ∀ t : Nat, 4 ≤ t -> t ≤ 5 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (5 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_4_5 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_6_17_bound :
    ∀ t : Nat, 6 ≤ t -> t ≤ 17 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (6 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (17 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_6_17 :
    ∀ t : Nat, 6 ≤ t -> t ≤ 17 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_18_24_bound :
    ∀ t : Nat, 18 ≤ t -> t ≤ 24 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (18 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (24 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_18_24 :
    ∀ t : Nat, 18 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_25_32_bound :
    ∀ t : Nat, 25 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (25 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_25_32 :
    ∀ t : Nat, 25 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_33_35_bound :
    ∀ t : Nat, 33 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_33_35 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_36_48_bound :
    ∀ t : Nat, 36 ≤ t -> t ≤ 48 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_36_48 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_49_52_bound :
    ∀ t : Nat, 49 ≤ t -> t ≤ 52 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (52 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_49_52 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_53_57_bound :
    ∀ t : Nat, 53 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (53 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_53_57 :
    ∀ t : Nat, 53 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_58_58_bound :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_58_58 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n59W17_prefix_segment_59_59_bound :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper, hlowR, hhighR]

theorem n59W17_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W17 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W17_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 4 ≤ t -> t ≤ 5 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 6 ≤ t -> t ≤ 17 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 18 ≤ t -> t ≤ 24 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 25 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 33 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 36 ≤ t -> t ≤ 48 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 49 ≤ t -> t ≤ 52 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 53 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) ∧
    (∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        10 * Real.sqrt (59 : ℝ) - 59) := by
  exact ⟨n59W17_prefix_segment_0_3_bound, ⟨n59W17_prefix_segment_4_5_bound, ⟨n59W17_prefix_segment_6_17_bound, ⟨n59W17_prefix_segment_18_24_bound, ⟨n59W17_prefix_segment_25_32_bound, ⟨n59W17_prefix_segment_33_35_bound, ⟨n59W17_prefix_segment_36_48_bound, ⟨n59W17_prefix_segment_49_52_bound, ⟨n59W17_prefix_segment_53_57_bound, ⟨n59W17_prefix_segment_58_58_bound, n59W17_prefix_segment_59_59_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59W17_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W17 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W17.card = 10 := by
    native_decide
  have hcard : (n59W17.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]

theorem n59W17_prefix_residual_segment_0_3_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_0_3_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_0_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_4_5_zero :
    ∀ t : Nat, 4 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_4_5_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_4_5 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_6_17_zero :
    ∀ t : Nat, 6 ≤ t -> t ≤ 17 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_6_17_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_6_17 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_18_24_zero :
    ∀ t : Nat, 18 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_18_24_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_18_24 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_25_32_zero :
    ∀ t : Nat, 25 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_25_32_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_25_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_33_35_zero :
    ∀ t : Nat, 33 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_33_35_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_33_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_36_48_zero :
    ∀ t : Nat, 36 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_36_48_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_36_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_49_52_zero :
    ∀ t : Nat, 49 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_49_52_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_49_52 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_53_57_zero :
    ∀ t : Nat, 53 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_53_57_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_53_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_58_58_zero :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_58_58_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_58_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_prefix_residual_segment_59_59_zero :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 t = 0 := by
  intro t hlow hhigh
  have hdev := n59W17_prefix_segment_59_59_bound t hlow hhigh
  have hcount := n59W17_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n59W17_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W17 := by
  intro t ht
  interval_cases t
  · exact n59W17_prefix_residual_segment_0_3_zero 0 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_0_3_zero 1 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_0_3_zero 2 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_0_3_zero 3 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_4_5_zero 4 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_4_5_zero 5 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 6 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 7 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 8 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 9 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 10 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 11 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 12 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 13 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 14 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 15 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 16 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_6_17_zero 17 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 18 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 19 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 20 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 21 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 22 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 23 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_18_24_zero 24 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 25 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 26 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 27 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 28 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 29 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 30 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 31 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_25_32_zero 32 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_33_35_zero 33 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_33_35_zero 34 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_33_35_zero 35 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 36 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 37 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 38 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 39 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 40 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 41 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 42 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 43 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 44 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 45 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 46 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 47 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_36_48_zero 48 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_49_52_zero 49 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_49_52_zero 50 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_49_52_zero 51 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_49_52_zero 52 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_53_57_zero 53 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_53_57_zero 54 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_53_57_zero 55 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_53_57_zero 56 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_53_57_zero 57 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_58_58_zero 58 (by norm_num) (by norm_num)
  · exact n59W17_prefix_residual_segment_59_59_zero 59 (by norm_num) (by norm_num)


lemma sqrt59_ge_59_div_10 : (59/10 : ℝ) ≤ Real.sqrt (59 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 59/10) (by norm_num : (0:ℝ) ≤ 59)]
  norm_num


theorem n59W0_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W0 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W0.card = 10 := by
    native_decide
  have hcard : (n59W0.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W1_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W1 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W1.card = 10 := by
    native_decide
  have hcard : (n59W1.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W2_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W2 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W2.card = 10 := by
    native_decide
  have hcard : (n59W2.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W3_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W3 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W3.card = 10 := by
    native_decide
  have hcard : (n59W3.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W4_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W4 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W4.card = 10 := by
    native_decide
  have hcard : (n59W4.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W6_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W6 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W6.card = 10 := by
    native_decide
  have hcard : (n59W6.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W8_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W8 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W8.card = 10 := by
    native_decide
  have hcard : (n59W8.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W9_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W9 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W9.card = 10 := by
    native_decide
  have hcard : (n59W9.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W12_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W12 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W12.card = 10 := by
    native_decide
  have hcard : (n59W12.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W13_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W13 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W13.card = 10 := by
    native_decide
  have hcard : (n59W13.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W14_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W14 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W14.card = 10 := by
    native_decide
  have hcard : (n59W14.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59W15_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 59 n59W15 =
      10 * Real.sqrt (59 : ℝ) - 59 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n59W15.card = 10 := by
    native_decide
  have hcard : (n59W15.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (59 : ℝ) := by
    nlinarith [sqrt59_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 59)]
  nlinarith [hsq]


theorem n59_exact_prefix_probe_5_zero :
    n59ExactPrefixProbe 5 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_7_zero :
    n59ExactPrefixProbe 7 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_10_zero :
    n59ExactPrefixProbe 10 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_11_zero :
    n59ExactPrefixProbe 11 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_16_zero :
    n59ExactPrefixProbe 16 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_17_zero :
    n59ExactPrefixProbe 17 = 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n59_exact_prefix_probe_0_positive :
    0 < n59ExactPrefixProbe 0 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_1_positive :
    0 < n59ExactPrefixProbe 1 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (23 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (59 : ℝ) < (9 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 9)]
    norm_num
  have hgap :
      0 < -((23 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) -
        (10 * Real.sqrt (59 : ℝ) - 59) := by
    nlinarith [hsqrtUpper]
  linarith [hgap]


theorem n59_exact_prefix_probe_2_positive :
    0 < n59ExactPrefixProbe 2 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_3_positive :
    0 < n59ExactPrefixProbe 3 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_4_positive :
    0 < n59ExactPrefixProbe 4 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (10 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (59 : ℝ) < (49/6 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 49/6)]
    norm_num
  have hgap :
      0 < -((10 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) -
        (10 * Real.sqrt (59 : ℝ) - 59) := by
    nlinarith [hsqrtUpper]
  linarith [hgap]


theorem n59_exact_prefix_probe_6_positive :
    0 < n59ExactPrefixProbe 6 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_8_positive :
    0 < n59ExactPrefixProbe 8 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_9_positive :
    0 < n59ExactPrefixProbe 9 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_12_positive :
    0 < n59ExactPrefixProbe 12 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_13_positive :
    0 < n59ExactPrefixProbe 13 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_14_positive :
    0 < n59ExactPrefixProbe 14 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_15_positive :
    0 < n59ExactPrefixProbe 15 := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_prefix_probe_positive_on_nonwinners :
    0 < n59ExactPrefixProbe 0 ∧
    0 < n59ExactPrefixProbe 1 ∧
    0 < n59ExactPrefixProbe 2 ∧
    0 < n59ExactPrefixProbe 3 ∧
    0 < n59ExactPrefixProbe 4 ∧
    0 < n59ExactPrefixProbe 6 ∧
    0 < n59ExactPrefixProbe 8 ∧
    0 < n59ExactPrefixProbe 9 ∧
    0 < n59ExactPrefixProbe 12 ∧
    0 < n59ExactPrefixProbe 13 ∧
    0 < n59ExactPrefixProbe 14 ∧
    0 < n59ExactPrefixProbe 15 := by
  exact ⟨n59_exact_prefix_probe_0_positive, ⟨n59_exact_prefix_probe_1_positive, ⟨n59_exact_prefix_probe_2_positive, ⟨n59_exact_prefix_probe_3_positive, ⟨n59_exact_prefix_probe_4_positive, ⟨n59_exact_prefix_probe_6_positive, ⟨n59_exact_prefix_probe_8_positive, ⟨n59_exact_prefix_probe_9_positive, ⟨n59_exact_prefix_probe_12_positive, ⟨n59_exact_prefix_probe_13_positive, ⟨n59_exact_prefix_probe_14_positive, n59_exact_prefix_probe_15_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩


theorem n59_exact_prefix_probe_0_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 0 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 (n59PacketPrefixLocation 0) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W0 (n59PacketPrefixLocation 0) =
        n59PacketPrefixCount 0 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_1_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 1 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 (n59PacketPrefixLocation 1) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W1 (n59PacketPrefixLocation 1) =
        n59PacketPrefixCount 1 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_2_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 2 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 (n59PacketPrefixLocation 2) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W2 (n59PacketPrefixLocation 2) =
        n59PacketPrefixCount 2 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_3_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 3 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 (n59PacketPrefixLocation 3) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W3 (n59PacketPrefixLocation 3) =
        n59PacketPrefixCount 3 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_4_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 4 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 (n59PacketPrefixLocation 4) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W4 (n59PacketPrefixLocation 4) =
        n59PacketPrefixCount 4 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_5_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 5 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W5 (n59PacketPrefixLocation 5) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W5 (n59PacketPrefixLocation 5) =
        n59PacketPrefixCount 5 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_6_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 6 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 (n59PacketPrefixLocation 6) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W6 (n59PacketPrefixLocation 6) =
        n59PacketPrefixCount 6 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_7_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 7 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W7 (n59PacketPrefixLocation 7) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W7 (n59PacketPrefixLocation 7) =
        n59PacketPrefixCount 7 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_8_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 8 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 (n59PacketPrefixLocation 8) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W8 (n59PacketPrefixLocation 8) =
        n59PacketPrefixCount 8 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_9_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 9 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 (n59PacketPrefixLocation 9) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W9 (n59PacketPrefixLocation 9) =
        n59PacketPrefixCount 9 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_10_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 10 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W10 (n59PacketPrefixLocation 10) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W10_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W10 (n59PacketPrefixLocation 10) =
        n59PacketPrefixCount 10 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_11_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 11 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W11 (n59PacketPrefixLocation 11) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W11_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W11 (n59PacketPrefixLocation 11) =
        n59PacketPrefixCount 11 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_12_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 12 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 (n59PacketPrefixLocation 12) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W12 (n59PacketPrefixLocation 12) =
        n59PacketPrefixCount 12 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_13_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 13 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 (n59PacketPrefixLocation 13) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W13 (n59PacketPrefixLocation 13) =
        n59PacketPrefixCount 13 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_14_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 14 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 (n59PacketPrefixLocation 14) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W14 (n59PacketPrefixLocation 14) =
        n59PacketPrefixCount 14 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_15_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 15 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 (n59PacketPrefixLocation 15) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W15 (n59PacketPrefixLocation 15) =
        n59PacketPrefixCount 15 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_16_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 16 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W16 (n59PacketPrefixLocation 16) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W16_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W16 (n59PacketPrefixLocation 16) =
        n59PacketPrefixCount 16 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_exact_prefix_probe_17_matches_prefix_residual_at_packet_location :
    n59ExactPrefixProbe 17 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W17 (n59PacketPrefixLocation 17) := by
  unfold n59ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W17_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n59W17 (n59PacketPrefixLocation 17) =
        n59PacketPrefixCount 17 := by
    native_decide
  rw [hcount]
  rfl


theorem n59_prefix_residual_at_packet_location_0_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 (n59PacketPrefixLocation 0) := by
  rw [← n59_exact_prefix_probe_0_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_0_positive


theorem n59_prefix_residual_at_packet_location_1_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 (n59PacketPrefixLocation 1) := by
  rw [← n59_exact_prefix_probe_1_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_1_positive


theorem n59_prefix_residual_at_packet_location_2_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 (n59PacketPrefixLocation 2) := by
  rw [← n59_exact_prefix_probe_2_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_2_positive


theorem n59_prefix_residual_at_packet_location_3_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 (n59PacketPrefixLocation 3) := by
  rw [← n59_exact_prefix_probe_3_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_3_positive


theorem n59_prefix_residual_at_packet_location_4_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 (n59PacketPrefixLocation 4) := by
  rw [← n59_exact_prefix_probe_4_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_4_positive


theorem n59_prefix_residual_at_packet_location_6_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 (n59PacketPrefixLocation 6) := by
  rw [← n59_exact_prefix_probe_6_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_6_positive


theorem n59_prefix_residual_at_packet_location_8_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 (n59PacketPrefixLocation 8) := by
  rw [← n59_exact_prefix_probe_8_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_8_positive


theorem n59_prefix_residual_at_packet_location_9_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 (n59PacketPrefixLocation 9) := by
  rw [← n59_exact_prefix_probe_9_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_9_positive


theorem n59_prefix_residual_at_packet_location_12_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 (n59PacketPrefixLocation 12) := by
  rw [← n59_exact_prefix_probe_12_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_12_positive


theorem n59_prefix_residual_at_packet_location_13_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 (n59PacketPrefixLocation 13) := by
  rw [← n59_exact_prefix_probe_13_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_13_positive


theorem n59_prefix_residual_at_packet_location_14_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 (n59PacketPrefixLocation 14) := by
  rw [← n59_exact_prefix_probe_14_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_14_positive


theorem n59_prefix_residual_at_packet_location_15_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 (n59PacketPrefixLocation 15) := by
  rw [← n59_exact_prefix_probe_15_matches_prefix_residual_at_packet_location]
  exact n59_exact_prefix_probe_15_positive


theorem n59_full_prefix_residual_zero_on_prefix_winners :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W5 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W7 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W10 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W11 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W16 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W17 := by
  exact ⟨n59W5_full_prefix_residual_zero, ⟨n59W7_full_prefix_residual_zero, ⟨n59W10_full_prefix_residual_zero, ⟨n59W11_full_prefix_residual_zero, ⟨n59W16_full_prefix_residual_zero, n59W17_full_prefix_residual_zero⟩⟩⟩⟩⟩


theorem n59_actual_prefix_residual_positive_on_probe_nonwinners :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 (n59PacketPrefixLocation 0) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 (n59PacketPrefixLocation 1) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 (n59PacketPrefixLocation 2) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 (n59PacketPrefixLocation 3) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 (n59PacketPrefixLocation 4) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 (n59PacketPrefixLocation 6) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 (n59PacketPrefixLocation 8) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 (n59PacketPrefixLocation 9) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 (n59PacketPrefixLocation 12) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 (n59PacketPrefixLocation 13) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 (n59PacketPrefixLocation 14) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 (n59PacketPrefixLocation 15) := by
  exact ⟨n59_prefix_residual_at_packet_location_0_positive, ⟨n59_prefix_residual_at_packet_location_1_positive, ⟨n59_prefix_residual_at_packet_location_2_positive, ⟨n59_prefix_residual_at_packet_location_3_positive, ⟨n59_prefix_residual_at_packet_location_4_positive, ⟨n59_prefix_residual_at_packet_location_6_positive, ⟨n59_prefix_residual_at_packet_location_8_positive, ⟨n59_prefix_residual_at_packet_location_9_positive, ⟨n59_prefix_residual_at_packet_location_12_positive, ⟨n59_prefix_residual_at_packet_location_13_positive, ⟨n59_prefix_residual_at_packet_location_14_positive, n59_prefix_residual_at_packet_location_15_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩


theorem n59_full_prefix_residual_max_5_zero :
    n59FullPrefixResidualMax 5 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W5_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_7_zero :
    n59FullPrefixResidualMax 7 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W7_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_10_zero :
    n59FullPrefixResidualMax 10 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W10_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_11_zero :
    n59FullPrefixResidualMax 11 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W11_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_16_zero :
    n59FullPrefixResidualMax 16 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W16_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_17_zero :
    n59FullPrefixResidualMax 17 = 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n59W17_full_prefix_residual_zero


theorem n59_full_prefix_residual_max_0_positive :
    0 < n59FullPrefixResidualMax 0 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 0 ≤ 59)
    n59_prefix_residual_at_packet_location_0_positive


theorem n59_full_prefix_residual_max_1_positive :
    0 < n59FullPrefixResidualMax 1 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 1 ≤ 59)
    n59_prefix_residual_at_packet_location_1_positive


theorem n59_full_prefix_residual_max_2_positive :
    0 < n59FullPrefixResidualMax 2 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 2 ≤ 59)
    n59_prefix_residual_at_packet_location_2_positive


theorem n59_full_prefix_residual_max_3_positive :
    0 < n59FullPrefixResidualMax 3 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 3 ≤ 59)
    n59_prefix_residual_at_packet_location_3_positive


theorem n59_full_prefix_residual_max_4_positive :
    0 < n59FullPrefixResidualMax 4 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 4 ≤ 59)
    n59_prefix_residual_at_packet_location_4_positive


theorem n59_full_prefix_residual_max_6_positive :
    0 < n59FullPrefixResidualMax 6 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 6 ≤ 59)
    n59_prefix_residual_at_packet_location_6_positive


theorem n59_full_prefix_residual_max_8_positive :
    0 < n59FullPrefixResidualMax 8 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 8 ≤ 59)
    n59_prefix_residual_at_packet_location_8_positive


theorem n59_full_prefix_residual_max_9_positive :
    0 < n59FullPrefixResidualMax 9 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 9 ≤ 59)
    n59_prefix_residual_at_packet_location_9_positive


theorem n59_full_prefix_residual_max_12_positive :
    0 < n59FullPrefixResidualMax 12 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 12 ≤ 59)
    n59_prefix_residual_at_packet_location_12_positive


theorem n59_full_prefix_residual_max_13_positive :
    0 < n59FullPrefixResidualMax 13 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 13 ≤ 59)
    n59_prefix_residual_at_packet_location_13_positive


theorem n59_full_prefix_residual_max_14_positive :
    0 < n59FullPrefixResidualMax 14 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 14 ≤ 59)
    n59_prefix_residual_at_packet_location_14_positive


theorem n59_full_prefix_residual_max_15_positive :
    0 < n59FullPrefixResidualMax 15 := by
  unfold n59FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n59PacketPrefixLocation 15 ≤ 59)
    n59_prefix_residual_at_packet_location_15_positive


theorem n59_full_prefix_residual_max_zero_on_prefix_winners :
    n59FullPrefixResidualMax 5 = 0 ∧
        n59FullPrefixResidualMax 7 = 0 ∧
        n59FullPrefixResidualMax 10 = 0 ∧
        n59FullPrefixResidualMax 11 = 0 ∧
        n59FullPrefixResidualMax 16 = 0 ∧
        n59FullPrefixResidualMax 17 = 0 := by
  exact ⟨n59_full_prefix_residual_max_5_zero, ⟨n59_full_prefix_residual_max_7_zero, ⟨n59_full_prefix_residual_max_10_zero, ⟨n59_full_prefix_residual_max_11_zero, ⟨n59_full_prefix_residual_max_16_zero, n59_full_prefix_residual_max_17_zero⟩⟩⟩⟩⟩


theorem n59_full_prefix_residual_max_positive_on_probe_nonwinners :
    0 < n59FullPrefixResidualMax 0 ∧
        0 < n59FullPrefixResidualMax 1 ∧
        0 < n59FullPrefixResidualMax 2 ∧
        0 < n59FullPrefixResidualMax 3 ∧
        0 < n59FullPrefixResidualMax 4 ∧
        0 < n59FullPrefixResidualMax 6 ∧
        0 < n59FullPrefixResidualMax 8 ∧
        0 < n59FullPrefixResidualMax 9 ∧
        0 < n59FullPrefixResidualMax 12 ∧
        0 < n59FullPrefixResidualMax 13 ∧
        0 < n59FullPrefixResidualMax 14 ∧
        0 < n59FullPrefixResidualMax 15 := by
  exact ⟨n59_full_prefix_residual_max_0_positive, ⟨n59_full_prefix_residual_max_1_positive, ⟨n59_full_prefix_residual_max_2_positive, ⟨n59_full_prefix_residual_max_3_positive, ⟨n59_full_prefix_residual_max_4_positive, ⟨n59_full_prefix_residual_max_6_positive, ⟨n59_full_prefix_residual_max_8_positive, ⟨n59_full_prefix_residual_max_9_positive, ⟨n59_full_prefix_residual_max_12_positive, ⟨n59_full_prefix_residual_max_13_positive, ⟨n59_full_prefix_residual_max_14_positive, n59_full_prefix_residual_max_15_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩


theorem n59_exact_prefix_probe_0_value :
    n59ExactPrefixProbe 0 = (4 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_0_value :
    n59ExactMassTwice 0 = 151 := by
  native_decide


theorem n59_exact_prefix_probe_1_value :
    n59ExactPrefixProbe 1 = ((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (23 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (59 : ℝ) < (9 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 9)]
    norm_num
  have hgap :
      0 < -((23 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) -
        (10 * Real.sqrt (59 : ℝ) - 59) := by
    nlinarith [hsqrtUpper]
  rw [max_eq_right (le_of_lt hgap)]
  ring


theorem n59_exact_mass_twice_1_value :
    n59ExactMassTwice 1 = 99 := by
  native_decide


theorem n59_exact_prefix_probe_2_value :
    n59ExactPrefixProbe 2 = (4 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_2_value :
    n59ExactMassTwice 2 = 47 := by
  native_decide


theorem n59_exact_prefix_probe_3_value :
    n59ExactPrefixProbe 3 = (1 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_3_value :
    n59ExactMassTwice 3 = 47 := by
  native_decide


theorem n59_exact_prefix_probe_4_value :
    n59ExactPrefixProbe 4 = ((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (10 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (59 : ℝ) < (49/6 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 49/6)]
    norm_num
  have hgap :
      0 < -((10 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) -
        (10 * Real.sqrt (59 : ℝ) - 59) := by
    nlinarith [hsqrtUpper]
  rw [max_eq_right (le_of_lt hgap)]
  ring


theorem n59_exact_mass_twice_4_value :
    n59ExactMassTwice 4 = 109 := by
  native_decide


theorem n59_exact_prefix_probe_5_value :
    n59ExactPrefixProbe 5 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_5_value :
    n59ExactMassTwice 5 = 19 := by
  native_decide


theorem n59_exact_prefix_probe_6_value :
    n59ExactPrefixProbe 6 = (1 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_6_value :
    n59ExactMassTwice 6 = 91 := by
  native_decide


theorem n59_exact_prefix_probe_7_value :
    n59ExactPrefixProbe 7 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_7_value :
    n59ExactMassTwice 7 = 9 := by
  native_decide


theorem n59_exact_prefix_probe_8_value :
    n59ExactPrefixProbe 8 = (3 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_8_value :
    n59ExactMassTwice 8 = 131 := by
  native_decide


theorem n59_exact_prefix_probe_9_value :
    n59ExactPrefixProbe 9 = (3 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_9_value :
    n59ExactMassTwice 9 = 27 := by
  native_decide


theorem n59_exact_prefix_probe_10_value :
    n59ExactPrefixProbe 10 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_10_value :
    n59ExactMassTwice 10 = 27 := by
  native_decide


theorem n59_exact_prefix_probe_11_value :
    n59ExactPrefixProbe 11 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_11_value :
    n59ExactMassTwice 11 = 71 := by
  native_decide


theorem n59_exact_prefix_probe_12_value :
    n59ExactPrefixProbe 12 = (2 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_12_value :
    n59ExactMassTwice 12 = 111 := by
  native_decide


theorem n59_exact_prefix_probe_13_value :
    n59ExactPrefixProbe 13 = (2 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_13_value :
    n59ExactMassTwice 13 = 7 := by
  native_decide


theorem n59_exact_prefix_probe_14_value :
    n59ExactPrefixProbe 14 = (1 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_14_value :
    n59ExactMassTwice 14 = 91 := by
  native_decide


theorem n59_exact_prefix_probe_15_value :
    n59ExactPrefixProbe 15 = (1 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_15_value :
    n59ExactMassTwice 15 = 13 := by
  native_decide


theorem n59_exact_prefix_probe_16_value :
    n59ExactPrefixProbe 16 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_16_value :
    n59ExactMassTwice 16 = 71 := by
  native_decide


theorem n59_exact_prefix_probe_17_value :
    n59ExactPrefixProbe 17 = (0 : ℝ) := by
  simp [n59ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n59PacketPrefixLocation, n59PacketPrefixCount]
  have hnonpos : (59 : ℝ) - 10 * Real.sqrt (59 : ℝ) ≤ 0 := by
    nlinarith [sqrt59_ge_59_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n59_exact_mass_twice_17_value :
    n59ExactMassTwice 17 = 33 := by
  native_decide


theorem n59W0_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_1_5 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_1_5_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (5 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_1_5 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_6_9 :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_6_9_le_exact_probe :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (6 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_6_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_10_22 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_10_22_le_exact_probe :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_10_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_23_25 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_23_25_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_23_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_26_33 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_26_33_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_26_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_34_40 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_34_40_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (40 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_34_40 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_41_52 :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_41_52_le_exact_probe :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (41 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (52 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_41_52 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_53_54 :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_53_54_le_exact_probe :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (53 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_53_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_prefix_count_segment_55_59 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W0 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W0_prefix_residual_segment_55_59_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W0_prefix_count_segment_55_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_0_positive
  · rw [n59_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W0_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 t ≤
        n59ExactPrefixProbe 0 := by
  intro t ht
  interval_cases t
  · exact n59W0_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_1_5_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_1_5_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_1_5_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_1_5_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_1_5_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_6_9_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_6_9_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_6_9_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_6_9_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_10_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_23_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_23_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_23_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_26_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_34_40_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_41_52_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_53_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_53_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_55_59_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_55_59_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_55_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_55_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W0_prefix_residual_segment_55_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W1_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_2_7 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_2_7_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_2_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_8_12 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_8_12_le_exact_probe :
    ∀ t : Nat, 8 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (12 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_8_12 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_13_21 :
    ∀ t : Nat, 13 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_13_21_le_exact_probe :
    ∀ t : Nat, 13 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (13 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_13_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_22_22 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_22_22_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_22_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_23_39 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 39 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_23_39_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 39 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (39 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_23_39 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_40_51 :
    ∀ t : Nat, 40 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_40_51_le_exact_probe :
    ∀ t : Nat, 40 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (40 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (51 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_40_51 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_52_55 :
    ∀ t : Nat, 52 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_52_55_le_exact_probe :
    ∀ t : Nat, 52 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (52 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_52_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_56_58 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_56_58_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_56_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W1 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W1_prefix_residual_segment_59_59_le_exact_probe :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((36 : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W1_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_1_positive
  · rw [n59_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W1_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 t ≤
        n59ExactPrefixProbe 1 := by
  intro t ht
  interval_cases t
  · exact n59W1_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_2_7_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_8_12_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_8_12_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_8_12_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_8_12_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_8_12_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_13_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_22_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_23_39_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_40_51_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_52_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_52_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_52_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_52_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_56_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_56_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_56_58_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W1_prefix_residual_segment_59_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W2_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_2_13 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_2_13_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_2_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_14_20 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_14_20_le_exact_probe :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_14_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_21_28 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_21_28_le_exact_probe :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_21_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_29_31 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_29_31_le_exact_probe :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_29_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_32_44 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_32_44_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_32_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_45_48 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_45_48_le_exact_probe :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_45_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_49_53 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_49_53_le_exact_probe :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_49_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_54_54 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_54_54_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_54_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_prefix_count_segment_55_59 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W2 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W2_prefix_residual_segment_55_59_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (4 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W2_prefix_count_segment_55_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_2_positive
  · rw [n59_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W2_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 t ≤
        n59ExactPrefixProbe 2 := by
  intro t ht
  interval_cases t
  · exact n59W2_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_2_13_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_14_20_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_21_28_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_29_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_29_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_29_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_32_44_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_45_48_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_45_48_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_45_48_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_45_48_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_49_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_49_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_49_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_49_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_49_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_54_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_55_59_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_55_59_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_55_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_55_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W2_prefix_residual_segment_55_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W3_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_2_14 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_2_14_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_2_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_15_20 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_15_20_le_exact_probe :
    ∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_15_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_21_21 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_21_21_le_exact_probe :
    ∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_21_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_22_31 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_22_31_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_22_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_32_45 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_32_45_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_32_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_46_49_le_exact_probe :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_50_54_le_exact_probe :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_55_57 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_55_57_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_55_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_prefix_count_segment_58_59 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W3 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W3_prefix_residual_segment_58_59_le_exact_probe :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W3_prefix_count_segment_58_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_3_positive
  · rw [n59_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W3_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 t ≤
        n59ExactPrefixProbe 3 := by
  intro t ht
  interval_cases t
  · exact n59W3_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_2_14_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_15_20_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_21_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_22_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_32_45_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_46_49_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_46_49_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_46_49_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_46_49_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_50_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_50_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_50_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_50_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_50_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_55_57_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_55_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_55_57_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_58_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W3_prefix_residual_segment_58_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W4_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_0_2_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_3_3 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_3_3_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_3_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_4_9 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_4_9_le_exact_probe :
    ∀ t : Nat, 4 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_4_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_10_21 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_10_21_le_exact_probe :
    ∀ t : Nat, 10 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_10_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_22_29 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_22_29_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_22_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_30_42 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_30_42_le_exact_probe :
    ∀ t : Nat, 30 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_30_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_43_44 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_43_44_le_exact_probe :
    ∀ t : Nat, 43 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_43_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_45_53 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_45_53_le_exact_probe :
    ∀ t : Nat, 45 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_45_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_54_58 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_54_58_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_54_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_prefix_count_segment_59_59 :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W4 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W4_prefix_residual_segment_59_59_le_exact_probe :
    ∀ t : Nat, 59 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (59 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (((49 : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W4_prefix_count_segment_59_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_4_positive
  · rw [n59_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W4_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 t ≤
        n59ExactPrefixProbe 4 := by
  intro t ht
  interval_cases t
  · exact n59W4_prefix_residual_segment_0_2_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_0_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_0_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_3_3_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_4_9_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_10_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_22_29_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_30_42_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_43_44_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_43_44_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_45_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_54_58_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_54_58_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_54_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_54_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_54_58_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W4_prefix_residual_segment_59_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W6_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_0_2_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_3_7 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_3_7_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_3_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_8_11 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_8_11_le_exact_probe :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (11 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_8_11 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_12_25 :
    ∀ t : Nat, 12 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_12_25_le_exact_probe :
    ∀ t : Nat, 12 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (12 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_12_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_26_35 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_26_35_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_26_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_36_36 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_36_36_le_exact_probe :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_36_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_37_42 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_37_42_le_exact_probe :
    ∀ t : Nat, 37 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_37_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_43_55 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_43_55_le_exact_probe :
    ∀ t : Nat, 43 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_43_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_56_57_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_prefix_count_segment_58_59 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W6 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W6_prefix_residual_segment_58_59_le_exact_probe :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W6_prefix_count_segment_58_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_6_positive
  · rw [n59_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W6_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 t ≤
        n59ExactPrefixProbe 6 := by
  intro t ht
  interval_cases t
  · exact n59W6_prefix_residual_segment_0_2_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_0_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_0_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_3_7_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_3_7_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_3_7_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_3_7_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_3_7_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_8_11_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_8_11_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_8_11_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_8_11_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_12_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_26_35_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_36_36_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_37_42_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_43_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_56_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_56_57_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_58_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W6_prefix_residual_segment_58_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W8_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_1_1 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_1_1_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_1_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_2_6 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_2_6_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (6 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_2_6 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_7_10 :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_7_10_le_exact_probe :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (7 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (10 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_7_10 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_11_23 :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_11_23_le_exact_probe :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (11 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_11_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_24_26 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_24_26_le_exact_probe :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_24_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_27_34 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_27_34_le_exact_probe :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_27_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_35_41 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_35_41_le_exact_probe :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (41 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_35_41 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_42_53 :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_42_53_le_exact_probe :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (42 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_42_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_54_55 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_54_55_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_54_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_prefix_count_segment_56_59 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W8 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W8_prefix_residual_segment_56_59_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W8_prefix_count_segment_56_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_8_positive
  · rw [n59_exact_prefix_probe_8_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W8_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 t ≤
        n59ExactPrefixProbe 8 := by
  intro t ht
  interval_cases t
  · exact n59W8_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_1_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_2_6_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_2_6_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_2_6_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_2_6_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_2_6_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_7_10_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_7_10_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_7_10_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_7_10_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_11_23_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_24_26_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_24_26_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_24_26_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_27_34_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_35_41_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_42_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_54_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_54_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_56_59_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_56_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_56_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W8_prefix_residual_segment_56_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W9_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_1_2 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_1_2_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_1_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_3_14 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_3_14_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_3_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_15_21 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_15_21_le_exact_probe :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_15_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_22_29 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_22_29_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_22_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_30_32 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_30_32_le_exact_probe :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_30_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_33_45 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_33_45_le_exact_probe :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_33_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_46_49_le_exact_probe :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_50_54_le_exact_probe :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_55_55 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_55_55_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_55_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_prefix_count_segment_56_59 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W9 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W9_prefix_residual_segment_56_59_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W9_prefix_count_segment_56_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_9_positive
  · rw [n59_exact_prefix_probe_9_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W9_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 t ≤
        n59ExactPrefixProbe 9 := by
  intro t ht
  interval_cases t
  · exact n59W9_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_1_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_1_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_3_14_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_15_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_22_29_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_30_32_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_30_32_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_30_32_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_33_45_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_46_49_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_46_49_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_46_49_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_46_49_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_50_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_50_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_50_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_50_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_50_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_55_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_56_59_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_56_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_56_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W9_prefix_residual_segment_56_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W12_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_2_2 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_2_2_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_2_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_3_7 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_3_7_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_3_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_8_11 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_8_11_le_exact_probe :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (11 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_8_11 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_12_24 :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_12_24_le_exact_probe :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (12 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (24 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_12_24 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_25_27 :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_25_27_le_exact_probe :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (25 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (27 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_25_27 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_28_35 :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_28_35_le_exact_probe :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (28 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_28_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_36_42 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_36_42_le_exact_probe :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_36_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_43_54 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_43_54_le_exact_probe :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_43_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_55_56 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_55_56_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_55_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_prefix_count_segment_57_59 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W12 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W12_prefix_residual_segment_57_59_le_exact_probe :
    ∀ t : Nat, 57 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W12_prefix_count_segment_57_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W12_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_12_positive
  · rw [n59_exact_prefix_probe_12_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W12_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 t ≤
        n59ExactPrefixProbe 12 := by
  intro t ht
  interval_cases t
  · exact n59W12_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_2_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_3_7_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_3_7_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_3_7_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_3_7_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_3_7_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_8_11_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_8_11_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_8_11_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_8_11_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_12_24_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_25_27_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_25_27_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_25_27_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_28_35_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_36_42_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_43_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_55_56_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_55_56_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_57_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_57_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W12_prefix_residual_segment_57_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W13_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_2_3 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_2_3_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_2_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_4_15 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_4_15_le_exact_probe :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (15 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_4_15 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_16_22 :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_16_22_le_exact_probe :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (16 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_16_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_23_30 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_23_30_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (30 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_23_30 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_31_33 :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_31_33_le_exact_probe :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (31 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_31_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_34_46 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_34_46_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (46 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_34_46 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_47_50 :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_47_50_le_exact_probe :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (47 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (50 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_47_50 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_51_55 :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_51_55_le_exact_probe :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (51 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_51_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_56_56 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_56_56_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_56_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_prefix_count_segment_57_59 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W13 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W13_prefix_residual_segment_57_59_le_exact_probe :
    ∀ t : Nat, 57 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W13_prefix_count_segment_57_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W13_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_13_positive
  · rw [n59_exact_prefix_probe_13_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W13_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 t ≤
        n59ExactPrefixProbe 13 := by
  intro t ht
  interval_cases t
  · exact n59W13_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_2_3_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_2_3_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_4_15_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_16_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_23_30_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_31_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_31_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_31_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_34_46_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_47_50_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_47_50_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_47_50_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_47_50_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_51_55_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_51_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_51_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_51_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_51_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_56_56_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_57_59_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_57_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W13_prefix_residual_segment_57_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W14_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_0_2_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_3_3 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_3_3_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_3_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_4_8 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_4_8_le_exact_probe :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (8 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_4_8 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_9_12 :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_9_12_le_exact_probe :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (9 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (12 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_9_12 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_13_25 :
    ∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_13_25_le_exact_probe :
    ∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (13 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_13_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_26_28 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_26_28_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_26_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_29_36 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_29_36_le_exact_probe :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_29_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_37_43 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_37_43_le_exact_probe :
    ∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (43 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_37_43 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_44_55 :
    ∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_44_55_le_exact_probe :
    ∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (44 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_44_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_56_57_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_prefix_count_segment_58_59 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W14 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W14_prefix_residual_segment_58_59_le_exact_probe :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W14_prefix_count_segment_58_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W14_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_14_positive
  · rw [n59_exact_prefix_probe_14_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W14_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 t ≤
        n59ExactPrefixProbe 14 := by
  intro t ht
  interval_cases t
  · exact n59W14_prefix_residual_segment_0_2_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_0_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_0_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_3_3_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_4_8_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_4_8_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_4_8_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_4_8_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_4_8_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_9_12_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_9_12_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_9_12_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_9_12_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_13_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_26_28_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_26_28_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_26_28_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_29_36_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_37_43_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_44_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_56_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_56_57_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_58_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W14_prefix_residual_segment_58_59_le_exact_probe 59 (by norm_num) (by norm_num)

theorem n59W15_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_0_2_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_3_4 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_3_4_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (4 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_3_4 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_5_16 :
    ∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_5_16_le_exact_probe :
    ∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (5 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (16 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_5_16 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_17_23 :
    ∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_17_23_le_exact_probe :
    ∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (17 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_17_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_24_31 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_24_31_le_exact_probe :
    ∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_24_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_32_34 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_32_34_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_32_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_35_47 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_35_47_le_exact_probe :
    ∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (47 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_35_47 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_48_51 :
    ∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_48_51_le_exact_probe :
    ∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (48 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (51 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_48_51 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_52_56 :
    ∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_52_56_le_exact_probe :
    ∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (52 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_52_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_57_57 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_57_57_le_exact_probe :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_57_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_prefix_count_segment_58_59 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixCount n59W15 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n59W15_prefix_residual_segment_58_59_le_exact_probe :
    ∀ t : Nat, 58 ≤ t -> t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (59 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (59 : ℝ)| ≤
        (10 * Real.sqrt (59 : ℝ) - 59) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper,
        sqrt59_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n59W15_prefix_count_segment_58_59 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n59W15_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n59_exact_prefix_probe_15_positive
  · rw [n59_exact_prefix_probe_15_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n59W15_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 59 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 t ≤
        n59ExactPrefixProbe 15 := by
  intro t ht
  interval_cases t
  · exact n59W15_prefix_residual_segment_0_2_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_0_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_0_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_3_4_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_3_4_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_5_16_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_17_23_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_24_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_32_34_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_32_34_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_32_34_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_35_47_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_48_51_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_48_51_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_48_51_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_48_51_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_52_56_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_52_56_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_52_56_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_52_56_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_52_56_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_57_57_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_58_59_le_exact_probe 58 (by norm_num) (by norm_num)
  · exact n59W15_prefix_residual_segment_58_59_le_exact_probe 59 (by norm_num) (by norm_num)


theorem n59_full_prefix_residual_max_0_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 0 = n59ExactPrefixProbe 0 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 0 ≤ 59)
    n59W0_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_0_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_1_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 1 = n59ExactPrefixProbe 1 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 1 ≤ 59)
    n59W1_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_1_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_2_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 2 = n59ExactPrefixProbe 2 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 2 ≤ 59)
    n59W2_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_2_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_3_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 3 = n59ExactPrefixProbe 3 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 3 ≤ 59)
    n59W3_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_3_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_4_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 4 = n59ExactPrefixProbe 4 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 4 ≤ 59)
    n59W4_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_4_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_5_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 5 = n59ExactPrefixProbe 5 := by
  rw [n59_full_prefix_residual_max_5_zero, n59_exact_prefix_probe_5_zero]


theorem n59_full_prefix_residual_max_6_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 6 = n59ExactPrefixProbe 6 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 6 ≤ 59)
    n59W6_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_6_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_7_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 7 = n59ExactPrefixProbe 7 := by
  rw [n59_full_prefix_residual_max_7_zero, n59_exact_prefix_probe_7_zero]


theorem n59_full_prefix_residual_max_8_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 8 = n59ExactPrefixProbe 8 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 8 ≤ 59)
    n59W8_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_8_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_9_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 9 = n59ExactPrefixProbe 9 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 9 ≤ 59)
    n59W9_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_9_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_10_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 10 = n59ExactPrefixProbe 10 := by
  rw [n59_full_prefix_residual_max_10_zero, n59_exact_prefix_probe_10_zero]


theorem n59_full_prefix_residual_max_11_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 11 = n59ExactPrefixProbe 11 := by
  rw [n59_full_prefix_residual_max_11_zero, n59_exact_prefix_probe_11_zero]


theorem n59_full_prefix_residual_max_12_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 12 = n59ExactPrefixProbe 12 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 12 ≤ 59)
    n59W12_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_12_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_13_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 13 = n59ExactPrefixProbe 13 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 13 ≤ 59)
    n59W13_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_13_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_14_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 14 = n59ExactPrefixProbe 14 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 14 ≤ 59)
    n59W14_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_14_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_15_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 15 = n59ExactPrefixProbe 15 := by
  unfold n59FullPrefixResidualMax n59WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n59PacketPrefixLocation 15 ≤ 59)
    n59W15_full_prefix_residual_le_exact_probe
    (n59_exact_prefix_probe_15_matches_prefix_residual_at_packet_location).symm


theorem n59_full_prefix_residual_max_16_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 16 = n59ExactPrefixProbe 16 := by
  rw [n59_full_prefix_residual_max_16_zero, n59_exact_prefix_probe_16_zero]


theorem n59_full_prefix_residual_max_17_matches_exact_prefix_probe :
    n59FullPrefixResidualMax 17 = n59ExactPrefixProbe 17 := by
  rw [n59_full_prefix_residual_max_17_zero, n59_exact_prefix_probe_17_zero]


theorem n59_full_prefix_residual_max_matches_exact_prefix_probe_on_face :
    ∀ i ∈ n59FaceIndices, n59FullPrefixResidualMax i = n59ExactPrefixProbe i := by
  intro i hi
  simp [n59FaceIndices] at hi
  interval_cases i
  · exact n59_full_prefix_residual_max_0_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_1_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_2_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_3_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_4_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_5_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_6_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_7_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_8_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_9_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_10_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_11_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_12_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_13_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_14_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_15_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_16_matches_exact_prefix_probe
  · exact n59_full_prefix_residual_max_17_matches_exact_prefix_probe


theorem n59_exact_joint_key_0_value :
    n59ExactJointKey 0 = (4 : ℝ) * Real.sqrt (59 : ℝ) + (151/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_0_value, n59_exact_mass_twice_0_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_0_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 0 = n59ExactJointKey 0 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_0_matches_exact_prefix_probe]


theorem n59_exact_joint_key_1_value :
    n59ExactJointKey 1 = (36 : ℝ) * Real.sqrt (59 : ℝ) - (373/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_1_value, n59_exact_mass_twice_1_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_1_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 1 = n59ExactJointKey 1 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_1_matches_exact_prefix_probe]


theorem n59_exact_joint_key_2_value :
    n59ExactJointKey 2 = (4 : ℝ) * Real.sqrt (59 : ℝ) + (47/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_2_value, n59_exact_mass_twice_2_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_2_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 2 = n59ExactJointKey 2 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_2_matches_exact_prefix_probe]


theorem n59_exact_joint_key_3_value :
    n59ExactJointKey 3 = (1 : ℝ) * Real.sqrt (59 : ℝ) + (47/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_3_value, n59_exact_mass_twice_3_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_3_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 3 = n59ExactJointKey 3 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_3_matches_exact_prefix_probe]


theorem n59_exact_joint_key_4_value :
    n59ExactJointKey 4 = (49 : ℝ) * Real.sqrt (59 : ℝ) - (599/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_4_value, n59_exact_mass_twice_4_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_4_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 4 = n59ExactJointKey 4 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_4_matches_exact_prefix_probe]


theorem n59_exact_joint_key_5_value :
    n59ExactJointKey 5 = (19/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_5_value, n59_exact_mass_twice_5_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_5_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 5 = n59ExactJointKey 5 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_5_matches_exact_prefix_probe]


theorem n59_exact_joint_key_6_value :
    n59ExactJointKey 6 = (1 : ℝ) * Real.sqrt (59 : ℝ) + (91/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_6_value, n59_exact_mass_twice_6_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_6_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 6 = n59ExactJointKey 6 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_6_matches_exact_prefix_probe]


theorem n59_exact_joint_key_7_value :
    n59ExactJointKey 7 = (9/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_7_value, n59_exact_mass_twice_7_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 7 = n59ExactJointKey 7 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_7_matches_exact_prefix_probe]


theorem n59_exact_joint_key_8_value :
    n59ExactJointKey 8 = (3 : ℝ) * Real.sqrt (59 : ℝ) + (131/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_8_value, n59_exact_mass_twice_8_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_8_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 8 = n59ExactJointKey 8 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_8_matches_exact_prefix_probe]


theorem n59_exact_joint_key_9_value :
    n59ExactJointKey 9 = (3 : ℝ) * Real.sqrt (59 : ℝ) + (27/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_9_value, n59_exact_mass_twice_9_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_9_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 9 = n59ExactJointKey 9 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_9_matches_exact_prefix_probe]


theorem n59_exact_joint_key_10_value :
    n59ExactJointKey 10 = (27/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_10_value, n59_exact_mass_twice_10_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_10_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 10 = n59ExactJointKey 10 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_10_matches_exact_prefix_probe]


theorem n59_exact_joint_key_11_value :
    n59ExactJointKey 11 = (71/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_11_value, n59_exact_mass_twice_11_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_11_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 11 = n59ExactJointKey 11 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_11_matches_exact_prefix_probe]


theorem n59_exact_joint_key_12_value :
    n59ExactJointKey 12 = (2 : ℝ) * Real.sqrt (59 : ℝ) + (111/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_12_value, n59_exact_mass_twice_12_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_12_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 12 = n59ExactJointKey 12 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_12_matches_exact_prefix_probe]


theorem n59_exact_joint_key_13_value :
    n59ExactJointKey 13 = (2 : ℝ) * Real.sqrt (59 : ℝ) + (7/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_13_value, n59_exact_mass_twice_13_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_13_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 13 = n59ExactJointKey 13 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_13_matches_exact_prefix_probe]


theorem n59_exact_joint_key_14_value :
    n59ExactJointKey 14 = (1 : ℝ) * Real.sqrt (59 : ℝ) + (91/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_14_value, n59_exact_mass_twice_14_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_14_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 14 = n59ExactJointKey 14 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_14_matches_exact_prefix_probe]


theorem n59_exact_joint_key_15_value :
    n59ExactJointKey 15 = (1 : ℝ) * Real.sqrt (59 : ℝ) + (13/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_15_value, n59_exact_mass_twice_15_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_15_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 15 = n59ExactJointKey 15 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_15_matches_exact_prefix_probe]


theorem n59_exact_joint_key_16_value :
    n59ExactJointKey 16 = (71/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_16_value, n59_exact_mass_twice_16_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_16_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 16 = n59ExactJointKey 16 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_16_matches_exact_prefix_probe]


theorem n59_exact_joint_key_17_value :
    n59ExactJointKey 17 = (33/2 : ℝ) := by
  unfold n59ExactJointKey
  rw [n59_exact_prefix_probe_17_value, n59_exact_mass_twice_17_value]
  have hsq : (Real.sqrt (59 : ℝ)) ^ 2 = (59 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 59)]
  nlinarith


theorem n59_scalar_full_prefix_joint_key_17_matches_exact_joint_key :
    n59ScalarFullPrefixJointKey 17 = n59ExactJointKey 17 := by
  unfold n59ScalarFullPrefixJointKey n59ExactJointKey
  rw [n59_full_prefix_residual_max_17_matches_exact_prefix_probe]


theorem n59_exact_joint_key_7_lt_0 :
    n59ExactJointKey 7 < n59ExactJointKey 0 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_0_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_0 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 0 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_0_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_0


theorem n59_exact_joint_key_7_lt_1 :
    n59ExactJointKey 7 < n59ExactJointKey 1 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_1_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_1 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 1 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_1_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_1


theorem n59_exact_joint_key_7_lt_2 :
    n59ExactJointKey 7 < n59ExactJointKey 2 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_2_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_2 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 2 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_2_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_2


theorem n59_exact_joint_key_7_lt_3 :
    n59ExactJointKey 7 < n59ExactJointKey 3 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_3_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_3 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 3 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_3_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_3


theorem n59_exact_joint_key_7_lt_4 :
    n59ExactJointKey 7 < n59ExactJointKey 4 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_4_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_4 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 4 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_4_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_4


theorem n59_exact_joint_key_7_lt_5 :
    n59ExactJointKey 7 < n59ExactJointKey 5 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_5_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_5 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 5 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_5_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_5


theorem n59_exact_joint_key_7_lt_6 :
    n59ExactJointKey 7 < n59ExactJointKey 6 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_6_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_6 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 6 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_6_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_6


theorem n59_exact_joint_key_7_lt_8 :
    n59ExactJointKey 7 < n59ExactJointKey 8 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_8_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_8 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 8 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_8_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_8


theorem n59_exact_joint_key_7_lt_9 :
    n59ExactJointKey 7 < n59ExactJointKey 9 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_9_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_9 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 9 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_9_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_9


theorem n59_exact_joint_key_7_lt_10 :
    n59ExactJointKey 7 < n59ExactJointKey 10 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_10_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_10 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 10 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_10_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_10


theorem n59_exact_joint_key_7_lt_11 :
    n59ExactJointKey 7 < n59ExactJointKey 11 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_11_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_11 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 11 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_11_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_11


theorem n59_exact_joint_key_7_lt_12 :
    n59ExactJointKey 7 < n59ExactJointKey 12 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_12_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_12 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 12 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_12_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_12


theorem n59_exact_joint_key_7_lt_13 :
    n59ExactJointKey 7 < n59ExactJointKey 13 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_13_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_13 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 13 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_13_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_13


theorem n59_exact_joint_key_7_lt_14 :
    n59ExactJointKey 7 < n59ExactJointKey 14 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_14_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_14 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 14 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_14_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_14


theorem n59_exact_joint_key_7_lt_15 :
    n59ExactJointKey 7 < n59ExactJointKey 15 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_15_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_15 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 15 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_15_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_15


theorem n59_exact_joint_key_7_lt_16 :
    n59ExactJointKey 7 < n59ExactJointKey 16 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_16_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_16 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 16 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_16_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_16


theorem n59_exact_joint_key_7_lt_17 :
    n59ExactJointKey 7 < n59ExactJointKey 17 := by
  rw [n59_exact_joint_key_7_value, n59_exact_joint_key_17_value]
  nlinarith [sqrt59_ge_prefix_lower, sqrt59_lt_prefix_upper]


theorem n59_scalar_full_prefix_joint_key_7_lt_17 :
    n59ScalarFullPrefixJointKey 7 <
      n59ScalarFullPrefixJointKey 17 := by
  rw [n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n59_scalar_full_prefix_joint_key_17_matches_exact_joint_key]
  exact n59_exact_joint_key_7_lt_17


theorem n59_exact_joint_winner_strict_certificate :
    n59ExactJointKey 7 < n59ExactJointKey 0 ∧
    n59ExactJointKey 7 < n59ExactJointKey 1 ∧
    n59ExactJointKey 7 < n59ExactJointKey 2 ∧
    n59ExactJointKey 7 < n59ExactJointKey 3 ∧
    n59ExactJointKey 7 < n59ExactJointKey 4 ∧
    n59ExactJointKey 7 < n59ExactJointKey 5 ∧
    n59ExactJointKey 7 < n59ExactJointKey 6 ∧
    n59ExactJointKey 7 < n59ExactJointKey 8 ∧
    n59ExactJointKey 7 < n59ExactJointKey 9 ∧
    n59ExactJointKey 7 < n59ExactJointKey 10 ∧
    n59ExactJointKey 7 < n59ExactJointKey 11 ∧
    n59ExactJointKey 7 < n59ExactJointKey 12 ∧
    n59ExactJointKey 7 < n59ExactJointKey 13 ∧
    n59ExactJointKey 7 < n59ExactJointKey 14 ∧
    n59ExactJointKey 7 < n59ExactJointKey 15 ∧
    n59ExactJointKey 7 < n59ExactJointKey 16 ∧
    n59ExactJointKey 7 < n59ExactJointKey 17 := by
  exact ⟨n59_exact_joint_key_7_lt_0, ⟨n59_exact_joint_key_7_lt_1, ⟨n59_exact_joint_key_7_lt_2, ⟨n59_exact_joint_key_7_lt_3, ⟨n59_exact_joint_key_7_lt_4, ⟨n59_exact_joint_key_7_lt_5, ⟨n59_exact_joint_key_7_lt_6, ⟨n59_exact_joint_key_7_lt_8, ⟨n59_exact_joint_key_7_lt_9, ⟨n59_exact_joint_key_7_lt_10, ⟨n59_exact_joint_key_7_lt_11, ⟨n59_exact_joint_key_7_lt_12, ⟨n59_exact_joint_key_7_lt_13, ⟨n59_exact_joint_key_7_lt_14, ⟨n59_exact_joint_key_7_lt_15, ⟨n59_exact_joint_key_7_lt_16, n59_exact_joint_key_7_lt_17⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩


theorem n59_scalar_full_prefix_joint_winner_strict_certificate :
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 0 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 1 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 2 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 3 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 4 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 5 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 6 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 8 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 9 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 10 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 11 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 12 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 13 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 14 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 15 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 16 ∧
    n59ScalarFullPrefixJointKey 7 < n59ScalarFullPrefixJointKey 17 := by
  exact ⟨n59_scalar_full_prefix_joint_key_7_lt_0, ⟨n59_scalar_full_prefix_joint_key_7_lt_1, ⟨n59_scalar_full_prefix_joint_key_7_lt_2, ⟨n59_scalar_full_prefix_joint_key_7_lt_3, ⟨n59_scalar_full_prefix_joint_key_7_lt_4, ⟨n59_scalar_full_prefix_joint_key_7_lt_5, ⟨n59_scalar_full_prefix_joint_key_7_lt_6, ⟨n59_scalar_full_prefix_joint_key_7_lt_8, ⟨n59_scalar_full_prefix_joint_key_7_lt_9, ⟨n59_scalar_full_prefix_joint_key_7_lt_10, ⟨n59_scalar_full_prefix_joint_key_7_lt_11, ⟨n59_scalar_full_prefix_joint_key_7_lt_12, ⟨n59_scalar_full_prefix_joint_key_7_lt_13, ⟨n59_scalar_full_prefix_joint_key_7_lt_14, ⟨n59_scalar_full_prefix_joint_key_7_lt_15, ⟨n59_scalar_full_prefix_joint_key_7_lt_16, n59_scalar_full_prefix_joint_key_7_lt_17⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩


theorem n59_prefix_min_5 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59PrefixRank 5 := by
  constructor
  · native_decide
  · intro y hy
    simp [n59FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n59_mass_min_13 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59MassRank 13 := by
  constructor
  · native_decide
  · intro y hy
    simp [n59FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n59_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59JointRank 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n59FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n59_exact_mass_min_13 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ExactMassTwice 13 := by
  constructor
  · native_decide
  · intro y hy
    simp [n59FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n59_exact_prefix_probe_min_5 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ExactPrefixProbe 5 := by
  constructor
  · native_decide
  · intro y hy
    rw [n59_exact_prefix_probe_5_zero]
    unfold n59ExactPrefixProbe Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
    exact le_max_left 0 _


theorem n59_full_prefix_residual_max_min_5 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59FullPrefixResidualMax 5 := by
  constructor
  · native_decide
  · intro y hy
    rw [n59_full_prefix_residual_max_5_zero]
    unfold n59FullPrefixResidualMax
    exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_nonneg 59 (n59WitnessOfIndex y)


theorem n59_exact_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ExactJointKey 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n59FaceIndices] at hy
    interval_cases y
    · exact le_of_lt n59_exact_joint_key_7_lt_0
    · exact le_of_lt n59_exact_joint_key_7_lt_1
    · exact le_of_lt n59_exact_joint_key_7_lt_2
    · exact le_of_lt n59_exact_joint_key_7_lt_3
    · exact le_of_lt n59_exact_joint_key_7_lt_4
    · exact le_of_lt n59_exact_joint_key_7_lt_5
    · exact le_of_lt n59_exact_joint_key_7_lt_6
    · exact le_rfl
    · exact le_of_lt n59_exact_joint_key_7_lt_8
    · exact le_of_lt n59_exact_joint_key_7_lt_9
    · exact le_of_lt n59_exact_joint_key_7_lt_10
    · exact le_of_lt n59_exact_joint_key_7_lt_11
    · exact le_of_lt n59_exact_joint_key_7_lt_12
    · exact le_of_lt n59_exact_joint_key_7_lt_13
    · exact le_of_lt n59_exact_joint_key_7_lt_14
    · exact le_of_lt n59_exact_joint_key_7_lt_15
    · exact le_of_lt n59_exact_joint_key_7_lt_16
    · exact le_of_lt n59_exact_joint_key_7_lt_17


theorem n59_scalar_full_prefix_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ScalarFullPrefixJointKey 7 := by
  exact Erdos.Collider.isFieldMinOn_weightedJoint_of_prefix_eq_on
    (F := n59FaceIndices)
    (x := 7)
    (prefixScalar := n59FullPrefixResidualMax)
    (prefixProbe := n59ExactPrefixProbe)
    (mass := fun i => (n59ExactMassTwice i : ℝ) / 2)
    (weight := Real.sqrt (59 : ℝ))
    n59_exact_joint_min_7
    n59_full_prefix_residual_max_matches_exact_prefix_probe_on_face


theorem n59_prefix_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59PrefixRank n59MassRank := by
  refine ⟨5, 13, n59_prefix_min_5, n59_mass_min_13, ?_⟩
  native_decide

theorem n59_prefix_mass_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_prefix_mass_field_split_full_exported_face


theorem n59_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59PrefixRank n59ExactMassTwice := by
  refine ⟨5, 13, n59_prefix_min_5, n59_exact_mass_min_13, ?_⟩
  native_decide

theorem n59_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_prefix_exact_mass_field_split_full_exported_face


theorem n59_exact_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59ExactPrefixProbe n59ExactMassTwice := by
  refine ⟨5, 13, n59_exact_prefix_probe_min_5, n59_exact_mass_min_13, ?_⟩
  native_decide

theorem n59_exact_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_exact_prefix_exact_mass_field_split_full_exported_face


theorem n59_full_prefix_residual_max_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59FullPrefixResidualMax n59ExactMassTwice := by
  refine ⟨5, 13, n59_full_prefix_residual_max_min_5, n59_exact_mass_min_13, ?_⟩
  native_decide

theorem n59_full_prefix_residual_max_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_full_prefix_residual_max_exact_mass_field_split_full_exported_face


theorem n59_mass_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59MassRank n59JointRank := by
  refine ⟨13, 7, n59_mass_min_13, n59_joint_min_7, ?_⟩
  native_decide

theorem n59_mass_joint_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_mass_joint_field_split_full_exported_face


theorem n59_exact_mass_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59ExactMassTwice n59JointRank := by
  refine ⟨13, 7, n59_exact_mass_min_13, n59_joint_min_7, ?_⟩
  native_decide

theorem n59_exact_mass_joint_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_exact_mass_joint_field_split_full_exported_face


theorem n59_prefix_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n59FaceIndices
      n59PrefixRank n59JointRank := by
  refine ⟨5, 7, n59_prefix_min_5, n59_joint_min_7, ?_⟩
  native_decide

theorem n59_prefix_joint_split_gives_two_exposed_indices :
    2 ≤ n59FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n59_prefix_joint_field_split_full_exported_face



theorem n59_micro_winner_table_matches_packet :
    n59PrefixWinnerIndicesByRank = [5, 7, 10, 11, 16, 17] ∧
    n59MassWinnerIndicesByRank = [13] ∧
    n59JointWinnerIndicesByRank = [7] ∧
    n59ParetoIndicesByRanks = [7, 13] := by
  exact ⟨n59_prefix_winner_match_packet, ⟨n59_mass_winner_match_packet, ⟨n59_joint_winner_match_packet, n59_pareto_minimal_match_packet⟩⟩⟩

theorem n59_micro_exact_mass_winner_table_matches_packet :
    n59MassWinnerIndicesByExactMass = [13] := by
  exact n59_exact_mass_winner_match_packet

theorem n59_micro_split_pattern_certificate :
    Erdos.Collider.FieldSplit n59FaceIndices n59PrefixRank n59MassRank ∧
    Erdos.Collider.FieldSplit n59FaceIndices n59MassRank n59JointRank ∧
    Erdos.Collider.FieldSplit n59FaceIndices n59PrefixRank n59JointRank := by
  exact ⟨n59_prefix_mass_field_split_full_exported_face, ⟨n59_mass_joint_field_split_full_exported_face, n59_prefix_joint_field_split_full_exported_face⟩⟩

theorem n59_micro_exact_prefix_probe_winner_table_matches_packet :
    n59ExactPrefixProbe 5 = 0 ∧
    n59ExactPrefixProbe 7 = 0 ∧
    n59ExactPrefixProbe 10 = 0 ∧
    n59ExactPrefixProbe 11 = 0 ∧
    n59ExactPrefixProbe 16 = 0 ∧
    n59ExactPrefixProbe 17 = 0 := by
  exact ⟨n59_exact_prefix_probe_5_zero, ⟨n59_exact_prefix_probe_7_zero, ⟨n59_exact_prefix_probe_10_zero, ⟨n59_exact_prefix_probe_11_zero, ⟨n59_exact_prefix_probe_16_zero, n59_exact_prefix_probe_17_zero⟩⟩⟩⟩⟩

theorem n59_micro_exact_prefix_probe_strict_nonwinner_certificate :
    0 < n59ExactPrefixProbe 0 ∧
    0 < n59ExactPrefixProbe 1 ∧
    0 < n59ExactPrefixProbe 2 ∧
    0 < n59ExactPrefixProbe 3 ∧
    0 < n59ExactPrefixProbe 4 ∧
    0 < n59ExactPrefixProbe 6 ∧
    0 < n59ExactPrefixProbe 8 ∧
    0 < n59ExactPrefixProbe 9 ∧
    0 < n59ExactPrefixProbe 12 ∧
    0 < n59ExactPrefixProbe 13 ∧
    0 < n59ExactPrefixProbe 14 ∧
    0 < n59ExactPrefixProbe 15 := by
  exact n59_exact_prefix_probe_positive_on_nonwinners

theorem n59_micro_full_prefix_zero_certificate :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W5 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W7 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W10 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W11 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W16 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 59 n59W17 := by
  exact ⟨n59W5_full_prefix_residual_zero, ⟨n59W7_full_prefix_residual_zero, ⟨n59W10_full_prefix_residual_zero, ⟨n59W11_full_prefix_residual_zero, ⟨n59W16_full_prefix_residual_zero, n59W17_full_prefix_residual_zero⟩⟩⟩⟩⟩

theorem n59_micro_packet_probe_separation_is_actual_residual_certificate :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W0 (n59PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W1 (n59PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W2 (n59PacketPrefixLocation 2) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W3 (n59PacketPrefixLocation 3) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W4 (n59PacketPrefixLocation 4) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W6 (n59PacketPrefixLocation 6) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W8 (n59PacketPrefixLocation 8) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W9 (n59PacketPrefixLocation 9) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W12 (n59PacketPrefixLocation 12) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W13 (n59PacketPrefixLocation 13) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W14 (n59PacketPrefixLocation 14) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 59 n59W15 (n59PacketPrefixLocation 15) := by
  exact ⟨n59_prefix_residual_at_packet_location_0_positive, ⟨n59_prefix_residual_at_packet_location_1_positive, ⟨n59_prefix_residual_at_packet_location_2_positive, ⟨n59_prefix_residual_at_packet_location_3_positive, ⟨n59_prefix_residual_at_packet_location_4_positive, ⟨n59_prefix_residual_at_packet_location_6_positive, ⟨n59_prefix_residual_at_packet_location_8_positive, ⟨n59_prefix_residual_at_packet_location_9_positive, ⟨n59_prefix_residual_at_packet_location_12_positive, ⟨n59_prefix_residual_at_packet_location_13_positive, ⟨n59_prefix_residual_at_packet_location_14_positive, n59_prefix_residual_at_packet_location_15_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_micro_full_prefix_scalar_zero_certificate :
    n59FullPrefixResidualMax 5 = 0 ∧
    n59FullPrefixResidualMax 7 = 0 ∧
    n59FullPrefixResidualMax 10 = 0 ∧
    n59FullPrefixResidualMax 11 = 0 ∧
    n59FullPrefixResidualMax 16 = 0 ∧
    n59FullPrefixResidualMax 17 = 0 := by
  exact ⟨n59_full_prefix_residual_max_5_zero, ⟨n59_full_prefix_residual_max_7_zero, ⟨n59_full_prefix_residual_max_10_zero, ⟨n59_full_prefix_residual_max_11_zero, ⟨n59_full_prefix_residual_max_16_zero, n59_full_prefix_residual_max_17_zero⟩⟩⟩⟩⟩

theorem n59_micro_full_prefix_scalar_positive_certificate :
    0 < n59FullPrefixResidualMax 0 ∧
    0 < n59FullPrefixResidualMax 1 ∧
    0 < n59FullPrefixResidualMax 2 ∧
    0 < n59FullPrefixResidualMax 3 ∧
    0 < n59FullPrefixResidualMax 4 ∧
    0 < n59FullPrefixResidualMax 6 ∧
    0 < n59FullPrefixResidualMax 8 ∧
    0 < n59FullPrefixResidualMax 9 ∧
    0 < n59FullPrefixResidualMax 12 ∧
    0 < n59FullPrefixResidualMax 13 ∧
    0 < n59FullPrefixResidualMax 14 ∧
    0 < n59FullPrefixResidualMax 15 := by
  exact ⟨n59_full_prefix_residual_max_0_positive, ⟨n59_full_prefix_residual_max_1_positive, ⟨n59_full_prefix_residual_max_2_positive, ⟨n59_full_prefix_residual_max_3_positive, ⟨n59_full_prefix_residual_max_4_positive, ⟨n59_full_prefix_residual_max_6_positive, ⟨n59_full_prefix_residual_max_8_positive, ⟨n59_full_prefix_residual_max_9_positive, ⟨n59_full_prefix_residual_max_12_positive, ⟨n59_full_prefix_residual_max_13_positive, ⟨n59_full_prefix_residual_max_14_positive, n59_full_prefix_residual_max_15_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_micro_full_prefix_scalar_matches_packet_probe_certificate :
    n59FullPrefixResidualMax 0 = n59ExactPrefixProbe 0 ∧
    n59FullPrefixResidualMax 1 = n59ExactPrefixProbe 1 ∧
    n59FullPrefixResidualMax 2 = n59ExactPrefixProbe 2 ∧
    n59FullPrefixResidualMax 3 = n59ExactPrefixProbe 3 ∧
    n59FullPrefixResidualMax 4 = n59ExactPrefixProbe 4 ∧
    n59FullPrefixResidualMax 5 = n59ExactPrefixProbe 5 ∧
    n59FullPrefixResidualMax 6 = n59ExactPrefixProbe 6 ∧
    n59FullPrefixResidualMax 7 = n59ExactPrefixProbe 7 ∧
    n59FullPrefixResidualMax 8 = n59ExactPrefixProbe 8 ∧
    n59FullPrefixResidualMax 9 = n59ExactPrefixProbe 9 ∧
    n59FullPrefixResidualMax 10 = n59ExactPrefixProbe 10 ∧
    n59FullPrefixResidualMax 11 = n59ExactPrefixProbe 11 ∧
    n59FullPrefixResidualMax 12 = n59ExactPrefixProbe 12 ∧
    n59FullPrefixResidualMax 13 = n59ExactPrefixProbe 13 ∧
    n59FullPrefixResidualMax 14 = n59ExactPrefixProbe 14 ∧
    n59FullPrefixResidualMax 15 = n59ExactPrefixProbe 15 ∧
    n59FullPrefixResidualMax 16 = n59ExactPrefixProbe 16 ∧
    n59FullPrefixResidualMax 17 = n59ExactPrefixProbe 17 := by
  exact ⟨n59_full_prefix_residual_max_0_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_1_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_2_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_3_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_4_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_5_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_6_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_7_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_8_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_9_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_10_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_11_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_12_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_13_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_14_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_15_matches_exact_prefix_probe, ⟨n59_full_prefix_residual_max_16_matches_exact_prefix_probe, n59_full_prefix_residual_max_17_matches_exact_prefix_probe⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_micro_scalar_joint_matches_probe_joint_certificate :
    n59ScalarFullPrefixJointKey 0 = n59ExactJointKey 0 ∧
    n59ScalarFullPrefixJointKey 1 = n59ExactJointKey 1 ∧
    n59ScalarFullPrefixJointKey 2 = n59ExactJointKey 2 ∧
    n59ScalarFullPrefixJointKey 3 = n59ExactJointKey 3 ∧
    n59ScalarFullPrefixJointKey 4 = n59ExactJointKey 4 ∧
    n59ScalarFullPrefixJointKey 5 = n59ExactJointKey 5 ∧
    n59ScalarFullPrefixJointKey 6 = n59ExactJointKey 6 ∧
    n59ScalarFullPrefixJointKey 7 = n59ExactJointKey 7 ∧
    n59ScalarFullPrefixJointKey 8 = n59ExactJointKey 8 ∧
    n59ScalarFullPrefixJointKey 9 = n59ExactJointKey 9 ∧
    n59ScalarFullPrefixJointKey 10 = n59ExactJointKey 10 ∧
    n59ScalarFullPrefixJointKey 11 = n59ExactJointKey 11 ∧
    n59ScalarFullPrefixJointKey 12 = n59ExactJointKey 12 ∧
    n59ScalarFullPrefixJointKey 13 = n59ExactJointKey 13 ∧
    n59ScalarFullPrefixJointKey 14 = n59ExactJointKey 14 ∧
    n59ScalarFullPrefixJointKey 15 = n59ExactJointKey 15 ∧
    n59ScalarFullPrefixJointKey 16 = n59ExactJointKey 16 ∧
    n59ScalarFullPrefixJointKey 17 = n59ExactJointKey 17 := by
  exact ⟨n59_scalar_full_prefix_joint_key_0_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_1_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_2_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_3_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_4_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_5_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_6_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_7_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_8_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_9_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_10_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_11_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_12_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_13_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_14_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_15_matches_exact_joint_key, ⟨n59_scalar_full_prefix_joint_key_16_matches_exact_joint_key, n59_scalar_full_prefix_joint_key_17_matches_exact_joint_key⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n59_micro_exact_joint_winner_certificate :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ExactJointKey 7 := by
  exact n59_exact_joint_min_7

theorem n59_micro_scalar_full_prefix_joint_winner_certificate :
    Erdos.Collider.IsFieldMinOn n59FaceIndices n59ScalarFullPrefixJointKey 7 := by
  exact n59_scalar_full_prefix_joint_min_7

end Erdos30FaceField59FullFaceCertificate
