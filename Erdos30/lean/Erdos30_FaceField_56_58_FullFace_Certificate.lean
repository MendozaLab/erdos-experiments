import Erdos30_FaceField
import Erdos30_FaceField_ExactObservables
import Erdos30_Sidon_Defs
import Mathlib

/-!
# Erdos #30 n = 56..58 Full Exported-Face Certificate

This file instantiates the face/field language over the complete exported
ground faces in the Maxwell handoff packet:

`EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.json`.

For each row `n = 56,57,58`, every exported face member is named, checked as a
Sidon set of size `h(n) = 10`, checked to lie in `[0,n]`, tied to a
rank-coded version of the packet observables, and checked against the exact
integer density-adjusted mass observable. The rank codes preserve the packet's
ordering of prefix and joint scores on the exported face; they are not raw
floating-point theorem statements. The mass winner theorem no longer depends
on packet floats.

This is a finite exported-face certificate only. It is not a proof of Erdos
#30, not a proof that the export was exhaustive beyond the source packet's
exact-count claim, and not an asymptotic theorem.
-/

open Finset Nat

namespace Erdos30FaceField5658FullFaceCertificate

/-! ## n = 56: complete exported face -/

def n56W0 : Finset Nat :=
  ([0, 1, 6, 10, 23, 26, 34, 41, 53, 55] : List Nat).toFinset

def n56W1 : Finset Nat :=
  ([0, 2, 14, 21, 29, 32, 45, 49, 54, 55] : List Nat).toFinset

def n56W2 : Finset Nat :=
  ([1, 2, 7, 11, 24, 27, 35, 42, 54, 56] : List Nat).toFinset

def n56W3 : Finset Nat :=
  ([1, 3, 15, 22, 30, 33, 46, 50, 55, 56] : List Nat).toFinset

def n56Face : Finset (Finset Nat) :=
  ([n56W0, n56W1, n56W2, n56W3] : List (Finset Nat)).toFinset

def n56FaceIndices : Finset Nat :=
  Finset.range 4

def n56IndexList : List Nat :=
  List.range 4

def n56WitnessOfIndex : Nat -> Finset Nat

  | 0 => n56W0

  | 1 => n56W1

  | 2 => n56W2

  | 3 => n56W3

  | _ => ∅

theorem n56_indexed_face_matches_exported_face :
    (n56IndexList.map n56WitnessOfIndex).toFinset = n56Face := by
  native_decide

theorem n56_exported_face_card :
    n56Face.card = 4 := by
  native_decide

theorem n56_face_indices_card :
    n56FaceIndices.card = 4 := by
  native_decide

theorem n56_all_exported_witnesses_have_card_h :
    n56W0.card = 10 ∧
    n56W1.card = 10 ∧
    n56W2.card = 10 ∧
    n56W3.card = 10 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩

theorem n56_all_exported_witnesses_in_range :
    n56W0 ⊆ Finset.range 57 ∧
    n56W1 ⊆ Finset.range 57 ∧
    n56W2 ⊆ Finset.range 57 ∧
    n56W3 ⊆ Finset.range 57 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩

theorem n56_all_exported_witnesses_are_sidon :
    Erdos.Sidon.IsSidonSet n56W0 ∧
    Erdos.Sidon.IsSidonSet n56W1 ∧
    Erdos.Sidon.IsSidonSet n56W2 ∧
    Erdos.Sidon.IsSidonSet n56W3 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩

def n56PrefixRank : Nat -> Nat
  | 0 => 3
  | 1 => 2
  | 2 => 1
  | 3 => 0
  | _ => 99

def n56MassRank : Nat -> Nat
  | 0 => 3
  | 1 => 1
  | 2 => 2
  | 3 => 0
  | _ => 99

def n56JointRank : Nat -> Nat
  | 0 => 3
  | 1 => 1
  | 2 => 2
  | 3 => 0
  | _ => 99

def n56ExactMassTwice (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.densityAdjustedMassTwice 56 (n56WitnessOfIndex i)

def n56MassWinnerIndicesByExactMass : List Nat :=
  n56IndexList.filter (fun i => n56ExactMassTwice i == 6)

def n56PacketPrefixLocation : Nat -> Nat
  | 0 => 10
  | 1 => 55
  | 2 => 11
  | 3 => 56
  | _ => 0

def n56PacketPrefixCount : Nat -> Nat
  | 0 => 4
  | 1 => 10
  | 2 => 4
  | 3 => 10
  | _ => 0

noncomputable def n56ExactPrefixProbe (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.prefixResidualProbeCard10 56 (n56PacketPrefixLocation i) (n56PacketPrefixCount i)

noncomputable def n56FullPrefixResidualMax (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.fullPrefixResidualMax 56 (n56WitnessOfIndex i)


noncomputable def n56ExactJointKey (i : Nat) : ℝ :=
  n56ExactPrefixProbe i * Real.sqrt (56 : ℝ) +
    (n56ExactMassTwice i : ℝ) / 2


noncomputable def n56ScalarFullPrefixJointKey (i : Nat) : ℝ :=
  n56FullPrefixResidualMax i * Real.sqrt (56 : ℝ) +
    (n56ExactMassTwice i : ℝ) / 2


def n56PrefixCountAtPacketLocation (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.prefixCount (n56WitnessOfIndex i) (n56PacketPrefixLocation i)

def n56PrefixWinnerIndicesByRank : List Nat :=
  n56IndexList.filter (fun i => n56PrefixRank i == 0)

def n56MassWinnerIndicesByRank : List Nat :=
  n56IndexList.filter (fun i => n56MassRank i == 0)

def n56JointWinnerIndicesByRank : List Nat :=
  n56IndexList.filter (fun i => n56JointRank i == 0)

def n56DominatesPrefixMass (i j : Nat) : Bool :=
  (n56PrefixRank i <= n56PrefixRank j) &&
  (n56MassRank i <= n56MassRank j) &&
  ((n56PrefixRank i < n56PrefixRank j) ||
   (n56MassRank i < n56MassRank j))

def n56ParetoIndicesByRanks : List Nat :=
  n56IndexList.filter (fun i => !(n56IndexList.any (fun j => n56DominatesPrefixMass j i)))

theorem n56_prefix_winner_match_packet :
    n56PrefixWinnerIndicesByRank = [3] := by
  native_decide

theorem n56_mass_winner_match_packet :
    n56MassWinnerIndicesByRank = [3] := by
  native_decide

theorem n56_joint_winner_match_packet :
    n56JointWinnerIndicesByRank = [3] := by
  native_decide

theorem n56_pareto_minimal_match_packet :
    n56ParetoIndicesByRanks = [3] := by
  native_decide

theorem n56_exact_mass_twice_table :
    n56IndexList.map n56ExactMassTwice = [118, 14, 98, 6] := by
  native_decide

theorem n56_exact_mass_winner_match_packet :
    n56MassWinnerIndicesByExactMass = [3] := by
  native_decide

theorem n56_exact_mass_winner_matches_rank_winner :
    n56MassWinnerIndicesByExactMass = n56MassWinnerIndicesByRank := by
  native_decide

theorem n56_packet_prefix_location_table :
    n56IndexList.map n56PacketPrefixLocation = [10, 55, 11, 56] := by
  native_decide

theorem n56_prefix_count_at_packet_location_table :
    n56IndexList.map n56PrefixCountAtPacketLocation = [4, 10, 4, 10] := by
  native_decide

lemma sqrt56_ge_prefix_lower : (7 : ℝ) ≤ Real.sqrt (56 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 7) (by norm_num : (0:ℝ) ≤ 56)]
  norm_num

lemma sqrt56_lt_prefix_upper : Real.sqrt (56 : ℝ) < (8 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 8)]
  norm_num

lemma sqrt56_lt_prefix_tight_upper : Real.sqrt (56 : ℝ) < (15/2 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 15/2)]
  norm_num

def n56W3PrefixCountTable : Nat -> Nat
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
  | 15 => 3
  | 16 => 3
  | 17 => 3
  | 18 => 3
  | 19 => 3
  | 20 => 3
  | 21 => 3
  | 22 => 4
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 4
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
  | 46 => 7
  | 47 => 7
  | 48 => 7
  | 49 => 7
  | 50 => 8
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 9
  | 56 => 10
  | _ => 10
theorem n56W3_prefix_count_table_matches_witness :
    (List.range 57).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n56W3 t) =
      [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

theorem n56W3_prefix_count_model_table :
    (List.range 57).map n56W3PrefixCountTable =
      [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

lemma n56W3_prefix_segment_0_0_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_1_2_bound :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_1_2 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_3_14_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_3_14 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_15_21_bound :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_15_21 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_22_29_bound :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_22_29 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_30_32_bound :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_30_32 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_33_45_bound :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_33_45 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_46_49_bound :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_50_54_bound :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_55_55_bound :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_55_55 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n56W3_prefix_segment_56_56_bound :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper, hlowR, hhighR]

theorem n56W3_prefix_count_segment_56_56 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W3 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W3_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) ∧
    (∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        10 * Real.sqrt (56 : ℝ) - 56) := by
  exact ⟨n56W3_prefix_segment_0_0_bound, ⟨n56W3_prefix_segment_1_2_bound, ⟨n56W3_prefix_segment_3_14_bound, ⟨n56W3_prefix_segment_15_21_bound, ⟨n56W3_prefix_segment_22_29_bound, ⟨n56W3_prefix_segment_30_32_bound, ⟨n56W3_prefix_segment_33_45_bound, ⟨n56W3_prefix_segment_46_49_bound, ⟨n56W3_prefix_segment_50_54_bound, ⟨n56W3_prefix_segment_55_55_bound, n56W3_prefix_segment_56_56_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n56W3_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 56 n56W3 =
      10 * Real.sqrt (56 : ℝ) - 56 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n56W3.card = 10 := by
    native_decide
  have hcard : (n56W3.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 56)]
  nlinarith [hsq]

theorem n56W3_prefix_residual_segment_0_0_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_0_0_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_1_2_zero :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_1_2_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_1_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_3_14_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_3_14_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_3_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_15_21_zero :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_15_21_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_15_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_22_29_zero :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_22_29_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_22_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_30_32_zero :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_30_32_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_30_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_33_45_zero :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_33_45_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_33_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_46_49_zero :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_46_49_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_50_54_zero :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_50_54_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_55_55_zero :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_55_55_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_55_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_prefix_residual_segment_56_56_zero :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 t = 0 := by
  intro t hlow hhigh
  have hdev := n56W3_prefix_segment_56_56_bound t hlow hhigh
  have hcount := n56W3_prefix_count_segment_56_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n56W3_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 56 n56W3 := by
  intro t ht
  interval_cases t
  · exact n56W3_prefix_residual_segment_0_0_zero 0 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_1_2_zero 1 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_1_2_zero 2 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 3 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 4 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 5 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 6 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 7 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 8 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 9 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 10 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 11 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 12 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 13 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_3_14_zero 14 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 15 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 16 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 17 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 18 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 19 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 20 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_15_21_zero 21 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 22 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 23 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 24 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 25 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 26 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 27 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 28 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_22_29_zero 29 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_30_32_zero 30 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_30_32_zero 31 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_30_32_zero 32 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 33 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 34 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 35 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 36 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 37 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 38 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 39 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 40 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 41 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 42 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 43 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 44 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_33_45_zero 45 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_46_49_zero 46 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_46_49_zero 47 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_46_49_zero 48 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_46_49_zero 49 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_50_54_zero 50 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_50_54_zero 51 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_50_54_zero 52 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_50_54_zero 53 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_50_54_zero 54 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_55_55_zero 55 (by norm_num) (by norm_num)
  · exact n56W3_prefix_residual_segment_56_56_zero 56 (by norm_num) (by norm_num)


lemma sqrt56_ge_56_div_10 : (56/10 : ℝ) ≤ Real.sqrt (56 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 56/10) (by norm_num : (0:ℝ) ≤ 56)]
  norm_num


theorem n56W0_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 56 n56W0 =
      10 * Real.sqrt (56 : ℝ) - 56 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n56W0.card = 10 := by
    native_decide
  have hcard : (n56W0.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 56)]
  nlinarith [hsq]


theorem n56W1_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 56 n56W1 =
      10 * Real.sqrt (56 : ℝ) - 56 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n56W1.card = 10 := by
    native_decide
  have hcard : (n56W1.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 56)]
  nlinarith [hsq]


theorem n56W2_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 56 n56W2 =
      10 * Real.sqrt (56 : ℝ) - 56 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n56W2.card = 10 := by
    native_decide
  have hcard : (n56W2.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (56 : ℝ) := by
    nlinarith [sqrt56_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 56)]
  nlinarith [hsq]


theorem n56_exact_prefix_probe_3_zero :
    n56ExactPrefixProbe 3 = 0 := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_56_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n56_exact_prefix_probe_0_positive :
    0 < n56ExactPrefixProbe 0 := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (10 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (56 : ℝ) < (23/3 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
    norm_num
  have hgap :
      0 < -((10 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)) -
        (10 * Real.sqrt (56 : ℝ) - 56) := by
    nlinarith [hsqrtUpper]
  linarith [hgap]


theorem n56_exact_prefix_probe_1_positive :
    0 < n56ExactPrefixProbe 1 := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_56_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n56_exact_prefix_probe_2_positive :
    0 < n56ExactPrefixProbe 2 := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (11 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (56 : ℝ) < (15/2 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 15/2)]
    norm_num
  have hgap :
      0 < -((11 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)) -
        (10 * Real.sqrt (56 : ℝ) - 56) := by
    nlinarith [hsqrtUpper]
  linarith [hgap]


theorem n56_exact_prefix_probe_positive_on_nonwinners :
    0 < n56ExactPrefixProbe 0 ∧
    0 < n56ExactPrefixProbe 1 ∧
    0 < n56ExactPrefixProbe 2 := by
  exact ⟨n56_exact_prefix_probe_0_positive, ⟨n56_exact_prefix_probe_1_positive, n56_exact_prefix_probe_2_positive⟩⟩


theorem n56_exact_prefix_probe_0_matches_prefix_residual_at_packet_location :
    n56ExactPrefixProbe 0 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 (n56PacketPrefixLocation 0) := by
  unfold n56ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n56W0 (n56PacketPrefixLocation 0) =
        n56PacketPrefixCount 0 := by
    native_decide
  rw [hcount]
  rfl


theorem n56_exact_prefix_probe_1_matches_prefix_residual_at_packet_location :
    n56ExactPrefixProbe 1 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 (n56PacketPrefixLocation 1) := by
  unfold n56ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n56W1 (n56PacketPrefixLocation 1) =
        n56PacketPrefixCount 1 := by
    native_decide
  rw [hcount]
  rfl


theorem n56_exact_prefix_probe_2_matches_prefix_residual_at_packet_location :
    n56ExactPrefixProbe 2 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 (n56PacketPrefixLocation 2) := by
  unfold n56ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n56W2 (n56PacketPrefixLocation 2) =
        n56PacketPrefixCount 2 := by
    native_decide
  rw [hcount]
  rfl


theorem n56_exact_prefix_probe_3_matches_prefix_residual_at_packet_location :
    n56ExactPrefixProbe 3 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W3 (n56PacketPrefixLocation 3) := by
  unfold n56ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n56W3 (n56PacketPrefixLocation 3) =
        n56PacketPrefixCount 3 := by
    native_decide
  rw [hcount]
  rfl


theorem n56_prefix_residual_at_packet_location_0_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 (n56PacketPrefixLocation 0) := by
  rw [← n56_exact_prefix_probe_0_matches_prefix_residual_at_packet_location]
  exact n56_exact_prefix_probe_0_positive


theorem n56_prefix_residual_at_packet_location_1_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 (n56PacketPrefixLocation 1) := by
  rw [← n56_exact_prefix_probe_1_matches_prefix_residual_at_packet_location]
  exact n56_exact_prefix_probe_1_positive


theorem n56_prefix_residual_at_packet_location_2_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 (n56PacketPrefixLocation 2) := by
  rw [← n56_exact_prefix_probe_2_matches_prefix_residual_at_packet_location]
  exact n56_exact_prefix_probe_2_positive


theorem n56_full_prefix_residual_zero_on_prefix_winners :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 56 n56W3 := by
  exact n56W3_full_prefix_residual_zero


theorem n56_actual_prefix_residual_positive_on_probe_nonwinners :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 (n56PacketPrefixLocation 0) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 (n56PacketPrefixLocation 1) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 (n56PacketPrefixLocation 2) := by
  exact ⟨n56_prefix_residual_at_packet_location_0_positive, ⟨n56_prefix_residual_at_packet_location_1_positive, n56_prefix_residual_at_packet_location_2_positive⟩⟩


theorem n56_full_prefix_residual_max_3_zero :
    n56FullPrefixResidualMax 3 = 0 := by
  unfold n56FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n56W3_full_prefix_residual_zero


theorem n56_full_prefix_residual_max_0_positive :
    0 < n56FullPrefixResidualMax 0 := by
  unfold n56FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n56PacketPrefixLocation 0 ≤ 56)
    n56_prefix_residual_at_packet_location_0_positive


theorem n56_full_prefix_residual_max_1_positive :
    0 < n56FullPrefixResidualMax 1 := by
  unfold n56FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n56PacketPrefixLocation 1 ≤ 56)
    n56_prefix_residual_at_packet_location_1_positive


theorem n56_full_prefix_residual_max_2_positive :
    0 < n56FullPrefixResidualMax 2 := by
  unfold n56FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n56PacketPrefixLocation 2 ≤ 56)
    n56_prefix_residual_at_packet_location_2_positive


theorem n56_full_prefix_residual_max_zero_on_prefix_winners :
    n56FullPrefixResidualMax 3 = 0 := by
  exact n56_full_prefix_residual_max_3_zero


theorem n56_full_prefix_residual_max_positive_on_probe_nonwinners :
    0 < n56FullPrefixResidualMax 0 ∧
        0 < n56FullPrefixResidualMax 1 ∧
        0 < n56FullPrefixResidualMax 2 := by
  exact ⟨n56_full_prefix_residual_max_0_positive, ⟨n56_full_prefix_residual_max_1_positive, n56_full_prefix_residual_max_2_positive⟩⟩


theorem n56_exact_prefix_probe_0_value :
    n56ExactPrefixProbe 0 = ((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (10 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (56 : ℝ) < (23/3 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
    norm_num
  have hgap :
      0 < -((10 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)) -
        (10 * Real.sqrt (56 : ℝ) - 56) := by
    nlinarith [hsqrtUpper]
  rw [max_eq_right (le_of_lt hgap)]
  ring


theorem n56_exact_mass_twice_0_value :
    n56ExactMassTwice 0 = 118 := by
  native_decide


theorem n56_exact_prefix_probe_1_value :
    n56ExactPrefixProbe 1 = (1 : ℝ) := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_56_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n56_exact_mass_twice_1_value :
    n56ExactMassTwice 1 = 14 := by
  native_decide


theorem n56_exact_prefix_probe_2_value :
    n56ExactPrefixProbe 2 = ((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (11 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (56 : ℝ) < (15/2 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 15/2)]
    norm_num
  have hgap :
      0 < -((11 : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)) -
        (10 * Real.sqrt (56 : ℝ) - 56) := by
    nlinarith [hsqrtUpper]
  rw [max_eq_right (le_of_lt hgap)]
  ring


theorem n56_exact_mass_twice_2_value :
    n56ExactMassTwice 2 = 98 := by
  native_decide


theorem n56_exact_prefix_probe_3_value :
    n56ExactPrefixProbe 3 = (0 : ℝ) := by
  simp [n56ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n56PacketPrefixLocation, n56PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (56 : ℝ) ≤ 0 := by
    nlinarith [sqrt56_ge_56_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n56_exact_mass_twice_3_value :
    n56ExactMassTwice 3 = 6 := by
  native_decide


theorem n56W0_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_1_5 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_1_5_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (5 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_1_5 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_6_9 :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_6_9_le_exact_probe :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (6 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_6_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_10_22 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_10_22_le_exact_probe :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_10_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_23_25 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_23_25_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_23_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_26_33 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_26_33_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_26_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_34_40 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_34_40_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (40 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_34_40 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_41_52 :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_41_52_le_exact_probe :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (41 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (52 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_41_52 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_53_54 :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_53_54_le_exact_probe :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (53 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_53_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_prefix_count_segment_55_56 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W0 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W0_prefix_residual_segment_55_56_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W0_prefix_count_segment_55_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_0_positive
  · rw [n56_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W0_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 t ≤
        n56ExactPrefixProbe 0 := by
  intro t ht
  interval_cases t
  · exact n56W0_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_1_5_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_1_5_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_1_5_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_1_5_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_1_5_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_6_9_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_6_9_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_6_9_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_6_9_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_10_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_23_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_23_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_23_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_26_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_34_40_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_41_52_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_53_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_53_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_55_56_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n56W0_prefix_residual_segment_55_56_le_exact_probe 56 (by norm_num) (by norm_num)

theorem n56W1_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_2_13 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_2_13_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_2_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_14_20 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_14_20_le_exact_probe :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_14_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_21_28 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_21_28_le_exact_probe :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_21_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_29_31 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_29_31_le_exact_probe :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_29_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_32_44 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_32_44_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_32_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_45_48 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_45_48_le_exact_probe :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_45_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_49_53 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_49_53_le_exact_probe :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_49_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_54_54 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_54_54_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_54_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_prefix_count_segment_55_56 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W1 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W1_prefix_residual_segment_55_56_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W1_prefix_count_segment_55_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_1_positive
  · rw [n56_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W1_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 t ≤
        n56ExactPrefixProbe 1 := by
  intro t ht
  interval_cases t
  · exact n56W1_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_2_13_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_14_20_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_21_28_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_29_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_29_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_29_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_32_44_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_45_48_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_45_48_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_45_48_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_45_48_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_49_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_49_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_49_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_49_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_49_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_54_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_55_56_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n56W1_prefix_residual_segment_55_56_le_exact_probe 56 (by norm_num) (by norm_num)

theorem n56W2_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_1_1 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_1_1_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_1_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_2_6 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_2_6_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (6 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_2_6 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_7_10 :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_7_10_le_exact_probe :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (7 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (10 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_7_10 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_11_23 :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_11_23_le_exact_probe :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (11 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_11_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_24_26 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_24_26_le_exact_probe :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_24_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_27_34 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_27_34_le_exact_probe :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_27_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_35_41 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_35_41_le_exact_probe :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (41 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_35_41 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_42_53 :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_42_53_le_exact_probe :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (42 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_42_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_54_55 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_54_55_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_54_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_prefix_count_segment_56_56 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n56W2 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n56W2_prefix_residual_segment_56_56_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (56 : ℝ)| ≤
        (10 * Real.sqrt (56 : ℝ) - 56) + (((45 : ℝ) - (6 : ℝ) * Real.sqrt (56 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper,
        sqrt56_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n56W2_prefix_count_segment_56_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n56W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n56_exact_prefix_probe_2_positive
  · rw [n56_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n56W2_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 t ≤
        n56ExactPrefixProbe 2 := by
  intro t ht
  interval_cases t
  · exact n56W2_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_1_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_2_6_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_2_6_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_2_6_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_2_6_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_2_6_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_7_10_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_7_10_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_7_10_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_7_10_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_11_23_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_24_26_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_24_26_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_24_26_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_27_34_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_35_41_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_42_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_54_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_54_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n56W2_prefix_residual_segment_56_56_le_exact_probe 56 (by norm_num) (by norm_num)


theorem n56_full_prefix_residual_max_0_matches_exact_prefix_probe :
    n56FullPrefixResidualMax 0 = n56ExactPrefixProbe 0 := by
  unfold n56FullPrefixResidualMax n56WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n56PacketPrefixLocation 0 ≤ 56)
    n56W0_full_prefix_residual_le_exact_probe
    (n56_exact_prefix_probe_0_matches_prefix_residual_at_packet_location).symm


theorem n56_full_prefix_residual_max_1_matches_exact_prefix_probe :
    n56FullPrefixResidualMax 1 = n56ExactPrefixProbe 1 := by
  unfold n56FullPrefixResidualMax n56WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n56PacketPrefixLocation 1 ≤ 56)
    n56W1_full_prefix_residual_le_exact_probe
    (n56_exact_prefix_probe_1_matches_prefix_residual_at_packet_location).symm


theorem n56_full_prefix_residual_max_2_matches_exact_prefix_probe :
    n56FullPrefixResidualMax 2 = n56ExactPrefixProbe 2 := by
  unfold n56FullPrefixResidualMax n56WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n56PacketPrefixLocation 2 ≤ 56)
    n56W2_full_prefix_residual_le_exact_probe
    (n56_exact_prefix_probe_2_matches_prefix_residual_at_packet_location).symm


theorem n56_full_prefix_residual_max_3_matches_exact_prefix_probe :
    n56FullPrefixResidualMax 3 = n56ExactPrefixProbe 3 := by
  rw [n56_full_prefix_residual_max_3_zero, n56_exact_prefix_probe_3_zero]


theorem n56_full_prefix_residual_max_matches_exact_prefix_probe_on_face :
    ∀ i ∈ n56FaceIndices, n56FullPrefixResidualMax i = n56ExactPrefixProbe i := by
  intro i hi
  simp [n56FaceIndices] at hi
  interval_cases i
  · exact n56_full_prefix_residual_max_0_matches_exact_prefix_probe
  · exact n56_full_prefix_residual_max_1_matches_exact_prefix_probe
  · exact n56_full_prefix_residual_max_2_matches_exact_prefix_probe
  · exact n56_full_prefix_residual_max_3_matches_exact_prefix_probe


theorem n56_exact_joint_key_0_value :
    n56ExactJointKey 0 = (46 : ℝ) * Real.sqrt (56 : ℝ) - (277 : ℝ) := by
  unfold n56ExactJointKey
  rw [n56_exact_prefix_probe_0_value, n56_exact_mass_twice_0_value]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 56)]
  nlinarith


theorem n56_scalar_full_prefix_joint_key_0_matches_exact_joint_key :
    n56ScalarFullPrefixJointKey 0 = n56ExactJointKey 0 := by
  unfold n56ScalarFullPrefixJointKey n56ExactJointKey
  rw [n56_full_prefix_residual_max_0_matches_exact_prefix_probe]


theorem n56_exact_joint_key_1_value :
    n56ExactJointKey 1 = (1 : ℝ) * Real.sqrt (56 : ℝ) + (7 : ℝ) := by
  unfold n56ExactJointKey
  rw [n56_exact_prefix_probe_1_value, n56_exact_mass_twice_1_value]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 56)]
  nlinarith


theorem n56_scalar_full_prefix_joint_key_1_matches_exact_joint_key :
    n56ScalarFullPrefixJointKey 1 = n56ExactJointKey 1 := by
  unfold n56ScalarFullPrefixJointKey n56ExactJointKey
  rw [n56_full_prefix_residual_max_1_matches_exact_prefix_probe]


theorem n56_exact_joint_key_2_value :
    n56ExactJointKey 2 = (45 : ℝ) * Real.sqrt (56 : ℝ) - (287 : ℝ) := by
  unfold n56ExactJointKey
  rw [n56_exact_prefix_probe_2_value, n56_exact_mass_twice_2_value]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 56)]
  nlinarith


theorem n56_scalar_full_prefix_joint_key_2_matches_exact_joint_key :
    n56ScalarFullPrefixJointKey 2 = n56ExactJointKey 2 := by
  unfold n56ScalarFullPrefixJointKey n56ExactJointKey
  rw [n56_full_prefix_residual_max_2_matches_exact_prefix_probe]


theorem n56_exact_joint_key_3_value :
    n56ExactJointKey 3 = (3 : ℝ) := by
  unfold n56ExactJointKey
  rw [n56_exact_prefix_probe_3_value, n56_exact_mass_twice_3_value]
  have hsq : (Real.sqrt (56 : ℝ)) ^ 2 = (56 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 56)]
  nlinarith


theorem n56_scalar_full_prefix_joint_key_3_matches_exact_joint_key :
    n56ScalarFullPrefixJointKey 3 = n56ExactJointKey 3 := by
  unfold n56ScalarFullPrefixJointKey n56ExactJointKey
  rw [n56_full_prefix_residual_max_3_matches_exact_prefix_probe]


theorem n56_exact_joint_key_3_lt_0 :
    n56ExactJointKey 3 < n56ExactJointKey 0 := by
  rw [n56_exact_joint_key_3_value, n56_exact_joint_key_0_value]
  nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper]


theorem n56_scalar_full_prefix_joint_key_3_lt_0 :
    n56ScalarFullPrefixJointKey 3 <
      n56ScalarFullPrefixJointKey 0 := by
  rw [n56_scalar_full_prefix_joint_key_3_matches_exact_joint_key,
    n56_scalar_full_prefix_joint_key_0_matches_exact_joint_key]
  exact n56_exact_joint_key_3_lt_0


theorem n56_exact_joint_key_3_lt_1 :
    n56ExactJointKey 3 < n56ExactJointKey 1 := by
  rw [n56_exact_joint_key_3_value, n56_exact_joint_key_1_value]
  nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper]


theorem n56_scalar_full_prefix_joint_key_3_lt_1 :
    n56ScalarFullPrefixJointKey 3 <
      n56ScalarFullPrefixJointKey 1 := by
  rw [n56_scalar_full_prefix_joint_key_3_matches_exact_joint_key,
    n56_scalar_full_prefix_joint_key_1_matches_exact_joint_key]
  exact n56_exact_joint_key_3_lt_1


theorem n56_exact_joint_key_3_lt_2 :
    n56ExactJointKey 3 < n56ExactJointKey 2 := by
  rw [n56_exact_joint_key_3_value, n56_exact_joint_key_2_value]
  nlinarith [sqrt56_ge_prefix_lower, sqrt56_lt_prefix_upper]


theorem n56_scalar_full_prefix_joint_key_3_lt_2 :
    n56ScalarFullPrefixJointKey 3 <
      n56ScalarFullPrefixJointKey 2 := by
  rw [n56_scalar_full_prefix_joint_key_3_matches_exact_joint_key,
    n56_scalar_full_prefix_joint_key_2_matches_exact_joint_key]
  exact n56_exact_joint_key_3_lt_2


theorem n56_exact_joint_winner_strict_certificate :
    n56ExactJointKey 3 < n56ExactJointKey 0 ∧
    n56ExactJointKey 3 < n56ExactJointKey 1 ∧
    n56ExactJointKey 3 < n56ExactJointKey 2 := by
  exact ⟨n56_exact_joint_key_3_lt_0, ⟨n56_exact_joint_key_3_lt_1, n56_exact_joint_key_3_lt_2⟩⟩


theorem n56_scalar_full_prefix_joint_winner_strict_certificate :
    n56ScalarFullPrefixJointKey 3 < n56ScalarFullPrefixJointKey 0 ∧
    n56ScalarFullPrefixJointKey 3 < n56ScalarFullPrefixJointKey 1 ∧
    n56ScalarFullPrefixJointKey 3 < n56ScalarFullPrefixJointKey 2 := by
  exact ⟨n56_scalar_full_prefix_joint_key_3_lt_0, ⟨n56_scalar_full_prefix_joint_key_3_lt_1, n56_scalar_full_prefix_joint_key_3_lt_2⟩⟩


theorem n56_prefix_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56PrefixRank 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n56FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n56_mass_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56MassRank 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n56FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n56_joint_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56JointRank 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n56FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n56_exact_mass_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ExactMassTwice 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n56FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n56_exact_prefix_probe_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ExactPrefixProbe 3 := by
  constructor
  · native_decide
  · intro y hy
    rw [n56_exact_prefix_probe_3_zero]
    unfold n56ExactPrefixProbe Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
    exact le_max_left 0 _


theorem n56_full_prefix_residual_max_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56FullPrefixResidualMax 3 := by
  constructor
  · native_decide
  · intro y hy
    rw [n56_full_prefix_residual_max_3_zero]
    unfold n56FullPrefixResidualMax
    exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_nonneg 56 (n56WitnessOfIndex y)


theorem n56_exact_joint_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ExactJointKey 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n56FaceIndices] at hy
    interval_cases y
    · exact le_of_lt n56_exact_joint_key_3_lt_0
    · exact le_of_lt n56_exact_joint_key_3_lt_1
    · exact le_of_lt n56_exact_joint_key_3_lt_2
    · exact le_rfl


theorem n56_scalar_full_prefix_joint_min_3 :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ScalarFullPrefixJointKey 3 := by
  exact Erdos.Collider.isFieldMinOn_weightedJoint_of_prefix_eq_on
    (F := n56FaceIndices)
    (x := 3)
    (prefixScalar := n56FullPrefixResidualMax)
    (prefixProbe := n56ExactPrefixProbe)
    (mass := fun i => (n56ExactMassTwice i : ℝ) / 2)
    (weight := Real.sqrt (56 : ℝ))
    n56_exact_joint_min_3
    n56_full_prefix_residual_max_matches_exact_prefix_probe_on_face

/-! ## n = 57: complete exported face -/

def n57W0 : Finset Nat :=
  ([0, 1, 6, 10, 23, 26, 34, 41, 53, 55] : List Nat).toFinset

def n57W1 : Finset Nat :=
  ([0, 2, 14, 21, 29, 32, 45, 49, 54, 55] : List Nat).toFinset

def n57W2 : Finset Nat :=
  ([1, 2, 7, 11, 24, 27, 35, 42, 54, 56] : List Nat).toFinset

def n57W3 : Finset Nat :=
  ([1, 3, 15, 22, 30, 33, 46, 50, 55, 56] : List Nat).toFinset

def n57W4 : Finset Nat :=
  ([2, 3, 8, 12, 25, 28, 36, 43, 55, 57] : List Nat).toFinset

def n57W5 : Finset Nat :=
  ([2, 4, 16, 23, 31, 34, 47, 51, 56, 57] : List Nat).toFinset

def n57Face : Finset (Finset Nat) :=
  ([n57W0, n57W1, n57W2, n57W3, n57W4, n57W5] : List (Finset Nat)).toFinset

def n57FaceIndices : Finset Nat :=
  Finset.range 6

def n57IndexList : List Nat :=
  List.range 6

def n57WitnessOfIndex : Nat -> Finset Nat

  | 0 => n57W0

  | 1 => n57W1

  | 2 => n57W2

  | 3 => n57W3

  | 4 => n57W4

  | 5 => n57W5

  | _ => ∅

theorem n57_indexed_face_matches_exported_face :
    (n57IndexList.map n57WitnessOfIndex).toFinset = n57Face := by
  native_decide

theorem n57_exported_face_card :
    n57Face.card = 6 := by
  native_decide

theorem n57_face_indices_card :
    n57FaceIndices.card = 6 := by
  native_decide

theorem n57_all_exported_witnesses_have_card_h :
    n57W0.card = 10 ∧
    n57W1.card = 10 ∧
    n57W2.card = 10 ∧
    n57W3.card = 10 ∧
    n57W4.card = 10 ∧
    n57W5.card = 10 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩

theorem n57_all_exported_witnesses_in_range :
    n57W0 ⊆ Finset.range 58 ∧
    n57W1 ⊆ Finset.range 58 ∧
    n57W2 ⊆ Finset.range 58 ∧
    n57W3 ⊆ Finset.range 58 ∧
    n57W4 ⊆ Finset.range 58 ∧
    n57W5 ⊆ Finset.range 58 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩

theorem n57_all_exported_witnesses_are_sidon :
    Erdos.Sidon.IsSidonSet n57W0 ∧
    Erdos.Sidon.IsSidonSet n57W1 ∧
    Erdos.Sidon.IsSidonSet n57W2 ∧
    Erdos.Sidon.IsSidonSet n57W3 ∧
    Erdos.Sidon.IsSidonSet n57W4 ∧
    Erdos.Sidon.IsSidonSet n57W5 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩

def n57PrefixRank : Nat -> Nat
  | 0 => 2
  | 1 => 2
  | 2 => 1
  | 3 => 1
  | 4 => 0
  | 5 => 0
  | _ => 99

def n57MassRank : Nat -> Nat
  | 0 => 5
  | 1 => 2
  | 2 => 4
  | 3 => 0
  | 4 => 3
  | 5 => 1
  | _ => 99

def n57JointRank : Nat -> Nat
  | 0 => 5
  | 1 => 2
  | 2 => 4
  | 3 => 1
  | 4 => 3
  | 5 => 0
  | _ => 99

def n57ExactMassTwice (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.densityAdjustedMassTwice 57 (n57WitnessOfIndex i)

def n57MassWinnerIndicesByExactMass : List Nat :=
  n57IndexList.filter (fun i => n57ExactMassTwice i == 5)

def n57PacketPrefixLocation : Nat -> Nat
  | 0 => 55
  | 1 => 55
  | 2 => 56
  | 3 => 56
  | 4 => 57
  | 5 => 57
  | _ => 0

def n57PacketPrefixCount : Nat -> Nat
  | 0 => 10
  | 1 => 10
  | 2 => 10
  | 3 => 10
  | 4 => 10
  | 5 => 10
  | _ => 0

noncomputable def n57ExactPrefixProbe (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.prefixResidualProbeCard10 57 (n57PacketPrefixLocation i) (n57PacketPrefixCount i)

noncomputable def n57FullPrefixResidualMax (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.fullPrefixResidualMax 57 (n57WitnessOfIndex i)


noncomputable def n57ExactJointKey (i : Nat) : ℝ :=
  n57ExactPrefixProbe i * Real.sqrt (57 : ℝ) +
    (n57ExactMassTwice i : ℝ) / 2


noncomputable def n57ScalarFullPrefixJointKey (i : Nat) : ℝ :=
  n57FullPrefixResidualMax i * Real.sqrt (57 : ℝ) +
    (n57ExactMassTwice i : ℝ) / 2


def n57PrefixCountAtPacketLocation (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.prefixCount (n57WitnessOfIndex i) (n57PacketPrefixLocation i)

def n57PrefixWinnerIndicesByRank : List Nat :=
  n57IndexList.filter (fun i => n57PrefixRank i == 0)

def n57MassWinnerIndicesByRank : List Nat :=
  n57IndexList.filter (fun i => n57MassRank i == 0)

def n57JointWinnerIndicesByRank : List Nat :=
  n57IndexList.filter (fun i => n57JointRank i == 0)

def n57DominatesPrefixMass (i j : Nat) : Bool :=
  (n57PrefixRank i <= n57PrefixRank j) &&
  (n57MassRank i <= n57MassRank j) &&
  ((n57PrefixRank i < n57PrefixRank j) ||
   (n57MassRank i < n57MassRank j))

def n57ParetoIndicesByRanks : List Nat :=
  n57IndexList.filter (fun i => !(n57IndexList.any (fun j => n57DominatesPrefixMass j i)))

theorem n57_prefix_winner_match_packet :
    n57PrefixWinnerIndicesByRank = [4, 5] := by
  native_decide

theorem n57_mass_winner_match_packet :
    n57MassWinnerIndicesByRank = [3] := by
  native_decide

theorem n57_joint_winner_match_packet :
    n57JointWinnerIndicesByRank = [5] := by
  native_decide

theorem n57_pareto_minimal_match_packet :
    n57ParetoIndicesByRanks = [3, 5] := by
  native_decide

theorem n57_exact_mass_twice_table :
    n57IndexList.map n57ExactMassTwice = [129, 25, 109, 5, 89, 15] := by
  native_decide

theorem n57_exact_mass_winner_match_packet :
    n57MassWinnerIndicesByExactMass = [3] := by
  native_decide

theorem n57_exact_mass_winner_matches_rank_winner :
    n57MassWinnerIndicesByExactMass = n57MassWinnerIndicesByRank := by
  native_decide

theorem n57_packet_prefix_location_table :
    n57IndexList.map n57PacketPrefixLocation = [55, 55, 56, 56, 57, 57] := by
  native_decide

theorem n57_prefix_count_at_packet_location_table :
    n57IndexList.map n57PrefixCountAtPacketLocation = [10, 10, 10, 10, 10, 10] := by
  native_decide

lemma sqrt57_ge_prefix_lower : (15/2 : ℝ) ≤ Real.sqrt (57 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 15/2) (by norm_num : (0:ℝ) ≤ 57)]
  norm_num

lemma sqrt57_lt_prefix_upper : Real.sqrt (57 : ℝ) < (8 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 8)]
  norm_num

lemma sqrt57_lt_prefix_tight_upper : Real.sqrt (57 : ℝ) < (23/3 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
  norm_num

def n57W4PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 1
  | 3 => 2
  | 4 => 2
  | 5 => 2
  | 6 => 2
  | 7 => 2
  | 8 => 3
  | 9 => 3
  | 10 => 3
  | 11 => 3
  | 12 => 4
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
  | 25 => 5
  | 26 => 5
  | 27 => 5
  | 28 => 6
  | 29 => 6
  | 30 => 6
  | 31 => 6
  | 32 => 6
  | 33 => 6
  | 34 => 6
  | 35 => 6
  | 36 => 7
  | 37 => 7
  | 38 => 7
  | 39 => 7
  | 40 => 7
  | 41 => 7
  | 42 => 7
  | 43 => 8
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
  | 55 => 9
  | 56 => 9
  | 57 => 10
  | _ => 10
theorem n57W4_prefix_count_table_matches_witness :
    (List.range 58).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n57W4 t) =
      [0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

theorem n57W4_prefix_count_model_table :
    (List.range 58).map n57W4PrefixCountTable =
      [0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

lemma n57W4_prefix_segment_0_1_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_2_2_bound :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_2_2 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_3_7_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_3_7 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_8_11_bound :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (11 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_8_11 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_12_24_bound :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (12 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (24 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_12_24 :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_25_27_bound :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (25 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (27 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_25_27 :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_28_35_bound :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (28 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_28_35 :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_36_42_bound :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_36_42 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_43_54_bound :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_43_54 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_55_56_bound :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_55_56 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W4_prefix_segment_57_57_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W4_prefix_count_segment_57_57 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W4 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W4_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) := by
  exact ⟨n57W4_prefix_segment_0_1_bound, ⟨n57W4_prefix_segment_2_2_bound, ⟨n57W4_prefix_segment_3_7_bound, ⟨n57W4_prefix_segment_8_11_bound, ⟨n57W4_prefix_segment_12_24_bound, ⟨n57W4_prefix_segment_25_27_bound, ⟨n57W4_prefix_segment_28_35_bound, ⟨n57W4_prefix_segment_36_42_bound, ⟨n57W4_prefix_segment_43_54_bound, ⟨n57W4_prefix_segment_55_56_bound, n57W4_prefix_segment_57_57_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n57W4_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W4 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W4.card = 10 := by
    native_decide
  have hcard : (n57W4.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]

theorem n57W4_prefix_residual_segment_0_1_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_0_1_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_2_2_zero :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_2_2_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_2_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_3_7_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_3_7_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_3_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_8_11_zero :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_8_11_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_8_11 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_12_24_zero :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_12_24_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_12_24 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_25_27_zero :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_25_27_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_25_27 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_28_35_zero :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_28_35_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_28_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_36_42_zero :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_36_42_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_36_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_43_54_zero :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_43_54_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_43_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_55_56_zero :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_55_56_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_55_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_prefix_residual_segment_57_57_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W4_prefix_segment_57_57_bound t hlow hhigh
  have hcount := n57W4_prefix_count_segment_57_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W4_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W4 := by
  intro t ht
  interval_cases t
  · exact n57W4_prefix_residual_segment_0_1_zero 0 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_0_1_zero 1 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_2_2_zero 2 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_3_7_zero 3 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_3_7_zero 4 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_3_7_zero 5 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_3_7_zero 6 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_3_7_zero 7 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_8_11_zero 8 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_8_11_zero 9 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_8_11_zero 10 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_8_11_zero 11 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 12 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 13 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 14 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 15 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 16 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 17 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 18 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 19 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 20 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 21 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 22 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 23 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_12_24_zero 24 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_25_27_zero 25 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_25_27_zero 26 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_25_27_zero 27 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 28 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 29 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 30 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 31 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 32 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 33 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 34 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_28_35_zero 35 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 36 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 37 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 38 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 39 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 40 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 41 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_36_42_zero 42 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 43 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 44 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 45 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 46 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 47 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 48 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 49 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 50 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 51 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 52 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 53 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_43_54_zero 54 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_55_56_zero 55 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_55_56_zero 56 (by norm_num) (by norm_num)
  · exact n57W4_prefix_residual_segment_57_57_zero 57 (by norm_num) (by norm_num)

def n57W5PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 1
  | 3 => 1
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
  | 22 => 3
  | 23 => 4
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 4
  | 30 => 4
  | 31 => 5
  | 32 => 5
  | 33 => 5
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
  | 57 => 10
  | _ => 10
theorem n57W5_prefix_count_table_matches_witness :
    (List.range 58).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n57W5 t) =
      [0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

theorem n57W5_prefix_count_model_table :
    (List.range 58).map n57W5PrefixCountTable =
      [0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

lemma n57W5_prefix_segment_0_1_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_2_3_bound :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_2_3 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_4_15_bound :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (15 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_4_15 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_16_22_bound :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (16 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_16_22 :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_23_30_bound :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (30 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_23_30 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_31_33_bound :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (31 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_31_33 :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_34_46_bound :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (46 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_34_46 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_47_50_bound :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (47 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (50 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_47_50 :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_51_55_bound :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (51 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_51_55 :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_56_56_bound :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_56_56 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n57W5_prefix_segment_57_57_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper, hlowR, hhighR]

theorem n57W5_prefix_count_segment_57_57 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W5 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W5_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        10 * Real.sqrt (57 : ℝ) - 57) := by
  exact ⟨n57W5_prefix_segment_0_1_bound, ⟨n57W5_prefix_segment_2_3_bound, ⟨n57W5_prefix_segment_4_15_bound, ⟨n57W5_prefix_segment_16_22_bound, ⟨n57W5_prefix_segment_23_30_bound, ⟨n57W5_prefix_segment_31_33_bound, ⟨n57W5_prefix_segment_34_46_bound, ⟨n57W5_prefix_segment_47_50_bound, ⟨n57W5_prefix_segment_51_55_bound, ⟨n57W5_prefix_segment_56_56_bound, n57W5_prefix_segment_57_57_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n57W5_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W5 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W5.card = 10 := by
    native_decide
  have hcard : (n57W5.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]

theorem n57W5_prefix_residual_segment_0_1_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_0_1_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_2_3_zero :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_2_3_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_2_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_4_15_zero :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_4_15_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_4_15 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_16_22_zero :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_16_22_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_16_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_23_30_zero :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_23_30_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_23_30 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_31_33_zero :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_31_33_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_31_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_34_46_zero :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_34_46_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_34_46 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_47_50_zero :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_47_50_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_47_50 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_51_55_zero :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_51_55_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_51_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_56_56_zero :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_56_56_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_56_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_prefix_residual_segment_57_57_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 t = 0 := by
  intro t hlow hhigh
  have hdev := n57W5_prefix_segment_57_57_bound t hlow hhigh
  have hcount := n57W5_prefix_count_segment_57_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n57W5_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W5 := by
  intro t ht
  interval_cases t
  · exact n57W5_prefix_residual_segment_0_1_zero 0 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_0_1_zero 1 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_2_3_zero 2 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_2_3_zero 3 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 4 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 5 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 6 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 7 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 8 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 9 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 10 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 11 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 12 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 13 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 14 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_4_15_zero 15 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 16 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 17 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 18 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 19 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 20 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 21 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_16_22_zero 22 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 23 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 24 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 25 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 26 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 27 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 28 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 29 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_23_30_zero 30 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_31_33_zero 31 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_31_33_zero 32 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_31_33_zero 33 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 34 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 35 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 36 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 37 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 38 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 39 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 40 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 41 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 42 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 43 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 44 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 45 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_34_46_zero 46 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_47_50_zero 47 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_47_50_zero 48 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_47_50_zero 49 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_47_50_zero 50 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_51_55_zero 51 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_51_55_zero 52 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_51_55_zero 53 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_51_55_zero 54 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_51_55_zero 55 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_56_56_zero 56 (by norm_num) (by norm_num)
  · exact n57W5_prefix_residual_segment_57_57_zero 57 (by norm_num) (by norm_num)


lemma sqrt57_ge_57_div_10 : (57/10 : ℝ) ≤ Real.sqrt (57 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 57/10) (by norm_num : (0:ℝ) ≤ 57)]
  norm_num


theorem n57W0_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W0 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W0.card = 10 := by
    native_decide
  have hcard : (n57W0.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]


theorem n57W1_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W1 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W1.card = 10 := by
    native_decide
  have hcard : (n57W1.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]


theorem n57W2_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W2 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W2.card = 10 := by
    native_decide
  have hcard : (n57W2.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]


theorem n57W3_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 57 n57W3 =
      10 * Real.sqrt (57 : ℝ) - 57 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n57W3.card = 10 := by
    native_decide
  have hcard : (n57W3.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (57 : ℝ) := by
    nlinarith [sqrt57_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 57)]
  nlinarith [hsq]


theorem n57_exact_prefix_probe_4_zero :
    n57ExactPrefixProbe 4 = 0 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n57_exact_prefix_probe_5_zero :
    n57ExactPrefixProbe 5 = 0 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n57_exact_prefix_probe_0_positive :
    0 < n57ExactPrefixProbe 0 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_prefix_probe_1_positive :
    0 < n57ExactPrefixProbe 1 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_prefix_probe_2_positive :
    0 < n57ExactPrefixProbe 2 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_prefix_probe_3_positive :
    0 < n57ExactPrefixProbe 3 := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_prefix_probe_positive_on_nonwinners :
    0 < n57ExactPrefixProbe 0 ∧
    0 < n57ExactPrefixProbe 1 ∧
    0 < n57ExactPrefixProbe 2 ∧
    0 < n57ExactPrefixProbe 3 := by
  exact ⟨n57_exact_prefix_probe_0_positive, ⟨n57_exact_prefix_probe_1_positive, ⟨n57_exact_prefix_probe_2_positive, n57_exact_prefix_probe_3_positive⟩⟩⟩


theorem n57_exact_prefix_probe_0_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 0 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 (n57PacketPrefixLocation 0) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W0 (n57PacketPrefixLocation 0) =
        n57PacketPrefixCount 0 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_exact_prefix_probe_1_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 1 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 (n57PacketPrefixLocation 1) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W1 (n57PacketPrefixLocation 1) =
        n57PacketPrefixCount 1 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_exact_prefix_probe_2_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 2 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 (n57PacketPrefixLocation 2) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W2 (n57PacketPrefixLocation 2) =
        n57PacketPrefixCount 2 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_exact_prefix_probe_3_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 3 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 (n57PacketPrefixLocation 3) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W3 (n57PacketPrefixLocation 3) =
        n57PacketPrefixCount 3 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_exact_prefix_probe_4_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 4 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W4 (n57PacketPrefixLocation 4) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W4 (n57PacketPrefixLocation 4) =
        n57PacketPrefixCount 4 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_exact_prefix_probe_5_matches_prefix_residual_at_packet_location :
    n57ExactPrefixProbe 5 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W5 (n57PacketPrefixLocation 5) := by
  unfold n57ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n57W5 (n57PacketPrefixLocation 5) =
        n57PacketPrefixCount 5 := by
    native_decide
  rw [hcount]
  rfl


theorem n57_prefix_residual_at_packet_location_0_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 (n57PacketPrefixLocation 0) := by
  rw [← n57_exact_prefix_probe_0_matches_prefix_residual_at_packet_location]
  exact n57_exact_prefix_probe_0_positive


theorem n57_prefix_residual_at_packet_location_1_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 (n57PacketPrefixLocation 1) := by
  rw [← n57_exact_prefix_probe_1_matches_prefix_residual_at_packet_location]
  exact n57_exact_prefix_probe_1_positive


theorem n57_prefix_residual_at_packet_location_2_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 (n57PacketPrefixLocation 2) := by
  rw [← n57_exact_prefix_probe_2_matches_prefix_residual_at_packet_location]
  exact n57_exact_prefix_probe_2_positive


theorem n57_prefix_residual_at_packet_location_3_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 (n57PacketPrefixLocation 3) := by
  rw [← n57_exact_prefix_probe_3_matches_prefix_residual_at_packet_location]
  exact n57_exact_prefix_probe_3_positive


theorem n57_full_prefix_residual_zero_on_prefix_winners :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W4 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W5 := by
  exact ⟨n57W4_full_prefix_residual_zero, n57W5_full_prefix_residual_zero⟩


theorem n57_actual_prefix_residual_positive_on_probe_nonwinners :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 (n57PacketPrefixLocation 0) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 (n57PacketPrefixLocation 1) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 (n57PacketPrefixLocation 2) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 (n57PacketPrefixLocation 3) := by
  exact ⟨n57_prefix_residual_at_packet_location_0_positive, ⟨n57_prefix_residual_at_packet_location_1_positive, ⟨n57_prefix_residual_at_packet_location_2_positive, n57_prefix_residual_at_packet_location_3_positive⟩⟩⟩


theorem n57_full_prefix_residual_max_4_zero :
    n57FullPrefixResidualMax 4 = 0 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n57W4_full_prefix_residual_zero


theorem n57_full_prefix_residual_max_5_zero :
    n57FullPrefixResidualMax 5 = 0 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n57W5_full_prefix_residual_zero


theorem n57_full_prefix_residual_max_0_positive :
    0 < n57FullPrefixResidualMax 0 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n57PacketPrefixLocation 0 ≤ 57)
    n57_prefix_residual_at_packet_location_0_positive


theorem n57_full_prefix_residual_max_1_positive :
    0 < n57FullPrefixResidualMax 1 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n57PacketPrefixLocation 1 ≤ 57)
    n57_prefix_residual_at_packet_location_1_positive


theorem n57_full_prefix_residual_max_2_positive :
    0 < n57FullPrefixResidualMax 2 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n57PacketPrefixLocation 2 ≤ 57)
    n57_prefix_residual_at_packet_location_2_positive


theorem n57_full_prefix_residual_max_3_positive :
    0 < n57FullPrefixResidualMax 3 := by
  unfold n57FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n57PacketPrefixLocation 3 ≤ 57)
    n57_prefix_residual_at_packet_location_3_positive


theorem n57_full_prefix_residual_max_zero_on_prefix_winners :
    n57FullPrefixResidualMax 4 = 0 ∧
        n57FullPrefixResidualMax 5 = 0 := by
  exact ⟨n57_full_prefix_residual_max_4_zero, n57_full_prefix_residual_max_5_zero⟩


theorem n57_full_prefix_residual_max_positive_on_probe_nonwinners :
    0 < n57FullPrefixResidualMax 0 ∧
        0 < n57FullPrefixResidualMax 1 ∧
        0 < n57FullPrefixResidualMax 2 ∧
        0 < n57FullPrefixResidualMax 3 := by
  exact ⟨n57_full_prefix_residual_max_0_positive, ⟨n57_full_prefix_residual_max_1_positive, ⟨n57_full_prefix_residual_max_2_positive, n57_full_prefix_residual_max_3_positive⟩⟩⟩


theorem n57_exact_prefix_probe_0_value :
    n57ExactPrefixProbe 0 = (2 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_0_value :
    n57ExactMassTwice 0 = 129 := by
  native_decide


theorem n57_exact_prefix_probe_1_value :
    n57ExactPrefixProbe 1 = (2 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_1_value :
    n57ExactMassTwice 1 = 25 := by
  native_decide


theorem n57_exact_prefix_probe_2_value :
    n57ExactPrefixProbe 2 = (1 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_2_value :
    n57ExactMassTwice 2 = 109 := by
  native_decide


theorem n57_exact_prefix_probe_3_value :
    n57ExactPrefixProbe 3 = (1 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_3_value :
    n57ExactMassTwice 3 = 5 := by
  native_decide


theorem n57_exact_prefix_probe_4_value :
    n57ExactPrefixProbe 4 = (0 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_4_value :
    n57ExactMassTwice 4 = 89 := by
  native_decide


theorem n57_exact_prefix_probe_5_value :
    n57ExactPrefixProbe 5 = (0 : ℝ) := by
  simp [n57ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n57PacketPrefixLocation, n57PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (57 : ℝ) ≤ 0 := by
    nlinarith [sqrt57_ge_57_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n57_exact_mass_twice_5_value :
    n57ExactMassTwice 5 = 15 := by
  native_decide


theorem n57W0_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_1_5 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_1_5_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (5 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_1_5 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_6_9 :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_6_9_le_exact_probe :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (6 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_6_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_10_22 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_10_22_le_exact_probe :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_10_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_23_25 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_23_25_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_23_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_26_33 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_26_33_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_26_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_34_40 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_34_40_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (40 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_34_40 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_41_52 :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_41_52_le_exact_probe :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (41 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (52 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_41_52 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_53_54 :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_53_54_le_exact_probe :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (53 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_53_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_prefix_count_segment_55_57 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W0 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W0_prefix_residual_segment_55_57_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W0_prefix_count_segment_55_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_0_positive
  · rw [n57_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W0_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 t ≤
        n57ExactPrefixProbe 0 := by
  intro t ht
  interval_cases t
  · exact n57W0_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_1_5_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_1_5_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_1_5_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_1_5_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_1_5_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_6_9_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_6_9_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_6_9_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_6_9_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_10_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_23_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_23_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_23_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_26_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_34_40_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_41_52_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_53_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_53_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_55_57_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_55_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n57W0_prefix_residual_segment_55_57_le_exact_probe 57 (by norm_num) (by norm_num)

theorem n57W1_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_2_13 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_2_13_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_2_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_14_20 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_14_20_le_exact_probe :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_14_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_21_28 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_21_28_le_exact_probe :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_21_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_29_31 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_29_31_le_exact_probe :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_29_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_32_44 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_32_44_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_32_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_45_48 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_45_48_le_exact_probe :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_45_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_49_53 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_49_53_le_exact_probe :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_49_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_54_54 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_54_54_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_54_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_prefix_count_segment_55_57 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W1 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W1_prefix_residual_segment_55_57_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W1_prefix_count_segment_55_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_1_positive
  · rw [n57_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W1_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 t ≤
        n57ExactPrefixProbe 1 := by
  intro t ht
  interval_cases t
  · exact n57W1_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_2_13_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_14_20_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_21_28_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_29_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_29_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_29_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_32_44_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_45_48_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_45_48_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_45_48_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_45_48_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_49_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_49_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_49_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_49_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_49_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_54_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_55_57_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_55_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n57W1_prefix_residual_segment_55_57_le_exact_probe 57 (by norm_num) (by norm_num)

theorem n57W2_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_1_1 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_1_1_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_1_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_2_6 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_2_6_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (6 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_2_6 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_7_10 :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_7_10_le_exact_probe :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (7 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (10 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_7_10 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_11_23 :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_11_23_le_exact_probe :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (11 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_11_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_24_26 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_24_26_le_exact_probe :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_24_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_27_34 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_27_34_le_exact_probe :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_27_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_35_41 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_35_41_le_exact_probe :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (41 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_35_41 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_42_53 :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_42_53_le_exact_probe :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (42 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_42_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_54_55 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_54_55_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_54_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W2 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W2_prefix_residual_segment_56_57_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W2_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_2_positive
  · rw [n57_exact_prefix_probe_2_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W2_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 t ≤
        n57ExactPrefixProbe 2 := by
  intro t ht
  interval_cases t
  · exact n57W2_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_1_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_2_6_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_2_6_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_2_6_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_2_6_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_2_6_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_7_10_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_7_10_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_7_10_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_7_10_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_11_23_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_24_26_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_24_26_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_24_26_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_27_34_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_35_41_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_42_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_54_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_54_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_56_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n57W2_prefix_residual_segment_56_57_le_exact_probe 57 (by norm_num) (by norm_num)

theorem n57W3_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_1_2 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_1_2_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_1_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_3_14 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_3_14_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_3_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_15_21 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_15_21_le_exact_probe :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_15_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_22_29 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_22_29_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_22_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_30_32 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_30_32_le_exact_probe :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_30_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_33_45 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_33_45_le_exact_probe :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_33_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_46_49_le_exact_probe :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_50_54_le_exact_probe :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_55_55 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_55_55_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_55_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n57W3 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n57W3_prefix_residual_segment_56_57_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (57 : ℝ)| ≤
        (10 * Real.sqrt (57 : ℝ) - 57) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper,
        sqrt57_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n57W3_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n57W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n57_exact_prefix_probe_3_positive
  · rw [n57_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n57W3_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 t ≤
        n57ExactPrefixProbe 3 := by
  intro t ht
  interval_cases t
  · exact n57W3_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_1_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_1_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_3_14_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_15_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_22_29_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_30_32_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_30_32_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_30_32_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_33_45_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_46_49_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_46_49_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_46_49_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_46_49_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_50_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_50_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_50_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_50_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_50_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_55_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_56_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n57W3_prefix_residual_segment_56_57_le_exact_probe 57 (by norm_num) (by norm_num)


theorem n57_full_prefix_residual_max_0_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 0 = n57ExactPrefixProbe 0 := by
  unfold n57FullPrefixResidualMax n57WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n57PacketPrefixLocation 0 ≤ 57)
    n57W0_full_prefix_residual_le_exact_probe
    (n57_exact_prefix_probe_0_matches_prefix_residual_at_packet_location).symm


theorem n57_full_prefix_residual_max_1_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 1 = n57ExactPrefixProbe 1 := by
  unfold n57FullPrefixResidualMax n57WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n57PacketPrefixLocation 1 ≤ 57)
    n57W1_full_prefix_residual_le_exact_probe
    (n57_exact_prefix_probe_1_matches_prefix_residual_at_packet_location).symm


theorem n57_full_prefix_residual_max_2_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 2 = n57ExactPrefixProbe 2 := by
  unfold n57FullPrefixResidualMax n57WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n57PacketPrefixLocation 2 ≤ 57)
    n57W2_full_prefix_residual_le_exact_probe
    (n57_exact_prefix_probe_2_matches_prefix_residual_at_packet_location).symm


theorem n57_full_prefix_residual_max_3_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 3 = n57ExactPrefixProbe 3 := by
  unfold n57FullPrefixResidualMax n57WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n57PacketPrefixLocation 3 ≤ 57)
    n57W3_full_prefix_residual_le_exact_probe
    (n57_exact_prefix_probe_3_matches_prefix_residual_at_packet_location).symm


theorem n57_full_prefix_residual_max_4_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 4 = n57ExactPrefixProbe 4 := by
  rw [n57_full_prefix_residual_max_4_zero, n57_exact_prefix_probe_4_zero]


theorem n57_full_prefix_residual_max_5_matches_exact_prefix_probe :
    n57FullPrefixResidualMax 5 = n57ExactPrefixProbe 5 := by
  rw [n57_full_prefix_residual_max_5_zero, n57_exact_prefix_probe_5_zero]


theorem n57_full_prefix_residual_max_matches_exact_prefix_probe_on_face :
    ∀ i ∈ n57FaceIndices, n57FullPrefixResidualMax i = n57ExactPrefixProbe i := by
  intro i hi
  simp [n57FaceIndices] at hi
  interval_cases i
  · exact n57_full_prefix_residual_max_0_matches_exact_prefix_probe
  · exact n57_full_prefix_residual_max_1_matches_exact_prefix_probe
  · exact n57_full_prefix_residual_max_2_matches_exact_prefix_probe
  · exact n57_full_prefix_residual_max_3_matches_exact_prefix_probe
  · exact n57_full_prefix_residual_max_4_matches_exact_prefix_probe
  · exact n57_full_prefix_residual_max_5_matches_exact_prefix_probe


theorem n57_exact_joint_key_0_value :
    n57ExactJointKey 0 = (2 : ℝ) * Real.sqrt (57 : ℝ) + (129/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_0_value, n57_exact_mass_twice_0_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_0_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 0 = n57ExactJointKey 0 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_0_matches_exact_prefix_probe]


theorem n57_exact_joint_key_1_value :
    n57ExactJointKey 1 = (2 : ℝ) * Real.sqrt (57 : ℝ) + (25/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_1_value, n57_exact_mass_twice_1_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_1_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 1 = n57ExactJointKey 1 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_1_matches_exact_prefix_probe]


theorem n57_exact_joint_key_2_value :
    n57ExactJointKey 2 = (1 : ℝ) * Real.sqrt (57 : ℝ) + (109/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_2_value, n57_exact_mass_twice_2_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_2_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 2 = n57ExactJointKey 2 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_2_matches_exact_prefix_probe]


theorem n57_exact_joint_key_3_value :
    n57ExactJointKey 3 = (1 : ℝ) * Real.sqrt (57 : ℝ) + (5/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_3_value, n57_exact_mass_twice_3_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_3_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 3 = n57ExactJointKey 3 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_3_matches_exact_prefix_probe]


theorem n57_exact_joint_key_4_value :
    n57ExactJointKey 4 = (89/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_4_value, n57_exact_mass_twice_4_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_4_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 4 = n57ExactJointKey 4 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_4_matches_exact_prefix_probe]


theorem n57_exact_joint_key_5_value :
    n57ExactJointKey 5 = (15/2 : ℝ) := by
  unfold n57ExactJointKey
  rw [n57_exact_prefix_probe_5_value, n57_exact_mass_twice_5_value]
  have hsq : (Real.sqrt (57 : ℝ)) ^ 2 = (57 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 57)]
  nlinarith


theorem n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key :
    n57ScalarFullPrefixJointKey 5 = n57ExactJointKey 5 := by
  unfold n57ScalarFullPrefixJointKey n57ExactJointKey
  rw [n57_full_prefix_residual_max_5_matches_exact_prefix_probe]


theorem n57_exact_joint_key_5_lt_0 :
    n57ExactJointKey 5 < n57ExactJointKey 0 := by
  rw [n57_exact_joint_key_5_value, n57_exact_joint_key_0_value]
  nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper]


theorem n57_scalar_full_prefix_joint_key_5_lt_0 :
    n57ScalarFullPrefixJointKey 5 <
      n57ScalarFullPrefixJointKey 0 := by
  rw [n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key,
    n57_scalar_full_prefix_joint_key_0_matches_exact_joint_key]
  exact n57_exact_joint_key_5_lt_0


theorem n57_exact_joint_key_5_lt_1 :
    n57ExactJointKey 5 < n57ExactJointKey 1 := by
  rw [n57_exact_joint_key_5_value, n57_exact_joint_key_1_value]
  nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper]


theorem n57_scalar_full_prefix_joint_key_5_lt_1 :
    n57ScalarFullPrefixJointKey 5 <
      n57ScalarFullPrefixJointKey 1 := by
  rw [n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key,
    n57_scalar_full_prefix_joint_key_1_matches_exact_joint_key]
  exact n57_exact_joint_key_5_lt_1


theorem n57_exact_joint_key_5_lt_2 :
    n57ExactJointKey 5 < n57ExactJointKey 2 := by
  rw [n57_exact_joint_key_5_value, n57_exact_joint_key_2_value]
  nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper]


theorem n57_scalar_full_prefix_joint_key_5_lt_2 :
    n57ScalarFullPrefixJointKey 5 <
      n57ScalarFullPrefixJointKey 2 := by
  rw [n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key,
    n57_scalar_full_prefix_joint_key_2_matches_exact_joint_key]
  exact n57_exact_joint_key_5_lt_2


theorem n57_exact_joint_key_5_lt_3 :
    n57ExactJointKey 5 < n57ExactJointKey 3 := by
  rw [n57_exact_joint_key_5_value, n57_exact_joint_key_3_value]
  nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper]


theorem n57_scalar_full_prefix_joint_key_5_lt_3 :
    n57ScalarFullPrefixJointKey 5 <
      n57ScalarFullPrefixJointKey 3 := by
  rw [n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key,
    n57_scalar_full_prefix_joint_key_3_matches_exact_joint_key]
  exact n57_exact_joint_key_5_lt_3


theorem n57_exact_joint_key_5_lt_4 :
    n57ExactJointKey 5 < n57ExactJointKey 4 := by
  rw [n57_exact_joint_key_5_value, n57_exact_joint_key_4_value]
  nlinarith [sqrt57_ge_prefix_lower, sqrt57_lt_prefix_upper]


theorem n57_scalar_full_prefix_joint_key_5_lt_4 :
    n57ScalarFullPrefixJointKey 5 <
      n57ScalarFullPrefixJointKey 4 := by
  rw [n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key,
    n57_scalar_full_prefix_joint_key_4_matches_exact_joint_key]
  exact n57_exact_joint_key_5_lt_4


theorem n57_exact_joint_winner_strict_certificate :
    n57ExactJointKey 5 < n57ExactJointKey 0 ∧
    n57ExactJointKey 5 < n57ExactJointKey 1 ∧
    n57ExactJointKey 5 < n57ExactJointKey 2 ∧
    n57ExactJointKey 5 < n57ExactJointKey 3 ∧
    n57ExactJointKey 5 < n57ExactJointKey 4 := by
  exact ⟨n57_exact_joint_key_5_lt_0, ⟨n57_exact_joint_key_5_lt_1, ⟨n57_exact_joint_key_5_lt_2, ⟨n57_exact_joint_key_5_lt_3, n57_exact_joint_key_5_lt_4⟩⟩⟩⟩


theorem n57_scalar_full_prefix_joint_winner_strict_certificate :
    n57ScalarFullPrefixJointKey 5 < n57ScalarFullPrefixJointKey 0 ∧
    n57ScalarFullPrefixJointKey 5 < n57ScalarFullPrefixJointKey 1 ∧
    n57ScalarFullPrefixJointKey 5 < n57ScalarFullPrefixJointKey 2 ∧
    n57ScalarFullPrefixJointKey 5 < n57ScalarFullPrefixJointKey 3 ∧
    n57ScalarFullPrefixJointKey 5 < n57ScalarFullPrefixJointKey 4 := by
  exact ⟨n57_scalar_full_prefix_joint_key_5_lt_0, ⟨n57_scalar_full_prefix_joint_key_5_lt_1, ⟨n57_scalar_full_prefix_joint_key_5_lt_2, ⟨n57_scalar_full_prefix_joint_key_5_lt_3, n57_scalar_full_prefix_joint_key_5_lt_4⟩⟩⟩⟩


theorem n57_prefix_min_4 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57PrefixRank 4 := by
  constructor
  · native_decide
  · intro y hy
    simp [n57FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n57_mass_min_3 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57MassRank 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n57FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n57_joint_min_5 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57JointRank 5 := by
  constructor
  · native_decide
  · intro y hy
    simp [n57FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n57_exact_mass_min_3 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ExactMassTwice 3 := by
  constructor
  · native_decide
  · intro y hy
    simp [n57FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n57_exact_prefix_probe_min_4 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ExactPrefixProbe 4 := by
  constructor
  · native_decide
  · intro y hy
    rw [n57_exact_prefix_probe_4_zero]
    unfold n57ExactPrefixProbe Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
    exact le_max_left 0 _


theorem n57_full_prefix_residual_max_min_4 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57FullPrefixResidualMax 4 := by
  constructor
  · native_decide
  · intro y hy
    rw [n57_full_prefix_residual_max_4_zero]
    unfold n57FullPrefixResidualMax
    exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_nonneg 57 (n57WitnessOfIndex y)


theorem n57_exact_joint_min_5 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ExactJointKey 5 := by
  constructor
  · native_decide
  · intro y hy
    simp [n57FaceIndices] at hy
    interval_cases y
    · exact le_of_lt n57_exact_joint_key_5_lt_0
    · exact le_of_lt n57_exact_joint_key_5_lt_1
    · exact le_of_lt n57_exact_joint_key_5_lt_2
    · exact le_of_lt n57_exact_joint_key_5_lt_3
    · exact le_of_lt n57_exact_joint_key_5_lt_4
    · exact le_rfl


theorem n57_scalar_full_prefix_joint_min_5 :
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ScalarFullPrefixJointKey 5 := by
  exact Erdos.Collider.isFieldMinOn_weightedJoint_of_prefix_eq_on
    (F := n57FaceIndices)
    (x := 5)
    (prefixScalar := n57FullPrefixResidualMax)
    (prefixProbe := n57ExactPrefixProbe)
    (mass := fun i => (n57ExactMassTwice i : ℝ) / 2)
    (weight := Real.sqrt (57 : ℝ))
    n57_exact_joint_min_5
    n57_full_prefix_residual_max_matches_exact_prefix_probe_on_face


theorem n57_prefix_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57PrefixRank n57MassRank := by
  refine ⟨4, 3, n57_prefix_min_4, n57_mass_min_3, ?_⟩
  native_decide

theorem n57_prefix_mass_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_prefix_mass_field_split_full_exported_face


theorem n57_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57PrefixRank n57ExactMassTwice := by
  refine ⟨4, 3, n57_prefix_min_4, n57_exact_mass_min_3, ?_⟩
  native_decide

theorem n57_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_prefix_exact_mass_field_split_full_exported_face


theorem n57_exact_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57ExactPrefixProbe n57ExactMassTwice := by
  refine ⟨4, 3, n57_exact_prefix_probe_min_4, n57_exact_mass_min_3, ?_⟩
  native_decide

theorem n57_exact_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_exact_prefix_exact_mass_field_split_full_exported_face


theorem n57_full_prefix_residual_max_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57FullPrefixResidualMax n57ExactMassTwice := by
  refine ⟨4, 3, n57_full_prefix_residual_max_min_4, n57_exact_mass_min_3, ?_⟩
  native_decide

theorem n57_full_prefix_residual_max_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_full_prefix_residual_max_exact_mass_field_split_full_exported_face


theorem n57_mass_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57MassRank n57JointRank := by
  refine ⟨3, 5, n57_mass_min_3, n57_joint_min_5, ?_⟩
  native_decide

theorem n57_mass_joint_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_mass_joint_field_split_full_exported_face


theorem n57_exact_mass_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57ExactMassTwice n57JointRank := by
  refine ⟨3, 5, n57_exact_mass_min_3, n57_joint_min_5, ?_⟩
  native_decide

theorem n57_exact_mass_joint_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_exact_mass_joint_field_split_full_exported_face


theorem n57_prefix_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n57FaceIndices
      n57PrefixRank n57JointRank := by
  refine ⟨4, 5, n57_prefix_min_4, n57_joint_min_5, ?_⟩
  native_decide

theorem n57_prefix_joint_split_gives_two_exposed_indices :
    2 ≤ n57FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n57_prefix_joint_field_split_full_exported_face

/-! ## n = 58: complete exported face -/

def n58W0 : Finset Nat :=
  ([0, 1, 6, 10, 23, 26, 34, 41, 53, 55] : List Nat).toFinset

def n58W1 : Finset Nat :=
  ([0, 2, 14, 21, 29, 32, 45, 49, 54, 55] : List Nat).toFinset

def n58W2 : Finset Nat :=
  ([0, 2, 15, 21, 22, 32, 46, 50, 55, 58] : List Nat).toFinset

def n58W3 : Finset Nat :=
  ([0, 3, 8, 12, 26, 36, 37, 43, 56, 58] : List Nat).toFinset

def n58W4 : Finset Nat :=
  ([1, 2, 7, 11, 24, 27, 35, 42, 54, 56] : List Nat).toFinset

def n58W5 : Finset Nat :=
  ([1, 3, 15, 22, 30, 33, 46, 50, 55, 56] : List Nat).toFinset

def n58W6 : Finset Nat :=
  ([2, 3, 8, 12, 25, 28, 36, 43, 55, 57] : List Nat).toFinset

def n58W7 : Finset Nat :=
  ([2, 4, 16, 23, 31, 34, 47, 51, 56, 57] : List Nat).toFinset

def n58W8 : Finset Nat :=
  ([3, 4, 9, 13, 26, 29, 37, 44, 56, 58] : List Nat).toFinset

def n58W9 : Finset Nat :=
  ([3, 5, 17, 24, 32, 35, 48, 52, 57, 58] : List Nat).toFinset

def n58Face : Finset (Finset Nat) :=
  ([n58W0, n58W1, n58W2, n58W3, n58W4, n58W5, n58W6, n58W7, n58W8, n58W9] : List (Finset Nat)).toFinset

def n58FaceIndices : Finset Nat :=
  Finset.range 10

def n58IndexList : List Nat :=
  List.range 10

def n58WitnessOfIndex : Nat -> Finset Nat

  | 0 => n58W0

  | 1 => n58W1

  | 2 => n58W2

  | 3 => n58W3

  | 4 => n58W4

  | 5 => n58W5

  | 6 => n58W6

  | 7 => n58W7

  | 8 => n58W8

  | 9 => n58W9

  | _ => ∅

theorem n58_indexed_face_matches_exported_face :
    (n58IndexList.map n58WitnessOfIndex).toFinset = n58Face := by
  native_decide

theorem n58_exported_face_card :
    n58Face.card = 10 := by
  native_decide

theorem n58_face_indices_card :
    n58FaceIndices.card = 10 := by
  native_decide

theorem n58_all_exported_witnesses_have_card_h :
    n58W0.card = 10 ∧
    n58W1.card = 10 ∧
    n58W2.card = 10 ∧
    n58W3.card = 10 ∧
    n58W4.card = 10 ∧
    n58W5.card = 10 ∧
    n58W6.card = 10 ∧
    n58W7.card = 10 ∧
    n58W8.card = 10 ∧
    n58W9.card = 10 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n58_all_exported_witnesses_in_range :
    n58W0 ⊆ Finset.range 59 ∧
    n58W1 ⊆ Finset.range 59 ∧
    n58W2 ⊆ Finset.range 59 ∧
    n58W3 ⊆ Finset.range 59 ∧
    n58W4 ⊆ Finset.range 59 ∧
    n58W5 ⊆ Finset.range 59 ∧
    n58W6 ⊆ Finset.range 59 ∧
    n58W7 ⊆ Finset.range 59 ∧
    n58W8 ⊆ Finset.range 59 ∧
    n58W9 ⊆ Finset.range 59 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n58_all_exported_witnesses_are_sidon :
    Erdos.Sidon.IsSidonSet n58W0 ∧
    Erdos.Sidon.IsSidonSet n58W1 ∧
    Erdos.Sidon.IsSidonSet n58W2 ∧
    Erdos.Sidon.IsSidonSet n58W3 ∧
    Erdos.Sidon.IsSidonSet n58W4 ∧
    Erdos.Sidon.IsSidonSet n58W5 ∧
    Erdos.Sidon.IsSidonSet n58W6 ∧
    Erdos.Sidon.IsSidonSet n58W7 ∧
    Erdos.Sidon.IsSidonSet n58W8 ∧
    Erdos.Sidon.IsSidonSet n58W9 := by
  exact ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, ⟨by native_decide, by native_decide⟩⟩⟩⟩⟩⟩⟩⟩⟩

def n58PrefixRank : Nat -> Nat
  | 0 => 4
  | 1 => 4
  | 2 => 0
  | 3 => 1
  | 4 => 3
  | 5 => 3
  | 6 => 2
  | 7 => 2
  | 8 => 0
  | 9 => 0
  | _ => 99

def n58MassRank : Nat -> Nat
  | 0 => 7
  | 1 => 3
  | 2 => 3
  | 3 => 4
  | 4 => 6
  | 5 => 1
  | 6 => 5
  | 7 => 0
  | 8 => 4
  | 9 => 2
  | _ => 99

def n58JointRank : Nat -> Nat
  | 0 => 9
  | 1 => 5
  | 2 => 2
  | 3 => 6
  | 4 => 8
  | 5 => 3
  | 6 => 7
  | 7 => 0
  | 8 => 4
  | 9 => 1
  | _ => 99

def n58ExactMassTwice (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.densityAdjustedMassTwice 58 (n58WitnessOfIndex i)

def n58MassWinnerIndicesByExactMass : List Nat :=
  n58IndexList.filter (fun i => n58ExactMassTwice i == 4)

def n58PacketPrefixLocation : Nat -> Nat
  | 0 => 55
  | 1 => 55
  | 2 => 58
  | 3 => 12
  | 4 => 56
  | 5 => 56
  | 6 => 57
  | 7 => 57
  | 8 => 58
  | 9 => 58
  | _ => 0

def n58PacketPrefixCount : Nat -> Nat
  | 0 => 10
  | 1 => 10
  | 2 => 10
  | 3 => 4
  | 4 => 10
  | 5 => 10
  | 6 => 10
  | 7 => 10
  | 8 => 10
  | 9 => 10
  | _ => 0

noncomputable def n58ExactPrefixProbe (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.prefixResidualProbeCard10 58 (n58PacketPrefixLocation i) (n58PacketPrefixCount i)

noncomputable def n58FullPrefixResidualMax (i : Nat) : ℝ :=
  Erdos30FaceFieldExactObservables.fullPrefixResidualMax 58 (n58WitnessOfIndex i)


noncomputable def n58ExactJointKey (i : Nat) : ℝ :=
  n58ExactPrefixProbe i * Real.sqrt (58 : ℝ) +
    (n58ExactMassTwice i : ℝ) / 2


noncomputable def n58ScalarFullPrefixJointKey (i : Nat) : ℝ :=
  n58FullPrefixResidualMax i * Real.sqrt (58 : ℝ) +
    (n58ExactMassTwice i : ℝ) / 2


def n58PrefixCountAtPacketLocation (i : Nat) : Nat :=
  Erdos30FaceFieldExactObservables.prefixCount (n58WitnessOfIndex i) (n58PacketPrefixLocation i)

def n58PrefixWinnerIndicesByRank : List Nat :=
  n58IndexList.filter (fun i => n58PrefixRank i == 0)

def n58MassWinnerIndicesByRank : List Nat :=
  n58IndexList.filter (fun i => n58MassRank i == 0)

def n58JointWinnerIndicesByRank : List Nat :=
  n58IndexList.filter (fun i => n58JointRank i == 0)

def n58DominatesPrefixMass (i j : Nat) : Bool :=
  (n58PrefixRank i <= n58PrefixRank j) &&
  (n58MassRank i <= n58MassRank j) &&
  ((n58PrefixRank i < n58PrefixRank j) ||
   (n58MassRank i < n58MassRank j))

def n58ParetoIndicesByRanks : List Nat :=
  n58IndexList.filter (fun i => !(n58IndexList.any (fun j => n58DominatesPrefixMass j i)))

theorem n58_prefix_winner_match_packet :
    n58PrefixWinnerIndicesByRank = [2, 8, 9] := by
  native_decide

theorem n58_mass_winner_match_packet :
    n58MassWinnerIndicesByRank = [7] := by
  native_decide

theorem n58_joint_winner_match_packet :
    n58JointWinnerIndicesByRank = [7] := by
  native_decide

theorem n58_pareto_minimal_match_packet :
    n58ParetoIndicesByRanks = [7, 9] := by
  native_decide

theorem n58_exact_mass_twice_table :
    n58IndexList.map n58ExactMassTwice = [140, 36, 36, 80, 120, 16, 100, 4, 80, 24] := by
  native_decide

theorem n58_exact_mass_winner_match_packet :
    n58MassWinnerIndicesByExactMass = [7] := by
  native_decide

theorem n58_exact_mass_winner_matches_rank_winner :
    n58MassWinnerIndicesByExactMass = n58MassWinnerIndicesByRank := by
  native_decide

theorem n58_packet_prefix_location_table :
    n58IndexList.map n58PacketPrefixLocation = [55, 55, 58, 12, 56, 56, 57, 57, 58, 58] := by
  native_decide

theorem n58_prefix_count_at_packet_location_table :
    n58IndexList.map n58PrefixCountAtPacketLocation = [10, 10, 10, 4, 10, 10, 10, 10, 10, 10] := by
  native_decide

lemma sqrt58_ge_prefix_lower : (15/2 : ℝ) ≤ Real.sqrt (58 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 15/2) (by norm_num : (0:ℝ) ≤ 58)]
  norm_num

lemma sqrt58_lt_prefix_upper : Real.sqrt (58 : ℝ) < (8 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 8)]
  norm_num

lemma sqrt58_lt_prefix_tight_upper : Real.sqrt (58 : ℝ) < (23/3 : ℝ) := by
  rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
  norm_num

def n58W2PrefixCountTable : Nat -> Nat
  | 0 => 1
  | 1 => 1
  | 2 => 2
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
  | 15 => 3
  | 16 => 3
  | 17 => 3
  | 18 => 3
  | 19 => 3
  | 20 => 3
  | 21 => 4
  | 22 => 5
  | 23 => 5
  | 24 => 5
  | 25 => 5
  | 26 => 5
  | 27 => 5
  | 28 => 5
  | 29 => 5
  | 30 => 5
  | 31 => 5
  | 32 => 6
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
  | 46 => 7
  | 47 => 7
  | 48 => 7
  | 49 => 7
  | 50 => 8
  | 51 => 8
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 9
  | 56 => 9
  | 57 => 9
  | 58 => 10
  | _ => 10
theorem n58W2_prefix_count_table_matches_witness :
    (List.range 59).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n58W2 t) =
      [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10] := by
  native_decide

theorem n58W2_prefix_count_model_table :
    (List.range 59).map n58W2PrefixCountTable =
      [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 10] := by
  native_decide

lemma n58W2_prefix_segment_0_1_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_2_14_bound :
    ∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_2_14 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_15_20_bound :
    ∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_15_20 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_21_21_bound :
    ∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_21_21 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_22_31_bound :
    ∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_22_31 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_32_45_bound :
    ∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_32_45 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_46_49_bound :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_50_54_bound :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_55_57_bound :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_55_57 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W2_prefix_segment_58_58_bound :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W2_prefix_count_segment_58_58 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W2 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W2_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) := by
  exact ⟨n58W2_prefix_segment_0_1_bound, ⟨n58W2_prefix_segment_2_14_bound, ⟨n58W2_prefix_segment_15_20_bound, ⟨n58W2_prefix_segment_21_21_bound, ⟨n58W2_prefix_segment_22_31_bound, ⟨n58W2_prefix_segment_32_45_bound, ⟨n58W2_prefix_segment_46_49_bound, ⟨n58W2_prefix_segment_50_54_bound, ⟨n58W2_prefix_segment_55_57_bound, n58W2_prefix_segment_58_58_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n58W2_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W2 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W2.card = 10 := by
    native_decide
  have hcard : (n58W2.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]

theorem n58W2_prefix_residual_segment_0_1_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_0_1_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_2_14_zero :
    ∀ t : Nat, 2 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_2_14_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_2_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_15_20_zero :
    ∀ t : Nat, 15 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_15_20_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_15_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_21_21_zero :
    ∀ t : Nat, 21 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_21_21_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_21_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_22_31_zero :
    ∀ t : Nat, 22 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_22_31_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_22_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_32_45_zero :
    ∀ t : Nat, 32 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_32_45_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_32_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_46_49_zero :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_46_49_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_50_54_zero :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_50_54_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_55_57_zero :
    ∀ t : Nat, 55 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_55_57_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_55_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_prefix_residual_segment_58_58_zero :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W2_prefix_segment_58_58_bound t hlow hhigh
  have hcount := n58W2_prefix_count_segment_58_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W2_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W2 := by
  intro t ht
  interval_cases t
  · exact n58W2_prefix_residual_segment_0_1_zero 0 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_0_1_zero 1 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 2 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 3 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 4 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 5 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 6 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 7 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 8 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 9 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 10 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 11 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 12 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 13 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_2_14_zero 14 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 15 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 16 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 17 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 18 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 19 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_15_20_zero 20 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_21_21_zero 21 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 22 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 23 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 24 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 25 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 26 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 27 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 28 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 29 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 30 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_22_31_zero 31 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 32 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 33 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 34 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 35 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 36 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 37 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 38 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 39 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 40 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 41 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 42 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 43 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 44 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_32_45_zero 45 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_46_49_zero 46 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_46_49_zero 47 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_46_49_zero 48 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_46_49_zero 49 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_50_54_zero 50 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_50_54_zero 51 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_50_54_zero 52 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_50_54_zero 53 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_50_54_zero 54 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_55_57_zero 55 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_55_57_zero 56 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_55_57_zero 57 (by norm_num) (by norm_num)
  · exact n58W2_prefix_residual_segment_58_58_zero 58 (by norm_num) (by norm_num)

def n58W8PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 0
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
  | 26 => 5
  | 27 => 5
  | 28 => 5
  | 29 => 6
  | 30 => 6
  | 31 => 6
  | 32 => 6
  | 33 => 6
  | 34 => 6
  | 35 => 6
  | 36 => 6
  | 37 => 7
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
  | 56 => 9
  | 57 => 9
  | 58 => 10
  | _ => 10
theorem n58W8_prefix_count_table_matches_witness :
    (List.range 59).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n58W8 t) =
      [0, 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

theorem n58W8_prefix_count_model_table :
    (List.range 59).map n58W8PrefixCountTable =
      [0, 0, 0, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 10] := by
  native_decide

lemma n58W8_prefix_segment_0_2_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_3_3_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_3_3 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_4_8_bound :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (8 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_4_8 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_9_12_bound :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (9 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (12 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_9_12 :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_13_25_bound :
    ∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (13 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_13_25 :
    ∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_26_28_bound :
    ∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_26_28 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_29_36_bound :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_29_36 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_37_43_bound :
    ∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (43 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_37_43 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_44_55_bound :
    ∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (44 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_44_55 :
    ∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_56_57_bound :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W8_prefix_segment_58_58_bound :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W8_prefix_count_segment_58_58 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W8 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W8_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) := by
  exact ⟨n58W8_prefix_segment_0_2_bound, ⟨n58W8_prefix_segment_3_3_bound, ⟨n58W8_prefix_segment_4_8_bound, ⟨n58W8_prefix_segment_9_12_bound, ⟨n58W8_prefix_segment_13_25_bound, ⟨n58W8_prefix_segment_26_28_bound, ⟨n58W8_prefix_segment_29_36_bound, ⟨n58W8_prefix_segment_37_43_bound, ⟨n58W8_prefix_segment_44_55_bound, ⟨n58W8_prefix_segment_56_57_bound, n58W8_prefix_segment_58_58_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n58W8_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W8 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W8.card = 10 := by
    native_decide
  have hcard : (n58W8.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]

theorem n58W8_prefix_residual_segment_0_2_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_0_2_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_3_3_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_3_3_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_3_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_4_8_zero :
    ∀ t : Nat, 4 ≤ t -> t ≤ 8 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_4_8_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_4_8 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_9_12_zero :
    ∀ t : Nat, 9 ≤ t -> t ≤ 12 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_9_12_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_9_12 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_13_25_zero :
    ∀ t : Nat, 13 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_13_25_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_13_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_26_28_zero :
    ∀ t : Nat, 26 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_26_28_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_26_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_29_36_zero :
    ∀ t : Nat, 29 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_29_36_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_29_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_37_43_zero :
    ∀ t : Nat, 37 ≤ t -> t ≤ 43 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_37_43_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_37_43 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_44_55_zero :
    ∀ t : Nat, 44 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_44_55_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_44_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_56_57_zero :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_56_57_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_prefix_residual_segment_58_58_zero :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W8_prefix_segment_58_58_bound t hlow hhigh
  have hcount := n58W8_prefix_count_segment_58_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W8_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W8 := by
  intro t ht
  interval_cases t
  · exact n58W8_prefix_residual_segment_0_2_zero 0 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_0_2_zero 1 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_0_2_zero 2 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_3_3_zero 3 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_4_8_zero 4 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_4_8_zero 5 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_4_8_zero 6 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_4_8_zero 7 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_4_8_zero 8 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_9_12_zero 9 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_9_12_zero 10 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_9_12_zero 11 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_9_12_zero 12 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 13 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 14 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 15 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 16 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 17 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 18 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 19 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 20 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 21 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 22 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 23 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 24 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_13_25_zero 25 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_26_28_zero 26 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_26_28_zero 27 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_26_28_zero 28 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 29 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 30 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 31 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 32 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 33 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 34 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 35 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_29_36_zero 36 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 37 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 38 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 39 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 40 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 41 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 42 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_37_43_zero 43 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 44 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 45 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 46 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 47 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 48 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 49 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 50 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 51 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 52 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 53 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 54 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_44_55_zero 55 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_56_57_zero 56 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_56_57_zero 57 (by norm_num) (by norm_num)
  · exact n58W8_prefix_residual_segment_58_58_zero 58 (by norm_num) (by norm_num)

def n58W9PrefixCountTable : Nat -> Nat
  | 0 => 0
  | 1 => 0
  | 2 => 0
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
  | 14 => 2
  | 15 => 2
  | 16 => 2
  | 17 => 3
  | 18 => 3
  | 19 => 3
  | 20 => 3
  | 21 => 3
  | 22 => 3
  | 23 => 3
  | 24 => 4
  | 25 => 4
  | 26 => 4
  | 27 => 4
  | 28 => 4
  | 29 => 4
  | 30 => 4
  | 31 => 4
  | 32 => 5
  | 33 => 5
  | 34 => 5
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
  | 47 => 6
  | 48 => 7
  | 49 => 7
  | 50 => 7
  | 51 => 7
  | 52 => 8
  | 53 => 8
  | 54 => 8
  | 55 => 8
  | 56 => 8
  | 57 => 9
  | 58 => 10
  | _ => 10
theorem n58W9_prefix_count_table_matches_witness :
    (List.range 59).map (fun t => Erdos30FaceFieldExactObservables.prefixCount n58W9 t) =
      [0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

theorem n58W9_prefix_count_model_table :
    (List.range 59).map n58W9PrefixCountTable =
      [0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 10] := by
  native_decide

lemma n58W9_prefix_segment_0_2_bound :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_3_4_bound :
    ∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (4 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_3_4 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_5_16_bound :
    ∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (5 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (16 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_5_16 :
    ∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_17_23_bound :
    ∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (17 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_17_23 :
    ∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_24_31_bound :
    ∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_24_31 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_32_34_bound :
    ∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_32_34 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_35_47_bound :
    ∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (47 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_35_47 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_48_51_bound :
    ∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (48 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (51 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_48_51 :
    ∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_52_56_bound :
    ∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (52 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_52_56 :
    ∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_57_57_bound :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_57_57 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

lemma n58W9_prefix_segment_58_58_bound :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  rw [abs_le]
  constructor <;> nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper, hlowR, hhighR]

theorem n58W9_prefix_count_segment_58_58 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W9 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W9_full_prefix_segment_bound_certificate :
    (∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) ∧
    (∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        10 * Real.sqrt (58 : ℝ) - 58) := by
  exact ⟨n58W9_prefix_segment_0_2_bound, ⟨n58W9_prefix_segment_3_4_bound, ⟨n58W9_prefix_segment_5_16_bound, ⟨n58W9_prefix_segment_17_23_bound, ⟨n58W9_prefix_segment_24_31_bound, ⟨n58W9_prefix_segment_32_34_bound, ⟨n58W9_prefix_segment_35_47_bound, ⟨n58W9_prefix_segment_48_51_bound, ⟨n58W9_prefix_segment_52_56_bound, ⟨n58W9_prefix_segment_57_57_bound, n58W9_prefix_segment_58_58_bound⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem n58W9_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W9 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W9.card = 10 := by
    native_decide
  have hcard : (n58W9.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]

theorem n58W9_prefix_residual_segment_0_2_zero :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_0_2_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_3_4_zero :
    ∀ t : Nat, 3 ≤ t -> t ≤ 4 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_3_4_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_3_4 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_5_16_zero :
    ∀ t : Nat, 5 ≤ t -> t ≤ 16 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_5_16_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_5_16 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_17_23_zero :
    ∀ t : Nat, 17 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_17_23_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_17_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_24_31_zero :
    ∀ t : Nat, 24 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_24_31_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_24_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_32_34_zero :
    ∀ t : Nat, 32 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_32_34_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_32_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_35_47_zero :
    ∀ t : Nat, 35 ≤ t -> t ≤ 47 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_35_47_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_35_47 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_48_51_zero :
    ∀ t : Nat, 48 ≤ t -> t ≤ 51 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_48_51_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_48_51 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_52_56_zero :
    ∀ t : Nat, 52 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_52_56_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_52_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_57_57_zero :
    ∀ t : Nat, 57 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_57_57_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_57_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_prefix_residual_segment_58_58_zero :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 t = 0 := by
  intro t hlow hhigh
  have hdev := n58W9_prefix_segment_58_58_bound t hlow hhigh
  have hcount := n58W9_prefix_count_segment_58_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  norm_num at hdev ⊢
  exact hdev

theorem n58W9_full_prefix_residual_zero :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W9 := by
  intro t ht
  interval_cases t
  · exact n58W9_prefix_residual_segment_0_2_zero 0 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_0_2_zero 1 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_0_2_zero 2 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_3_4_zero 3 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_3_4_zero 4 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 5 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 6 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 7 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 8 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 9 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 10 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 11 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 12 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 13 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 14 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 15 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_5_16_zero 16 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 17 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 18 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 19 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 20 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 21 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 22 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_17_23_zero 23 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 24 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 25 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 26 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 27 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 28 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 29 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 30 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_24_31_zero 31 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_32_34_zero 32 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_32_34_zero 33 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_32_34_zero 34 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 35 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 36 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 37 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 38 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 39 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 40 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 41 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 42 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 43 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 44 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 45 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 46 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_35_47_zero 47 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_48_51_zero 48 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_48_51_zero 49 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_48_51_zero 50 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_48_51_zero 51 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_52_56_zero 52 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_52_56_zero 53 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_52_56_zero 54 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_52_56_zero 55 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_52_56_zero 56 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_57_57_zero 57 (by norm_num) (by norm_num)
  · exact n58W9_prefix_residual_segment_58_58_zero 58 (by norm_num) (by norm_num)


lemma sqrt58_ge_58_div_10 : (58/10 : ℝ) ≤ Real.sqrt (58 : ℝ) := by
  rw [Real.le_sqrt (by norm_num : (0:ℝ) ≤ 58/10) (by norm_num : (0:ℝ) ≤ 58)]
  norm_num


theorem n58W0_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W0 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W0.card = 10 := by
    native_decide
  have hcard : (n58W0.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W1_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W1 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W1.card = 10 := by
    native_decide
  have hcard : (n58W1.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W3_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W3 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W3.card = 10 := by
    native_decide
  have hcard : (n58W3.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W4_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W4 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W4.card = 10 := by
    native_decide
  have hcard : (n58W4.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W5_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W5 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W5.card = 10 := by
    native_decide
  have hcard : (n58W5.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W6_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W6 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W6.card = 10 := by
    native_decide
  have hcard : (n58W6.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58W7_prefix_drift_eq :
    Erdos30FaceFieldExactObservables.prefixDrift 58 n58W7 =
      10 * Real.sqrt (58 : ℝ) - 58 := by
  unfold Erdos30FaceFieldExactObservables.prefixDrift
  have hcardNat : n58W7.card = 10 := by
    native_decide
  have hcard : (n58W7.card : ℝ) = (10 : ℝ) := by
    exact_mod_cast hcardNat
  rw [hcard]
  have hnonneg : 0 ≤ (10 : ℝ) - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  norm_num at ⊢
  rw [abs_of_nonneg hnonneg]
  have hge1 : (1 : ℝ) ≤ 10 - Real.sqrt (58 : ℝ) := by
    nlinarith [sqrt58_lt_prefix_upper]
  rw [max_eq_left hge1]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 58)]
  nlinarith [hsq]


theorem n58_exact_prefix_probe_2_zero :
    n58ExactPrefixProbe 2 = 0 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n58_exact_prefix_probe_8_zero :
    n58ExactPrefixProbe 8 = 0 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n58_exact_prefix_probe_9_zero :
    n58ExactPrefixProbe 9 = 0 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  linarith


theorem n58_exact_prefix_probe_0_positive :
    0 < n58ExactPrefixProbe 0 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_1_positive :
    0 < n58ExactPrefixProbe 1 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_3_positive :
    0 < n58ExactPrefixProbe 3 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (12 : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (58 : ℝ) < (23/3 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
    norm_num
  have hgap :
      0 < -((12 : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)) -
        (10 * Real.sqrt (58 : ℝ) - 58) := by
    nlinarith [hsqrtUpper]
  linarith [hgap]


theorem n58_exact_prefix_probe_4_positive :
    0 < n58ExactPrefixProbe 4 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_5_positive :
    0 < n58ExactPrefixProbe 5 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_6_positive :
    0 < n58ExactPrefixProbe 6 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_7_positive :
    0 < n58ExactPrefixProbe 7 := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_prefix_probe_positive_on_nonwinners :
    0 < n58ExactPrefixProbe 0 ∧
    0 < n58ExactPrefixProbe 1 ∧
    0 < n58ExactPrefixProbe 3 ∧
    0 < n58ExactPrefixProbe 4 ∧
    0 < n58ExactPrefixProbe 5 ∧
    0 < n58ExactPrefixProbe 6 ∧
    0 < n58ExactPrefixProbe 7 := by
  exact ⟨n58_exact_prefix_probe_0_positive, ⟨n58_exact_prefix_probe_1_positive, ⟨n58_exact_prefix_probe_3_positive, ⟨n58_exact_prefix_probe_4_positive, ⟨n58_exact_prefix_probe_5_positive, ⟨n58_exact_prefix_probe_6_positive, n58_exact_prefix_probe_7_positive⟩⟩⟩⟩⟩⟩


theorem n58_exact_prefix_probe_0_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 0 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 (n58PacketPrefixLocation 0) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W0 (n58PacketPrefixLocation 0) =
        n58PacketPrefixCount 0 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_1_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 1 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 (n58PacketPrefixLocation 1) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W1 (n58PacketPrefixLocation 1) =
        n58PacketPrefixCount 1 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_2_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 2 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W2 (n58PacketPrefixLocation 2) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W2_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W2 (n58PacketPrefixLocation 2) =
        n58PacketPrefixCount 2 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_3_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 3 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 (n58PacketPrefixLocation 3) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W3 (n58PacketPrefixLocation 3) =
        n58PacketPrefixCount 3 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_4_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 4 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 (n58PacketPrefixLocation 4) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W4 (n58PacketPrefixLocation 4) =
        n58PacketPrefixCount 4 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_5_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 5 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 (n58PacketPrefixLocation 5) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W5 (n58PacketPrefixLocation 5) =
        n58PacketPrefixCount 5 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_6_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 6 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 (n58PacketPrefixLocation 6) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W6 (n58PacketPrefixLocation 6) =
        n58PacketPrefixCount 6 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_7_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 7 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 (n58PacketPrefixLocation 7) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W7 (n58PacketPrefixLocation 7) =
        n58PacketPrefixCount 7 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_8_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 8 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W8 (n58PacketPrefixLocation 8) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W8_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W8 (n58PacketPrefixLocation 8) =
        n58PacketPrefixCount 8 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_exact_prefix_probe_9_matches_prefix_residual_at_packet_location :
    n58ExactPrefixProbe 9 =
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W9 (n58PacketPrefixLocation 9) := by
  unfold n58ExactPrefixProbe
  unfold Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W9_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  have hcount :
      Erdos30FaceFieldExactObservables.prefixCount n58W9 (n58PacketPrefixLocation 9) =
        n58PacketPrefixCount 9 := by
    native_decide
  rw [hcount]
  rfl


theorem n58_prefix_residual_at_packet_location_0_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 (n58PacketPrefixLocation 0) := by
  rw [← n58_exact_prefix_probe_0_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_0_positive


theorem n58_prefix_residual_at_packet_location_1_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 (n58PacketPrefixLocation 1) := by
  rw [← n58_exact_prefix_probe_1_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_1_positive


theorem n58_prefix_residual_at_packet_location_3_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 (n58PacketPrefixLocation 3) := by
  rw [← n58_exact_prefix_probe_3_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_3_positive


theorem n58_prefix_residual_at_packet_location_4_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 (n58PacketPrefixLocation 4) := by
  rw [← n58_exact_prefix_probe_4_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_4_positive


theorem n58_prefix_residual_at_packet_location_5_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 (n58PacketPrefixLocation 5) := by
  rw [← n58_exact_prefix_probe_5_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_5_positive


theorem n58_prefix_residual_at_packet_location_6_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 (n58PacketPrefixLocation 6) := by
  rw [← n58_exact_prefix_probe_6_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_6_positive


theorem n58_prefix_residual_at_packet_location_7_positive :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 (n58PacketPrefixLocation 7) := by
  rw [← n58_exact_prefix_probe_7_matches_prefix_residual_at_packet_location]
  exact n58_exact_prefix_probe_7_positive


theorem n58_full_prefix_residual_zero_on_prefix_winners :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W2 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W8 ∧
        Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W9 := by
  exact ⟨n58W2_full_prefix_residual_zero, ⟨n58W8_full_prefix_residual_zero, n58W9_full_prefix_residual_zero⟩⟩


theorem n58_actual_prefix_residual_positive_on_probe_nonwinners :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 (n58PacketPrefixLocation 0) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 (n58PacketPrefixLocation 1) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 (n58PacketPrefixLocation 3) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 (n58PacketPrefixLocation 4) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 (n58PacketPrefixLocation 5) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 (n58PacketPrefixLocation 6) ∧
        0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 (n58PacketPrefixLocation 7) := by
  exact ⟨n58_prefix_residual_at_packet_location_0_positive, ⟨n58_prefix_residual_at_packet_location_1_positive, ⟨n58_prefix_residual_at_packet_location_3_positive, ⟨n58_prefix_residual_at_packet_location_4_positive, ⟨n58_prefix_residual_at_packet_location_5_positive, ⟨n58_prefix_residual_at_packet_location_6_positive, n58_prefix_residual_at_packet_location_7_positive⟩⟩⟩⟩⟩⟩


theorem n58_full_prefix_residual_max_2_zero :
    n58FullPrefixResidualMax 2 = 0 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n58W2_full_prefix_residual_zero


theorem n58_full_prefix_residual_max_8_zero :
    n58FullPrefixResidualMax 8 = 0 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n58W8_full_prefix_residual_zero


theorem n58_full_prefix_residual_max_9_zero :
    n58FullPrefixResidualMax 9 = 0 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero
    n58W9_full_prefix_residual_zero


theorem n58_full_prefix_residual_max_0_positive :
    0 < n58FullPrefixResidualMax 0 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 0 ≤ 58)
    n58_prefix_residual_at_packet_location_0_positive


theorem n58_full_prefix_residual_max_1_positive :
    0 < n58FullPrefixResidualMax 1 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 1 ≤ 58)
    n58_prefix_residual_at_packet_location_1_positive


theorem n58_full_prefix_residual_max_3_positive :
    0 < n58FullPrefixResidualMax 3 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 3 ≤ 58)
    n58_prefix_residual_at_packet_location_3_positive


theorem n58_full_prefix_residual_max_4_positive :
    0 < n58FullPrefixResidualMax 4 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 4 ≤ 58)
    n58_prefix_residual_at_packet_location_4_positive


theorem n58_full_prefix_residual_max_5_positive :
    0 < n58FullPrefixResidualMax 5 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 5 ≤ 58)
    n58_prefix_residual_at_packet_location_5_positive


theorem n58_full_prefix_residual_max_6_positive :
    0 < n58FullPrefixResidualMax 6 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 6 ≤ 58)
    n58_prefix_residual_at_packet_location_6_positive


theorem n58_full_prefix_residual_max_7_positive :
    0 < n58FullPrefixResidualMax 7 := by
  unfold n58FullPrefixResidualMax
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_pos_of_pos_at
    (by native_decide : n58PacketPrefixLocation 7 ≤ 58)
    n58_prefix_residual_at_packet_location_7_positive


theorem n58_full_prefix_residual_max_zero_on_prefix_winners :
    n58FullPrefixResidualMax 2 = 0 ∧
        n58FullPrefixResidualMax 8 = 0 ∧
        n58FullPrefixResidualMax 9 = 0 := by
  exact ⟨n58_full_prefix_residual_max_2_zero, ⟨n58_full_prefix_residual_max_8_zero, n58_full_prefix_residual_max_9_zero⟩⟩


theorem n58_full_prefix_residual_max_positive_on_probe_nonwinners :
    0 < n58FullPrefixResidualMax 0 ∧
        0 < n58FullPrefixResidualMax 1 ∧
        0 < n58FullPrefixResidualMax 3 ∧
        0 < n58FullPrefixResidualMax 4 ∧
        0 < n58FullPrefixResidualMax 5 ∧
        0 < n58FullPrefixResidualMax 6 ∧
        0 < n58FullPrefixResidualMax 7 := by
  exact ⟨n58_full_prefix_residual_max_0_positive, ⟨n58_full_prefix_residual_max_1_positive, ⟨n58_full_prefix_residual_max_3_positive, ⟨n58_full_prefix_residual_max_4_positive, ⟨n58_full_prefix_residual_max_5_positive, ⟨n58_full_prefix_residual_max_6_positive, n58_full_prefix_residual_max_7_positive⟩⟩⟩⟩⟩⟩


theorem n58_exact_prefix_probe_0_value :
    n58ExactPrefixProbe 0 = (3 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_0_value :
    n58ExactMassTwice 0 = 140 := by
  native_decide


theorem n58_exact_prefix_probe_1_value :
    n58ExactPrefixProbe 1 = (3 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (55 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_1_value :
    n58ExactMassTwice 1 = 36 := by
  native_decide


theorem n58_exact_prefix_probe_2_value :
    n58ExactPrefixProbe 2 = (0 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_2_value :
    n58ExactMassTwice 2 = 36 := by
  native_decide


theorem n58_exact_prefix_probe_3_value :
    n58ExactPrefixProbe 3 = ((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (12 : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_prefix_lower]
  rw [abs_of_nonpos hnonpos]
  have hsqrtUpper : Real.sqrt (58 : ℝ) < (23/3 : ℝ) := by
    rw [Real.sqrt_lt' (by norm_num : (0:ℝ) < 23/3)]
    norm_num
  have hgap :
      0 < -((12 : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)) -
        (10 * Real.sqrt (58 : ℝ) - 58) := by
    nlinarith [hsqrtUpper]
  rw [max_eq_right (le_of_lt hgap)]
  ring


theorem n58_exact_mass_twice_3_value :
    n58ExactMassTwice 3 = 80 := by
  native_decide


theorem n58_exact_prefix_probe_4_value :
    n58ExactPrefixProbe 4 = (2 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_4_value :
    n58ExactMassTwice 4 = 120 := by
  native_decide


theorem n58_exact_prefix_probe_5_value :
    n58ExactPrefixProbe 5 = (2 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (56 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_5_value :
    n58ExactMassTwice 5 = 16 := by
  native_decide


theorem n58_exact_prefix_probe_6_value :
    n58ExactPrefixProbe 6 = (1 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_6_value :
    n58ExactMassTwice 6 = 100 := by
  native_decide


theorem n58_exact_prefix_probe_7_value :
    n58ExactPrefixProbe 7 = (1 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (57 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_7_value :
    n58ExactMassTwice 7 = 4 := by
  native_decide


theorem n58_exact_prefix_probe_8_value :
    n58ExactPrefixProbe 8 = (0 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_8_value :
    n58ExactMassTwice 8 = 80 := by
  native_decide


theorem n58_exact_prefix_probe_9_value :
    n58ExactPrefixProbe 9 = (0 : ℝ) := by
  simp [n58ExactPrefixProbe, Erdos30FaceFieldExactObservables.prefixResidualProbeCard10,
    n58PacketPrefixLocation, n58PacketPrefixCount]
  have hnonpos : (58 : ℝ) - 10 * Real.sqrt (58 : ℝ) ≤ 0 := by
    nlinarith [sqrt58_ge_58_div_10]
  rw [abs_of_nonpos hnonpos]
  norm_num


theorem n58_exact_mass_twice_9_value :
    n58ExactMassTwice 9 = 24 := by
  native_decide


theorem n58W0_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_1_5 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_1_5_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 5 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (5 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_1_5 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_6_9 :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_6_9_le_exact_probe :
    ∀ t : Nat, 6 ≤ t -> t ≤ 9 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (6 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (9 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_6_9 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_10_22 :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_10_22_le_exact_probe :
    ∀ t : Nat, 10 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (10 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_10_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_23_25 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_23_25_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_23_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_26_33 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_26_33_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_26_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_34_40 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_34_40_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 40 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (40 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_34_40 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_41_52 :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_41_52_le_exact_probe :
    ∀ t : Nat, 41 ≤ t -> t ≤ 52 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (41 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (52 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_41_52 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_53_54 :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_53_54_le_exact_probe :
    ∀ t : Nat, 53 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (53 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_53_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_prefix_count_segment_55_58 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W0 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W0_prefix_residual_segment_55_58_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W0_prefix_count_segment_55_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W0_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_0_positive
  · rw [n58_exact_prefix_probe_0_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W0_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 t ≤
        n58ExactPrefixProbe 0 := by
  intro t ht
  interval_cases t
  · exact n58W0_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_1_5_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_1_5_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_1_5_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_1_5_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_1_5_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_6_9_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_6_9_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_6_9_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_6_9_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_10_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_23_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_23_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_23_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_26_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_34_40_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_41_52_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_53_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_53_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_55_58_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_55_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_55_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W0_prefix_residual_segment_55_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W1_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_2_13 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_2_13_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 13 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (13 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_2_13 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_14_20 :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_14_20_le_exact_probe :
    ∀ t : Nat, 14 ≤ t -> t ≤ 20 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (14 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (20 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_14_20 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_21_28 :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_21_28_le_exact_probe :
    ∀ t : Nat, 21 ≤ t -> t ≤ 28 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (21 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (28 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_21_28 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_29_31 :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_29_31_le_exact_probe :
    ∀ t : Nat, 29 ≤ t -> t ≤ 31 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (29 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (31 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_29_31 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_32_44 :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_32_44_le_exact_probe :
    ∀ t : Nat, 32 ≤ t -> t ≤ 44 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (32 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (44 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_32_44 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_45_48 :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_45_48_le_exact_probe :
    ∀ t : Nat, 45 ≤ t -> t ≤ 48 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (45 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (48 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_45_48 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_49_53 :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_49_53_le_exact_probe :
    ∀ t : Nat, 49 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (49 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_49_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_54_54 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_54_54_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_54_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_prefix_count_segment_55_58 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W1 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W1_prefix_residual_segment_55_58_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (3 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W1_prefix_count_segment_55_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W1_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_1_positive
  · rw [n58_exact_prefix_probe_1_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W1_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 t ≤
        n58ExactPrefixProbe 1 := by
  intro t ht
  interval_cases t
  · exact n58W1_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_2_13_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_14_20_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_21_28_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_29_31_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_29_31_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_29_31_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_32_44_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_45_48_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_45_48_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_45_48_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_45_48_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_49_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_49_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_49_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_49_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_49_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_54_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_55_58_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_55_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_55_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W1_prefix_residual_segment_55_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W3_prefix_count_segment_0_2 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_0_2_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_0_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_3_7 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_3_7_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_3_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_8_11 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_8_11_le_exact_probe :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (11 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_8_11 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_12_25 :
    ∀ t : Nat, 12 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_12_25_le_exact_probe :
    ∀ t : Nat, 12 ≤ t -> t ≤ 25 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (12 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (25 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_12_25 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_26_35 :
    ∀ t : Nat, 26 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_26_35_le_exact_probe :
    ∀ t : Nat, 26 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (26 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_26_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_36_36 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_36_36_le_exact_probe :
    ∀ t : Nat, 36 ≤ t -> t ≤ 36 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (36 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_36_36 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_37_42 :
    ∀ t : Nat, 37 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_37_42_le_exact_probe :
    ∀ t : Nat, 37 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (37 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_37_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_43_55 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_43_55_le_exact_probe :
    ∀ t : Nat, 43 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_43_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_56_57 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_56_57_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 57 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (57 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_56_57 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_prefix_count_segment_58_58 :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W3 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W3_prefix_residual_segment_58_58_le_exact_probe :
    ∀ t : Nat, 58 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t hlow hhigh
  have hlowR : (58 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (((46 : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)) : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W3_prefix_count_segment_58_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W3_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_3_positive
  · rw [n58_exact_prefix_probe_3_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W3_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 t ≤
        n58ExactPrefixProbe 3 := by
  intro t ht
  interval_cases t
  · exact n58W3_prefix_residual_segment_0_2_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_0_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_0_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_3_7_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_3_7_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_3_7_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_3_7_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_3_7_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_8_11_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_8_11_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_8_11_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_8_11_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_12_25_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_26_35_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_36_36_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_37_42_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_43_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_56_57_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_56_57_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W3_prefix_residual_segment_58_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W4_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_1_1 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_1_1_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_1_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_2_6 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_2_6_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 6 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (6 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_2_6 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_7_10 :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_7_10_le_exact_probe :
    ∀ t : Nat, 7 ≤ t -> t ≤ 10 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (7 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (10 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_7_10 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_11_23 :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_11_23_le_exact_probe :
    ∀ t : Nat, 11 ≤ t -> t ≤ 23 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (11 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (23 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_11_23 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_24_26 :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_24_26_le_exact_probe :
    ∀ t : Nat, 24 ≤ t -> t ≤ 26 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (24 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (26 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_24_26 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_27_34 :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_27_34_le_exact_probe :
    ∀ t : Nat, 27 ≤ t -> t ≤ 34 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (27 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (34 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_27_34 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_35_41 :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_35_41_le_exact_probe :
    ∀ t : Nat, 35 ≤ t -> t ≤ 41 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (35 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (41 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_35_41 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_42_53 :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_42_53_le_exact_probe :
    ∀ t : Nat, 42 ≤ t -> t ≤ 53 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (42 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (53 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_42_53 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_54_55 :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_54_55_le_exact_probe :
    ∀ t : Nat, 54 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (54 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_54_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_prefix_count_segment_56_58 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W4 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W4_prefix_residual_segment_56_58_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W4_prefix_count_segment_56_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W4_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_4_positive
  · rw [n58_exact_prefix_probe_4_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W4_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 t ≤
        n58ExactPrefixProbe 4 := by
  intro t ht
  interval_cases t
  · exact n58W4_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_1_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_2_6_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_2_6_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_2_6_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_2_6_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_2_6_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_7_10_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_7_10_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_7_10_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_7_10_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_11_23_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_24_26_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_24_26_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_24_26_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_27_34_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_35_41_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_42_53_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_54_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_54_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_56_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_56_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W4_prefix_residual_segment_56_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W5_prefix_count_segment_0_0 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_0_0_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 0 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (0 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_0_0 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_1_2 :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_1_2_le_exact_probe :
    ∀ t : Nat, 1 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (1 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_1_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_3_14 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_3_14_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 14 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (14 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_3_14 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_15_21 :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_15_21_le_exact_probe :
    ∀ t : Nat, 15 ≤ t -> t ≤ 21 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (15 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (21 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_15_21 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_22_29 :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_22_29_le_exact_probe :
    ∀ t : Nat, 22 ≤ t -> t ≤ 29 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (22 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (29 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_22_29 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_30_32 :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_30_32_le_exact_probe :
    ∀ t : Nat, 30 ≤ t -> t ≤ 32 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (30 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (32 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_30_32 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_33_45 :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_33_45_le_exact_probe :
    ∀ t : Nat, 33 ≤ t -> t ≤ 45 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (33 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (45 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_33_45 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_46_49 :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_46_49_le_exact_probe :
    ∀ t : Nat, 46 ≤ t -> t ≤ 49 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (46 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (49 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_46_49 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_50_54 :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_50_54_le_exact_probe :
    ∀ t : Nat, 50 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (50 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_50_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_55_55 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_55_55_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_55_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_prefix_count_segment_56_58 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W5 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W5_prefix_residual_segment_56_58_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (2 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W5_prefix_count_segment_56_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W5_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_5_positive
  · rw [n58_exact_prefix_probe_5_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W5_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 t ≤
        n58ExactPrefixProbe 5 := by
  intro t ht
  interval_cases t
  · exact n58W5_prefix_residual_segment_0_0_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_1_2_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_1_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_3_14_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_15_21_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_22_29_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_30_32_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_30_32_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_30_32_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_33_45_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_46_49_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_46_49_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_46_49_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_46_49_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_50_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_50_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_50_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_50_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_50_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_55_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_56_58_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_56_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W5_prefix_residual_segment_56_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W6_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_2_2 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_2_2_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 2 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (2 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_2_2 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_3_7 :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_3_7_le_exact_probe :
    ∀ t : Nat, 3 ≤ t -> t ≤ 7 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (3 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (7 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_3_7 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_8_11 :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_8_11_le_exact_probe :
    ∀ t : Nat, 8 ≤ t -> t ≤ 11 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (8 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (11 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_8_11 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_12_24 :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_12_24_le_exact_probe :
    ∀ t : Nat, 12 ≤ t -> t ≤ 24 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (12 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (24 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_12_24 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_25_27 :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_25_27_le_exact_probe :
    ∀ t : Nat, 25 ≤ t -> t ≤ 27 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (25 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (27 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_25_27 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_28_35 :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_28_35_le_exact_probe :
    ∀ t : Nat, 28 ≤ t -> t ≤ 35 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (28 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (35 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_28_35 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_36_42 :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_36_42_le_exact_probe :
    ∀ t : Nat, 36 ≤ t -> t ≤ 42 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (36 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (42 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_36_42 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_43_54 :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_43_54_le_exact_probe :
    ∀ t : Nat, 43 ≤ t -> t ≤ 54 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (43 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (54 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_43_54 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_55_56 :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_55_56_le_exact_probe :
    ∀ t : Nat, 55 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (55 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_55_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_prefix_count_segment_57_58 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W6 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W6_prefix_residual_segment_57_58_le_exact_probe :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W6_prefix_count_segment_57_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W6_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_6_positive
  · rw [n58_exact_prefix_probe_6_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W6_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 t ≤
        n58ExactPrefixProbe 6 := by
  intro t ht
  interval_cases t
  · exact n58W6_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_2_2_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_3_7_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_3_7_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_3_7_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_3_7_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_3_7_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_8_11_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_8_11_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_8_11_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_8_11_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_12_24_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_25_27_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_25_27_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_25_27_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_28_35_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_36_42_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_43_54_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_55_56_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_55_56_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_57_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W6_prefix_residual_segment_57_58_le_exact_probe 58 (by norm_num) (by norm_num)

theorem n58W7_prefix_count_segment_0_1 :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 0 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_0_1_le_exact_probe :
    ∀ t : Nat, 0 ≤ t -> t ≤ 1 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (0 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (1 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (0 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_0_1 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_2_3 :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 1 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_2_3_le_exact_probe :
    ∀ t : Nat, 2 ≤ t -> t ≤ 3 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (2 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (3 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (1 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_2_3 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_4_15 :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 2 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_4_15_le_exact_probe :
    ∀ t : Nat, 4 ≤ t -> t ≤ 15 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (4 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (15 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (2 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_4_15 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_16_22 :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 3 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_16_22_le_exact_probe :
    ∀ t : Nat, 16 ≤ t -> t ≤ 22 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (16 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (22 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (3 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_16_22 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_23_30 :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 4 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_23_30_le_exact_probe :
    ∀ t : Nat, 23 ≤ t -> t ≤ 30 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (23 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (30 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (4 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_23_30 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_31_33 :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 5 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_31_33_le_exact_probe :
    ∀ t : Nat, 31 ≤ t -> t ≤ 33 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (31 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (33 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (5 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_31_33 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_34_46 :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 6 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_34_46_le_exact_probe :
    ∀ t : Nat, 34 ≤ t -> t ≤ 46 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (34 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (46 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (6 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_34_46 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_47_50 :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 7 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_47_50_le_exact_probe :
    ∀ t : Nat, 47 ≤ t -> t ≤ 50 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (47 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (50 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (7 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_47_50 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_51_55 :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 8 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_51_55_le_exact_probe :
    ∀ t : Nat, 51 ≤ t -> t ≤ 55 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (51 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (55 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (8 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_51_55 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_56_56 :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 9 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_56_56_le_exact_probe :
    ∀ t : Nat, 56 ≤ t -> t ≤ 56 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (56 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (56 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (9 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_56_56 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_prefix_count_segment_57_58 :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixCount n58W7 t = 10 := by
  intro t hlow hhigh
  interval_cases t <;> native_decide

theorem n58W7_prefix_residual_segment_57_58_le_exact_probe :
    ∀ t : Nat, 57 ≤ t -> t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t hlow hhigh
  have hlowR : (57 : ℝ) ≤ (t : ℝ) := by exact_mod_cast hlow
  have hhighR : (t : ℝ) ≤ (58 : ℝ) := by exact_mod_cast hhigh
  have hdev :
      |(t : ℝ) - (10 : ℝ) * Real.sqrt (58 : ℝ)| ≤
        (10 * Real.sqrt (58 : ℝ) - 58) + (1 : ℝ) := by
    rw [abs_le]
    constructor
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
    · norm_num at ⊢
      nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper,
        sqrt58_lt_prefix_tight_upper, hlowR, hhighR]
  have hcount := n58W7_prefix_count_segment_57_58 t hlow hhigh
  unfold Erdos30FaceFieldExactObservables.prefixResidualAt
  rw [n58W7_prefix_drift_eq]
  unfold Erdos30FaceFieldExactObservables.prefixDeviationAt
  rw [hcount]
  apply max_le
  · exact le_of_lt n58_exact_prefix_probe_7_positive
  · rw [n58_exact_prefix_probe_7_value]
    norm_num at hdev ⊢
    nlinarith [hdev]

theorem n58W7_full_prefix_residual_le_exact_probe :
    ∀ t : Nat, t ≤ 58 ->
      Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 t ≤
        n58ExactPrefixProbe 7 := by
  intro t ht
  interval_cases t
  · exact n58W7_prefix_residual_segment_0_1_le_exact_probe 0 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_0_1_le_exact_probe 1 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_2_3_le_exact_probe 2 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_2_3_le_exact_probe 3 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 4 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 5 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 6 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 7 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 8 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 9 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 10 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 11 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 12 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 13 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 14 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_4_15_le_exact_probe 15 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 16 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 17 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 18 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 19 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 20 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 21 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_16_22_le_exact_probe 22 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 23 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 24 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 25 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 26 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 27 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 28 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 29 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_23_30_le_exact_probe 30 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_31_33_le_exact_probe 31 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_31_33_le_exact_probe 32 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_31_33_le_exact_probe 33 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 34 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 35 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 36 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 37 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 38 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 39 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 40 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 41 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 42 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 43 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 44 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 45 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_34_46_le_exact_probe 46 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_47_50_le_exact_probe 47 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_47_50_le_exact_probe 48 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_47_50_le_exact_probe 49 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_47_50_le_exact_probe 50 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_51_55_le_exact_probe 51 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_51_55_le_exact_probe 52 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_51_55_le_exact_probe 53 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_51_55_le_exact_probe 54 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_51_55_le_exact_probe 55 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_56_56_le_exact_probe 56 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_57_58_le_exact_probe 57 (by norm_num) (by norm_num)
  · exact n58W7_prefix_residual_segment_57_58_le_exact_probe 58 (by norm_num) (by norm_num)


theorem n58_full_prefix_residual_max_0_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 0 = n58ExactPrefixProbe 0 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 0 ≤ 58)
    n58W0_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_0_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_1_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 1 = n58ExactPrefixProbe 1 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 1 ≤ 58)
    n58W1_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_1_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_2_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 2 = n58ExactPrefixProbe 2 := by
  rw [n58_full_prefix_residual_max_2_zero, n58_exact_prefix_probe_2_zero]


theorem n58_full_prefix_residual_max_3_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 3 = n58ExactPrefixProbe 3 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 3 ≤ 58)
    n58W3_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_3_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_4_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 4 = n58ExactPrefixProbe 4 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 4 ≤ 58)
    n58W4_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_4_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_5_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 5 = n58ExactPrefixProbe 5 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 5 ≤ 58)
    n58W5_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_5_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_6_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 6 = n58ExactPrefixProbe 6 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 6 ≤ 58)
    n58W6_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_6_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_7_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 7 = n58ExactPrefixProbe 7 := by
  unfold n58FullPrefixResidualMax n58WitnessOfIndex
  exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at
    (by native_decide : n58PacketPrefixLocation 7 ≤ 58)
    n58W7_full_prefix_residual_le_exact_probe
    (n58_exact_prefix_probe_7_matches_prefix_residual_at_packet_location).symm


theorem n58_full_prefix_residual_max_8_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 8 = n58ExactPrefixProbe 8 := by
  rw [n58_full_prefix_residual_max_8_zero, n58_exact_prefix_probe_8_zero]


theorem n58_full_prefix_residual_max_9_matches_exact_prefix_probe :
    n58FullPrefixResidualMax 9 = n58ExactPrefixProbe 9 := by
  rw [n58_full_prefix_residual_max_9_zero, n58_exact_prefix_probe_9_zero]


theorem n58_full_prefix_residual_max_matches_exact_prefix_probe_on_face :
    ∀ i ∈ n58FaceIndices, n58FullPrefixResidualMax i = n58ExactPrefixProbe i := by
  intro i hi
  simp [n58FaceIndices] at hi
  interval_cases i
  · exact n58_full_prefix_residual_max_0_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_1_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_2_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_3_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_4_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_5_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_6_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_7_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_8_matches_exact_prefix_probe
  · exact n58_full_prefix_residual_max_9_matches_exact_prefix_probe


theorem n58_exact_joint_key_0_value :
    n58ExactJointKey 0 = (3 : ℝ) * Real.sqrt (58 : ℝ) + (70 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_0_value, n58_exact_mass_twice_0_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_0_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 0 = n58ExactJointKey 0 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_0_matches_exact_prefix_probe]


theorem n58_exact_joint_key_1_value :
    n58ExactJointKey 1 = (3 : ℝ) * Real.sqrt (58 : ℝ) + (18 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_1_value, n58_exact_mass_twice_1_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_1_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 1 = n58ExactJointKey 1 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_1_matches_exact_prefix_probe]


theorem n58_exact_joint_key_2_value :
    n58ExactJointKey 2 = (18 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_2_value, n58_exact_mass_twice_2_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_2_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 2 = n58ExactJointKey 2 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_2_matches_exact_prefix_probe]


theorem n58_exact_joint_key_3_value :
    n58ExactJointKey 3 = (46 : ℝ) * Real.sqrt (58 : ℝ) - (308 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_3_value, n58_exact_mass_twice_3_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_3_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 3 = n58ExactJointKey 3 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_3_matches_exact_prefix_probe]


theorem n58_exact_joint_key_4_value :
    n58ExactJointKey 4 = (2 : ℝ) * Real.sqrt (58 : ℝ) + (60 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_4_value, n58_exact_mass_twice_4_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_4_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 4 = n58ExactJointKey 4 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_4_matches_exact_prefix_probe]


theorem n58_exact_joint_key_5_value :
    n58ExactJointKey 5 = (2 : ℝ) * Real.sqrt (58 : ℝ) + (8 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_5_value, n58_exact_mass_twice_5_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_5_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 5 = n58ExactJointKey 5 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_5_matches_exact_prefix_probe]


theorem n58_exact_joint_key_6_value :
    n58ExactJointKey 6 = (1 : ℝ) * Real.sqrt (58 : ℝ) + (50 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_6_value, n58_exact_mass_twice_6_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_6_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 6 = n58ExactJointKey 6 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_6_matches_exact_prefix_probe]


theorem n58_exact_joint_key_7_value :
    n58ExactJointKey 7 = (1 : ℝ) * Real.sqrt (58 : ℝ) + (2 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_7_value, n58_exact_mass_twice_7_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 7 = n58ExactJointKey 7 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_7_matches_exact_prefix_probe]


theorem n58_exact_joint_key_8_value :
    n58ExactJointKey 8 = (40 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_8_value, n58_exact_mass_twice_8_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_8_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 8 = n58ExactJointKey 8 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_8_matches_exact_prefix_probe]


theorem n58_exact_joint_key_9_value :
    n58ExactJointKey 9 = (12 : ℝ) := by
  unfold n58ExactJointKey
  rw [n58_exact_prefix_probe_9_value, n58_exact_mass_twice_9_value]
  have hsq : (Real.sqrt (58 : ℝ)) ^ 2 = (58 : ℝ) := by
    rw [Real.sq_sqrt (by norm_num : (0:ℝ) ≤ 58)]
  nlinarith


theorem n58_scalar_full_prefix_joint_key_9_matches_exact_joint_key :
    n58ScalarFullPrefixJointKey 9 = n58ExactJointKey 9 := by
  unfold n58ScalarFullPrefixJointKey n58ExactJointKey
  rw [n58_full_prefix_residual_max_9_matches_exact_prefix_probe]


theorem n58_exact_joint_key_7_lt_0 :
    n58ExactJointKey 7 < n58ExactJointKey 0 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_0_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_0 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 0 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_0_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_0


theorem n58_exact_joint_key_7_lt_1 :
    n58ExactJointKey 7 < n58ExactJointKey 1 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_1_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_1 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 1 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_1_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_1


theorem n58_exact_joint_key_7_lt_2 :
    n58ExactJointKey 7 < n58ExactJointKey 2 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_2_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_2 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 2 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_2_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_2


theorem n58_exact_joint_key_7_lt_3 :
    n58ExactJointKey 7 < n58ExactJointKey 3 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_3_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_3 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 3 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_3_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_3


theorem n58_exact_joint_key_7_lt_4 :
    n58ExactJointKey 7 < n58ExactJointKey 4 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_4_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_4 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 4 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_4_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_4


theorem n58_exact_joint_key_7_lt_5 :
    n58ExactJointKey 7 < n58ExactJointKey 5 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_5_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_5 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 5 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_5_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_5


theorem n58_exact_joint_key_7_lt_6 :
    n58ExactJointKey 7 < n58ExactJointKey 6 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_6_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_6 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 6 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_6_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_6


theorem n58_exact_joint_key_7_lt_8 :
    n58ExactJointKey 7 < n58ExactJointKey 8 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_8_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_8 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 8 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_8_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_8


theorem n58_exact_joint_key_7_lt_9 :
    n58ExactJointKey 7 < n58ExactJointKey 9 := by
  rw [n58_exact_joint_key_7_value, n58_exact_joint_key_9_value]
  nlinarith [sqrt58_ge_prefix_lower, sqrt58_lt_prefix_upper]


theorem n58_scalar_full_prefix_joint_key_7_lt_9 :
    n58ScalarFullPrefixJointKey 7 <
      n58ScalarFullPrefixJointKey 9 := by
  rw [n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key,
    n58_scalar_full_prefix_joint_key_9_matches_exact_joint_key]
  exact n58_exact_joint_key_7_lt_9


theorem n58_exact_joint_winner_strict_certificate :
    n58ExactJointKey 7 < n58ExactJointKey 0 ∧
    n58ExactJointKey 7 < n58ExactJointKey 1 ∧
    n58ExactJointKey 7 < n58ExactJointKey 2 ∧
    n58ExactJointKey 7 < n58ExactJointKey 3 ∧
    n58ExactJointKey 7 < n58ExactJointKey 4 ∧
    n58ExactJointKey 7 < n58ExactJointKey 5 ∧
    n58ExactJointKey 7 < n58ExactJointKey 6 ∧
    n58ExactJointKey 7 < n58ExactJointKey 8 ∧
    n58ExactJointKey 7 < n58ExactJointKey 9 := by
  exact ⟨n58_exact_joint_key_7_lt_0, ⟨n58_exact_joint_key_7_lt_1, ⟨n58_exact_joint_key_7_lt_2, ⟨n58_exact_joint_key_7_lt_3, ⟨n58_exact_joint_key_7_lt_4, ⟨n58_exact_joint_key_7_lt_5, ⟨n58_exact_joint_key_7_lt_6, ⟨n58_exact_joint_key_7_lt_8, n58_exact_joint_key_7_lt_9⟩⟩⟩⟩⟩⟩⟩⟩


theorem n58_scalar_full_prefix_joint_winner_strict_certificate :
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 0 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 1 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 2 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 3 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 4 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 5 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 6 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 8 ∧
    n58ScalarFullPrefixJointKey 7 < n58ScalarFullPrefixJointKey 9 := by
  exact ⟨n58_scalar_full_prefix_joint_key_7_lt_0, ⟨n58_scalar_full_prefix_joint_key_7_lt_1, ⟨n58_scalar_full_prefix_joint_key_7_lt_2, ⟨n58_scalar_full_prefix_joint_key_7_lt_3, ⟨n58_scalar_full_prefix_joint_key_7_lt_4, ⟨n58_scalar_full_prefix_joint_key_7_lt_5, ⟨n58_scalar_full_prefix_joint_key_7_lt_6, ⟨n58_scalar_full_prefix_joint_key_7_lt_8, n58_scalar_full_prefix_joint_key_7_lt_9⟩⟩⟩⟩⟩⟩⟩⟩


theorem n58_prefix_min_2 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58PrefixRank 2 := by
  constructor
  · native_decide
  · intro y hy
    simp [n58FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n58_mass_min_7 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58MassRank 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n58FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n58_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58JointRank 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n58FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n58_exact_mass_min_7 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ExactMassTwice 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n58FaceIndices] at hy
    interval_cases y <;> native_decide


theorem n58_exact_prefix_probe_min_2 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ExactPrefixProbe 2 := by
  constructor
  · native_decide
  · intro y hy
    rw [n58_exact_prefix_probe_2_zero]
    unfold n58ExactPrefixProbe Erdos30FaceFieldExactObservables.prefixResidualProbeCard10
    exact le_max_left 0 _


theorem n58_full_prefix_residual_max_min_2 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58FullPrefixResidualMax 2 := by
  constructor
  · native_decide
  · intro y hy
    rw [n58_full_prefix_residual_max_2_zero]
    unfold n58FullPrefixResidualMax
    exact Erdos30FaceFieldExactObservables.fullPrefixResidualMax_nonneg 58 (n58WitnessOfIndex y)


theorem n58_exact_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ExactJointKey 7 := by
  constructor
  · native_decide
  · intro y hy
    simp [n58FaceIndices] at hy
    interval_cases y
    · exact le_of_lt n58_exact_joint_key_7_lt_0
    · exact le_of_lt n58_exact_joint_key_7_lt_1
    · exact le_of_lt n58_exact_joint_key_7_lt_2
    · exact le_of_lt n58_exact_joint_key_7_lt_3
    · exact le_of_lt n58_exact_joint_key_7_lt_4
    · exact le_of_lt n58_exact_joint_key_7_lt_5
    · exact le_of_lt n58_exact_joint_key_7_lt_6
    · exact le_rfl
    · exact le_of_lt n58_exact_joint_key_7_lt_8
    · exact le_of_lt n58_exact_joint_key_7_lt_9


theorem n58_scalar_full_prefix_joint_min_7 :
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ScalarFullPrefixJointKey 7 := by
  exact Erdos.Collider.isFieldMinOn_weightedJoint_of_prefix_eq_on
    (F := n58FaceIndices)
    (x := 7)
    (prefixScalar := n58FullPrefixResidualMax)
    (prefixProbe := n58ExactPrefixProbe)
    (mass := fun i => (n58ExactMassTwice i : ℝ) / 2)
    (weight := Real.sqrt (58 : ℝ))
    n58_exact_joint_min_7
    n58_full_prefix_residual_max_matches_exact_prefix_probe_on_face


theorem n58_prefix_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n58FaceIndices
      n58PrefixRank n58MassRank := by
  refine ⟨2, 7, n58_prefix_min_2, n58_mass_min_7, ?_⟩
  native_decide

theorem n58_prefix_mass_split_gives_two_exposed_indices :
    2 ≤ n58FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n58_prefix_mass_field_split_full_exported_face


theorem n58_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n58FaceIndices
      n58PrefixRank n58ExactMassTwice := by
  refine ⟨2, 7, n58_prefix_min_2, n58_exact_mass_min_7, ?_⟩
  native_decide

theorem n58_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n58FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n58_prefix_exact_mass_field_split_full_exported_face


theorem n58_exact_prefix_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n58FaceIndices
      n58ExactPrefixProbe n58ExactMassTwice := by
  refine ⟨2, 7, n58_exact_prefix_probe_min_2, n58_exact_mass_min_7, ?_⟩
  native_decide

theorem n58_exact_prefix_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n58FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n58_exact_prefix_exact_mass_field_split_full_exported_face


theorem n58_full_prefix_residual_max_exact_mass_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n58FaceIndices
      n58FullPrefixResidualMax n58ExactMassTwice := by
  refine ⟨2, 7, n58_full_prefix_residual_max_min_2, n58_exact_mass_min_7, ?_⟩
  native_decide

theorem n58_full_prefix_residual_max_exact_mass_split_gives_two_exposed_indices :
    2 ≤ n58FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n58_full_prefix_residual_max_exact_mass_field_split_full_exported_face


theorem n58_prefix_joint_field_split_full_exported_face :
    Erdos.Collider.FieldSplit n58FaceIndices
      n58PrefixRank n58JointRank := by
  refine ⟨2, 7, n58_prefix_min_2, n58_joint_min_7, ?_⟩
  native_decide

theorem n58_prefix_joint_split_gives_two_exposed_indices :
    2 ≤ n58FaceIndices.card :=
  Erdos.Collider.fieldSplit_card_two_le n58_prefix_joint_field_split_full_exported_face


theorem face_handoff_56_58_winner_table_matches_packet :
    n56PrefixWinnerIndicesByRank = [3] ∧
    n56MassWinnerIndicesByRank = [3] ∧
    n56JointWinnerIndicesByRank = [3] ∧
    n56ParetoIndicesByRanks = [3] ∧
    n57PrefixWinnerIndicesByRank = [4, 5] ∧
    n57MassWinnerIndicesByRank = [3] ∧
    n57JointWinnerIndicesByRank = [5] ∧
    n57ParetoIndicesByRanks = [3, 5] ∧
    n58PrefixWinnerIndicesByRank = [2, 8, 9] ∧
    n58MassWinnerIndicesByRank = [7] ∧
    n58JointWinnerIndicesByRank = [7] ∧
    n58ParetoIndicesByRanks = [7, 9] := by
  native_decide

theorem face_handoff_56_58_exact_mass_winner_table_matches_packet :
    n56MassWinnerIndicesByExactMass = [3] ∧
    n57MassWinnerIndicesByExactMass = [3] ∧
    n58MassWinnerIndicesByExactMass = [7] := by
  native_decide

theorem face_handoff_56_58_exact_mass_matches_rank_mass_winners :
    n56MassWinnerIndicesByExactMass = n56MassWinnerIndicesByRank ∧
    n57MassWinnerIndicesByExactMass = n57MassWinnerIndicesByRank ∧
    n58MassWinnerIndicesByExactMass = n58MassWinnerIndicesByRank := by
  native_decide

theorem face_handoff_56_58_split_pattern_certificate :
    n56PrefixWinnerIndicesByRank = [3] ∧
    n56MassWinnerIndicesByRank = [3] ∧
    n56JointWinnerIndicesByRank = [3] ∧
    Erdos.Collider.FieldSplit n57FaceIndices n57PrefixRank n57MassRank ∧
    Erdos.Collider.FieldSplit n57FaceIndices n57MassRank n57JointRank ∧
    Erdos.Collider.FieldSplit n58FaceIndices n58PrefixRank n58MassRank := by
  exact ⟨n56_prefix_winner_match_packet, ⟨n56_mass_winner_match_packet, ⟨n56_joint_winner_match_packet, ⟨n57_prefix_mass_field_split_full_exported_face, ⟨n57_mass_joint_field_split_full_exported_face, n58_prefix_mass_field_split_full_exported_face⟩⟩⟩⟩⟩

theorem face_handoff_56_58_exact_prefix_probe_winner_table_matches_packet :
    n56ExactPrefixProbe 3 = 0 ∧
    n57ExactPrefixProbe 4 = 0 ∧
    n57ExactPrefixProbe 5 = 0 ∧
    n58ExactPrefixProbe 2 = 0 ∧
    n58ExactPrefixProbe 8 = 0 ∧
    n58ExactPrefixProbe 9 = 0 := by
  exact ⟨n56_exact_prefix_probe_3_zero, ⟨n57_exact_prefix_probe_4_zero, ⟨n57_exact_prefix_probe_5_zero, ⟨n58_exact_prefix_probe_2_zero, ⟨n58_exact_prefix_probe_8_zero, n58_exact_prefix_probe_9_zero⟩⟩⟩⟩⟩

theorem face_handoff_56_58_exact_prefix_probe_strict_nonwinner_certificate :
    (0 < n56ExactPrefixProbe 0 ∧
    0 < n56ExactPrefixProbe 1 ∧
    0 < n56ExactPrefixProbe 2) ∧
    (0 < n57ExactPrefixProbe 0 ∧
    0 < n57ExactPrefixProbe 1 ∧
    0 < n57ExactPrefixProbe 2 ∧
    0 < n57ExactPrefixProbe 3) ∧
    (0 < n58ExactPrefixProbe 0 ∧
    0 < n58ExactPrefixProbe 1 ∧
    0 < n58ExactPrefixProbe 3 ∧
    0 < n58ExactPrefixProbe 4 ∧
    0 < n58ExactPrefixProbe 5 ∧
    0 < n58ExactPrefixProbe 6 ∧
    0 < n58ExactPrefixProbe 7) := by
  exact ⟨n56_exact_prefix_probe_positive_on_nonwinners,
    ⟨n57_exact_prefix_probe_positive_on_nonwinners,
      n58_exact_prefix_probe_positive_on_nonwinners⟩⟩

theorem face_handoff_56_58_full_prefix_zero_certificate :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 56 n56W3 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W4 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W5 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W2 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W8 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W9 := by
  exact ⟨n56W3_full_prefix_residual_zero, ⟨n57W4_full_prefix_residual_zero, ⟨n57W5_full_prefix_residual_zero, ⟨n58W2_full_prefix_residual_zero, ⟨n58W8_full_prefix_residual_zero, n58W9_full_prefix_residual_zero⟩⟩⟩⟩⟩

theorem face_handoff_56_58_packet_probe_separation_is_actual_residual_certificate :
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 (n56PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 (n56PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 (n56PacketPrefixLocation 2) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 (n57PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 (n57PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 (n57PacketPrefixLocation 2) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 (n57PacketPrefixLocation 3) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 (n58PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 (n58PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 (n58PacketPrefixLocation 3) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 (n58PacketPrefixLocation 4) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 (n58PacketPrefixLocation 5) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 (n58PacketPrefixLocation 6) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 (n58PacketPrefixLocation 7) := by
  exact ⟨n56_prefix_residual_at_packet_location_0_positive, ⟨n56_prefix_residual_at_packet_location_1_positive, ⟨n56_prefix_residual_at_packet_location_2_positive, ⟨n57_prefix_residual_at_packet_location_0_positive, ⟨n57_prefix_residual_at_packet_location_1_positive, ⟨n57_prefix_residual_at_packet_location_2_positive, ⟨n57_prefix_residual_at_packet_location_3_positive, ⟨n58_prefix_residual_at_packet_location_0_positive, ⟨n58_prefix_residual_at_packet_location_1_positive, ⟨n58_prefix_residual_at_packet_location_3_positive, ⟨n58_prefix_residual_at_packet_location_4_positive, ⟨n58_prefix_residual_at_packet_location_5_positive, ⟨n58_prefix_residual_at_packet_location_6_positive, n58_prefix_residual_at_packet_location_7_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem face_handoff_56_58_full_prefix_zero_and_probe_separation_certificate :
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 56 n56W3 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W4 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 57 n57W5 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W2 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W8 ∧
    Erdos30FaceFieldExactObservables.fullPrefixResidualZero 58 n58W9 ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W0 (n56PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W1 (n56PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 56 n56W2 (n56PacketPrefixLocation 2) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W0 (n57PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W1 (n57PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W2 (n57PacketPrefixLocation 2) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 57 n57W3 (n57PacketPrefixLocation 3) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W0 (n58PacketPrefixLocation 0) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W1 (n58PacketPrefixLocation 1) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W3 (n58PacketPrefixLocation 3) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W4 (n58PacketPrefixLocation 4) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W5 (n58PacketPrefixLocation 5) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W6 (n58PacketPrefixLocation 6) ∧
    0 < Erdos30FaceFieldExactObservables.prefixResidualAt 58 n58W7 (n58PacketPrefixLocation 7) := by
  exact ⟨n56W3_full_prefix_residual_zero, ⟨n57W4_full_prefix_residual_zero, ⟨n57W5_full_prefix_residual_zero, ⟨n58W2_full_prefix_residual_zero, ⟨n58W8_full_prefix_residual_zero, ⟨n58W9_full_prefix_residual_zero, ⟨n56_prefix_residual_at_packet_location_0_positive, ⟨n56_prefix_residual_at_packet_location_1_positive, ⟨n56_prefix_residual_at_packet_location_2_positive, ⟨n57_prefix_residual_at_packet_location_0_positive, ⟨n57_prefix_residual_at_packet_location_1_positive, ⟨n57_prefix_residual_at_packet_location_2_positive, ⟨n57_prefix_residual_at_packet_location_3_positive, ⟨n58_prefix_residual_at_packet_location_0_positive, ⟨n58_prefix_residual_at_packet_location_1_positive, ⟨n58_prefix_residual_at_packet_location_3_positive, ⟨n58_prefix_residual_at_packet_location_4_positive, ⟨n58_prefix_residual_at_packet_location_5_positive, ⟨n58_prefix_residual_at_packet_location_6_positive, n58_prefix_residual_at_packet_location_7_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem face_handoff_56_58_full_prefix_scalar_zero_certificate :
    n56FullPrefixResidualMax 3 = 0 ∧
    n57FullPrefixResidualMax 4 = 0 ∧
    n57FullPrefixResidualMax 5 = 0 ∧
    n58FullPrefixResidualMax 2 = 0 ∧
    n58FullPrefixResidualMax 8 = 0 ∧
    n58FullPrefixResidualMax 9 = 0 := by
  exact ⟨n56_full_prefix_residual_max_3_zero, ⟨n57_full_prefix_residual_max_4_zero, ⟨n57_full_prefix_residual_max_5_zero, ⟨n58_full_prefix_residual_max_2_zero, ⟨n58_full_prefix_residual_max_8_zero, n58_full_prefix_residual_max_9_zero⟩⟩⟩⟩⟩

theorem face_handoff_56_58_full_prefix_scalar_positive_certificate :
    0 < n56FullPrefixResidualMax 0 ∧
    0 < n56FullPrefixResidualMax 1 ∧
    0 < n56FullPrefixResidualMax 2 ∧
    0 < n57FullPrefixResidualMax 0 ∧
    0 < n57FullPrefixResidualMax 1 ∧
    0 < n57FullPrefixResidualMax 2 ∧
    0 < n57FullPrefixResidualMax 3 ∧
    0 < n58FullPrefixResidualMax 0 ∧
    0 < n58FullPrefixResidualMax 1 ∧
    0 < n58FullPrefixResidualMax 3 ∧
    0 < n58FullPrefixResidualMax 4 ∧
    0 < n58FullPrefixResidualMax 5 ∧
    0 < n58FullPrefixResidualMax 6 ∧
    0 < n58FullPrefixResidualMax 7 := by
  exact ⟨n56_full_prefix_residual_max_0_positive, ⟨n56_full_prefix_residual_max_1_positive, ⟨n56_full_prefix_residual_max_2_positive, ⟨n57_full_prefix_residual_max_0_positive, ⟨n57_full_prefix_residual_max_1_positive, ⟨n57_full_prefix_residual_max_2_positive, ⟨n57_full_prefix_residual_max_3_positive, ⟨n58_full_prefix_residual_max_0_positive, ⟨n58_full_prefix_residual_max_1_positive, ⟨n58_full_prefix_residual_max_3_positive, ⟨n58_full_prefix_residual_max_4_positive, ⟨n58_full_prefix_residual_max_5_positive, ⟨n58_full_prefix_residual_max_6_positive, n58_full_prefix_residual_max_7_positive⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem face_handoff_56_58_full_prefix_scalar_min_certificate :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56FullPrefixResidualMax 3 ∧
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57FullPrefixResidualMax 4 ∧
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58FullPrefixResidualMax 2 := by
  exact ⟨n56_full_prefix_residual_max_min_3, ⟨n57_full_prefix_residual_max_min_4, n58_full_prefix_residual_max_min_2⟩⟩

theorem face_handoff_56_58_full_prefix_scalar_exact_mass_split_certificate :
    Erdos.Collider.FieldSplit n57FaceIndices n57FullPrefixResidualMax n57ExactMassTwice ∧
    Erdos.Collider.FieldSplit n58FaceIndices n58FullPrefixResidualMax n58ExactMassTwice := by
  exact ⟨n57_full_prefix_residual_max_exact_mass_field_split_full_exported_face, n58_full_prefix_residual_max_exact_mass_field_split_full_exported_face⟩

theorem face_handoff_56_58_full_prefix_scalar_matches_packet_probe_certificate :
    n56FullPrefixResidualMax 0 = n56ExactPrefixProbe 0 ∧
    n56FullPrefixResidualMax 1 = n56ExactPrefixProbe 1 ∧
    n56FullPrefixResidualMax 2 = n56ExactPrefixProbe 2 ∧
    n56FullPrefixResidualMax 3 = n56ExactPrefixProbe 3 ∧
    n57FullPrefixResidualMax 0 = n57ExactPrefixProbe 0 ∧
    n57FullPrefixResidualMax 1 = n57ExactPrefixProbe 1 ∧
    n57FullPrefixResidualMax 2 = n57ExactPrefixProbe 2 ∧
    n57FullPrefixResidualMax 3 = n57ExactPrefixProbe 3 ∧
    n57FullPrefixResidualMax 4 = n57ExactPrefixProbe 4 ∧
    n57FullPrefixResidualMax 5 = n57ExactPrefixProbe 5 ∧
    n58FullPrefixResidualMax 0 = n58ExactPrefixProbe 0 ∧
    n58FullPrefixResidualMax 1 = n58ExactPrefixProbe 1 ∧
    n58FullPrefixResidualMax 2 = n58ExactPrefixProbe 2 ∧
    n58FullPrefixResidualMax 3 = n58ExactPrefixProbe 3 ∧
    n58FullPrefixResidualMax 4 = n58ExactPrefixProbe 4 ∧
    n58FullPrefixResidualMax 5 = n58ExactPrefixProbe 5 ∧
    n58FullPrefixResidualMax 6 = n58ExactPrefixProbe 6 ∧
    n58FullPrefixResidualMax 7 = n58ExactPrefixProbe 7 ∧
    n58FullPrefixResidualMax 8 = n58ExactPrefixProbe 8 ∧
    n58FullPrefixResidualMax 9 = n58ExactPrefixProbe 9 := by
  exact ⟨n56_full_prefix_residual_max_0_matches_exact_prefix_probe, ⟨n56_full_prefix_residual_max_1_matches_exact_prefix_probe, ⟨n56_full_prefix_residual_max_2_matches_exact_prefix_probe, ⟨n56_full_prefix_residual_max_3_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_0_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_1_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_2_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_3_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_4_matches_exact_prefix_probe, ⟨n57_full_prefix_residual_max_5_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_0_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_1_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_2_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_3_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_4_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_5_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_6_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_7_matches_exact_prefix_probe, ⟨n58_full_prefix_residual_max_8_matches_exact_prefix_probe, n58_full_prefix_residual_max_9_matches_exact_prefix_probe⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem face_handoff_56_58_scalar_joint_matches_probe_joint_certificate :
    n56ScalarFullPrefixJointKey 0 = n56ExactJointKey 0 ∧
    n56ScalarFullPrefixJointKey 1 = n56ExactJointKey 1 ∧
    n56ScalarFullPrefixJointKey 2 = n56ExactJointKey 2 ∧
    n56ScalarFullPrefixJointKey 3 = n56ExactJointKey 3 ∧
    n57ScalarFullPrefixJointKey 0 = n57ExactJointKey 0 ∧
    n57ScalarFullPrefixJointKey 1 = n57ExactJointKey 1 ∧
    n57ScalarFullPrefixJointKey 2 = n57ExactJointKey 2 ∧
    n57ScalarFullPrefixJointKey 3 = n57ExactJointKey 3 ∧
    n57ScalarFullPrefixJointKey 4 = n57ExactJointKey 4 ∧
    n57ScalarFullPrefixJointKey 5 = n57ExactJointKey 5 ∧
    n58ScalarFullPrefixJointKey 0 = n58ExactJointKey 0 ∧
    n58ScalarFullPrefixJointKey 1 = n58ExactJointKey 1 ∧
    n58ScalarFullPrefixJointKey 2 = n58ExactJointKey 2 ∧
    n58ScalarFullPrefixJointKey 3 = n58ExactJointKey 3 ∧
    n58ScalarFullPrefixJointKey 4 = n58ExactJointKey 4 ∧
    n58ScalarFullPrefixJointKey 5 = n58ExactJointKey 5 ∧
    n58ScalarFullPrefixJointKey 6 = n58ExactJointKey 6 ∧
    n58ScalarFullPrefixJointKey 7 = n58ExactJointKey 7 ∧
    n58ScalarFullPrefixJointKey 8 = n58ExactJointKey 8 ∧
    n58ScalarFullPrefixJointKey 9 = n58ExactJointKey 9 := by
  exact ⟨n56_scalar_full_prefix_joint_key_0_matches_exact_joint_key, ⟨n56_scalar_full_prefix_joint_key_1_matches_exact_joint_key, ⟨n56_scalar_full_prefix_joint_key_2_matches_exact_joint_key, ⟨n56_scalar_full_prefix_joint_key_3_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_0_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_1_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_2_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_3_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_4_matches_exact_joint_key, ⟨n57_scalar_full_prefix_joint_key_5_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_0_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_1_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_2_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_3_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_4_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_5_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_6_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_7_matches_exact_joint_key, ⟨n58_scalar_full_prefix_joint_key_8_matches_exact_joint_key, n58_scalar_full_prefix_joint_key_9_matches_exact_joint_key⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩⟩

theorem face_handoff_56_58_exact_joint_winner_certificate :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ExactJointKey 3 ∧
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ExactJointKey 5 ∧
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ExactJointKey 7 := by
  exact ⟨n56_exact_joint_min_3,
    ⟨n57_exact_joint_min_5,
      n58_exact_joint_min_7⟩⟩

theorem face_handoff_56_58_scalar_full_prefix_joint_winner_certificate :
    Erdos.Collider.IsFieldMinOn n56FaceIndices n56ScalarFullPrefixJointKey 3 ∧
    Erdos.Collider.IsFieldMinOn n57FaceIndices n57ScalarFullPrefixJointKey 5 ∧
    Erdos.Collider.IsFieldMinOn n58FaceIndices n58ScalarFullPrefixJointKey 7 := by
  exact ⟨n56_scalar_full_prefix_joint_min_3, ⟨n57_scalar_full_prefix_joint_min_5, n58_scalar_full_prefix_joint_min_7⟩⟩

theorem face_handoff_56_58_exact_mass_split_pattern_certificate :
    n56MassWinnerIndicesByExactMass = [3] ∧
    Erdos.Collider.FieldSplit n57FaceIndices n57PrefixRank n57ExactMassTwice ∧
    Erdos.Collider.FieldSplit n57FaceIndices n57ExactMassTwice n57JointRank ∧
    Erdos.Collider.FieldSplit n58FaceIndices n58PrefixRank n58ExactMassTwice := by
  exact ⟨n56_exact_mass_winner_match_packet, ⟨n57_prefix_exact_mass_field_split_full_exported_face, ⟨n57_exact_mass_joint_field_split_full_exported_face, n58_prefix_exact_mass_field_split_full_exported_face⟩⟩⟩

theorem face_handoff_56_58_exact_prefix_probe_exact_mass_split_certificate :
    n56MassWinnerIndicesByExactMass = [3] ∧
    Erdos.Collider.FieldSplit n57FaceIndices n57ExactPrefixProbe n57ExactMassTwice ∧
    Erdos.Collider.FieldSplit n58FaceIndices n58ExactPrefixProbe n58ExactMassTwice := by
  exact ⟨n56_exact_mass_winner_match_packet, ⟨n57_exact_prefix_exact_mass_field_split_full_exported_face, n58_exact_prefix_exact_mass_field_split_full_exported_face⟩⟩

end Erdos30FaceField5658FullFaceCertificate
