/-
  Erdos #30 -- Singer-modulus finite certificate
  =================================================

  This file formalizes finite, packet-backed checks from:

  * EXP-MM-030-PMF-SINGER57-CERTIFICATE-V2-2026-04-30
  * EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30

  It is deliberately not an asymptotic theorem about Erdos #30.  The claims
  below are finite exact-face certificates checked by decidable enumeration.
-/

import Mathlib

namespace Erdos30Singer57Certificate

def containsNat (x : Nat) : List Nat -> Bool
  | [] => false
  | y :: ys => if x == y then true else containsNat x ys

def countNat (x : Nat) : List Nat -> Nat
  | [] => 0
  | y :: ys => (if x == y then 1 else 0) + countNat x ys

def noDupsNat : List Nat -> Bool
  | [] => true
  | x :: xs => if containsNat x xs then false else noDupsNat xs

def allNat (p : Nat -> Bool) : List Nat -> Bool
  | [] => true
  | x :: xs => p x && allNat p xs

def allListNat (p : List Nat -> Bool) : List (List Nat) -> Bool
  | [] => true
  | x :: xs => p x && allListNat p xs

def anyNat (p : Nat -> Bool) : List Nat -> Bool
  | [] => false
  | x :: xs => p x || anyNat p xs

def anyListNat (p : List Nat -> Bool) : List (List Nat) -> Bool
  | [] => false
  | x :: xs => p x || anyListNat p xs

def positiveDiffs : List Nat -> List Nat
  | [] => []
  | a :: rest => rest.map (fun b => b - a) ++ positiveDiffs rest

def orderedModDiffsFrom (m : Nat) (a : Nat) : List Nat -> List Nat
  | [] => []
  | b :: bs =>
      if a == b then
        orderedModDiffsFrom m a bs
      else
        ((a + m - b) % m) :: orderedModDiffsFrom m a bs

def orderedModDiffsAux (m : Nat) (full : List Nat) : List Nat -> List Nat
  | [] => []
  | a :: rest => orderedModDiffsFrom m a full ++ orderedModDiffsAux m full rest

def orderedModDiffs (m : Nat) (xs : List Nat) : List Nat :=
  orderedModDiffsAux m xs xs

def nonzeroResidues (m : Nat) : List Nat :=
  (List.range m).drop 1

def sameNatSet (xs : List Nat) (ys : List Nat) : Bool :=
  allNat (fun x => containsNat x ys) xs && allNat (fun y => containsNat y xs) ys

def intervalSidon10 (xs : List Nat) : Bool :=
  (xs.length == 10) && ((positiveDiffs xs).length == 45) && noDupsNat (positiveDiffs xs)

def coversAllNonzeroResidues (m : Nat) (xs : List Nat) : Bool :=
  allNat (fun r => containsNat r (orderedModDiffs m xs)) (nonzeroResidues m)

def perfectDifferenceSet (m : Nat) (xs : List Nat) : Bool :=
  noDupsNat xs &&
  ((orderedModDiffs m xs).length == m - 1) &&
  allNat (fun r => countNat r (orderedModDiffs m xs) == 1) (nonzeroResidues m)

def translate (shift : Nat) (xs : List Nat) : List Nat :=
  xs.map (fun x => x + shift)

def residueImage (m : Nat) (multiplier : Nat) (shift : Nat) (xs : List Nat) : List Nat :=
  xs.map (fun x => (multiplier * x + shift) % m)

def subsetResidues (needle : List Nat) (haystack : List Nat) : Bool :=
  allNat (fun x => containsNat x (haystack.map (fun y => y % 57))) needle

def hasAffinePdsImageInWitness (pds : List Nat) (witness : List Nat) : Bool :=
  anyNat
    (fun multiplier =>
      if Nat.gcd multiplier 57 == 1 then
        anyNat
          (fun shift => subsetResidues (residueImage 57 multiplier shift pds) witness)
          (List.range 57)
      else
        false)
    (List.range 57)

def anyAffinePdsImageContained (pdsSets : List (List Nat)) (witnesses : List (List Nat)) : Bool :=
  anyListNat
    (fun pds =>
      anyListNat (fun witness => hasAffinePdsImageInWitness pds witness) witnesses)
    pdsSets

/-! ## Packet-backed witnesses at n = 57 -/

def N57_0 : List Nat := [0, 1, 6, 10, 23, 26, 34, 41, 53, 55]
def N57_1 : List Nat := [0, 2, 14, 21, 29, 32, 45, 49, 54, 55]
def N57_2 : List Nat := [1, 2, 7, 11, 24, 27, 35, 42, 54, 56]
def N57_3 : List Nat := [1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
def N57_4 : List Nat := [2, 3, 8, 12, 25, 28, 36, 43, 55, 57]
def N57_5 : List Nat := [2, 4, 16, 23, 31, 34, 47, 51, 56, 57]

def n57Witnesses : List (List Nat) := [N57_0, N57_1, N57_2, N57_3, N57_4, N57_5]

/-! ## Cited Singer PDS representatives, recorded in the V2 packet -/

def PDS_A : List Nat := [0, 11, 19, 20, 24, 26, 36, 54]
def PDS_B : List Nat := [0, 1, 5, 7, 17, 35, 38, 49]
def PDS_C : List Nat := [16, 19, 30, 38, 39, 43, 45, 55]

def citedPdsSets : List (List Nat) := [PDS_A, PDS_B, PDS_C]

/-! ## Packet-backed witnesses at n = 58 -/

def N58_0 : List Nat := [0, 1, 6, 10, 23, 26, 34, 41, 53, 55]
def N58_1 : List Nat := [0, 2, 14, 21, 29, 32, 45, 49, 54, 55]
def N58_2 : List Nat := [0, 2, 15, 21, 22, 32, 46, 50, 55, 58]
def N58_3 : List Nat := [0, 3, 8, 12, 26, 36, 37, 43, 56, 58]
def N58_4 : List Nat := [1, 2, 7, 11, 24, 27, 35, 42, 54, 56]
def N58_5 : List Nat := [1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
def N58_6 : List Nat := [2, 3, 8, 12, 25, 28, 36, 43, 55, 57]
def N58_7 : List Nat := [2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
def N58_8 : List Nat := [3, 4, 9, 13, 26, 29, 37, 44, 56, 58]
def N58_9 : List Nat := [3, 5, 17, 24, 32, 35, 48, 52, 57, 58]

def n58Witnesses : List (List Nat) :=
  [N58_0, N58_1, N58_2, N58_3, N58_4, N58_5, N58_6, N58_7, N58_8, N58_9]

def n58PrefixWinners : List Nat := [2, 8, 9]
def n58MassWinners : List Nat := [7]
def n58JointWinners : List Nat := [7]
def n58ParetoWinners : List Nat := [7, 9]

/-! ## C1-C6 finite certificate checks -/

theorem C1_n57_witnesses_are_interval_sidon :
    allListNat intervalSidon10 n57Witnesses = true := by
  native_decide

theorem C2_n57_witnesses_share_positive_difference_skeleton :
    allListNat (fun xs => sameNatSet (positiveDiffs xs) (positiveDiffs N57_0)) n57Witnesses = true := by
  native_decide

theorem C3_n57_representative_covers_all_nonzero_residues_mod57 :
    coversAllNonzeroResidues 57 N57_3 = true := by
  native_decide

theorem C4_n57_joint_winner_is_mass_winner_plus_one :
    N57_5 = translate 1 N57_3 := by
  native_decide

theorem C5_cited_representatives_are_57_8_1_pds :
    allListNat (perfectDifferenceSet 57) citedPdsSets = true := by
  native_decide

theorem C5_no_cited_pds_affine_image_is_contained_in_n57_witness :
    anyAffinePdsImageContained citedPdsSets n57Witnesses = false := by
  native_decide

theorem C6_n58_exported_witnesses_are_interval_sidon :
    allListNat intervalSidon10 n58Witnesses = true := by
  native_decide

theorem C6_n58_branch_is_prefix_side_not_pareto_side :
    sameNatSet (positiveDiffs N58_7) (positiveDiffs N58_9) &&
    (N58_9 == translate 1 N58_7) &&
    sameNatSet (positiveDiffs N58_2) (positiveDiffs N58_3) &&
    (sameNatSet (positiveDiffs N58_7) (positiveDiffs N58_2) == false) &&
    containsNat 2 n58PrefixWinners &&
    ((containsNat 2 n58ParetoWinners) == false) &&
    (n58MassWinners == [7]) &&
    (n58JointWinners == [7]) &&
    (n58ParetoWinners == [7, 9]) = true := by
  native_decide

theorem C6_n58_pareto_chain_persists_from_n57 :
    (N58_7 == N57_5) &&
    (N58_9 == translate 1 N57_5) = true := by
  native_decide

end Erdos30Singer57Certificate
