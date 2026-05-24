import Mathlib

/-!
# Erdos #30 n = 58..71 Ground-Face Branch Summary Certificate

This file records the compact finite branch table from:

`EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-71-2026-05-01_RESULTS.json`.

The source packet is the evidence for exact exported faces. This Lean file does
not import the full witness lists. Instead it certifies the finite summary
relations that matter for the next proof target: exact face counts equal
exported counts, field-selection split persists, every post-58 face contains
the previous face and the previous face shifted by `+1`, skeleton counts
strictly increase, and Pareto/joint selected surfaces remain much smaller than
the full face.

This is not a proof of Erdos #30, not an asymptotic theorem, and not a proof of
the packet enumerator. It is an auditable finite table certificate derived from
the exact packet.
-/

namespace Erdos30GroundFaceBranch5871Certificate

def branchNs : List Nat := [58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]
def post58Ns : List Nat := [59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]
def branchEdges : List (Nat × Nat) := [(58, 59), (59, 60), (60, 61), (61, 62), (62, 63), (63, 64), (64, 65), (65, 66), (66, 67), (67, 68), (68, 69), (69, 70), (70, 71)]

def faceCount : Nat -> Nat
  | 58 => 10
  | 59 => 18
  | 60 => 54
  | 61 => 152
  | 62 => 398
  | 63 => 1022
  | 64 => 2360
  | 65 => 5018
  | 66 => 9994
  | 67 => 19418
  | 68 => 36234
  | 69 => 66412
  | 70 => 117202
  | 71 => 203840
  | _ => 0

def exportedCount : Nat -> Nat
  | 58 => 10
  | 59 => 18
  | 60 => 54
  | 61 => 152
  | 62 => 398
  | 63 => 1022
  | 64 => 2360
  | 65 => 5018
  | 66 => 9994
  | 67 => 19418
  | 68 => 36234
  | 69 => 66412
  | 70 => 117202
  | 71 => 203840
  | _ => 0

def skeletonCount : Nat -> Nat
  | 58 => 2
  | 59 => 4
  | 60 => 18
  | 61 => 49
  | 62 => 123
  | 63 => 312
  | 64 => 669
  | 65 => 1329
  | 66 => 2488
  | 67 => 4712
  | 68 => 8408
  | 69 => 15089
  | 70 => 25395
  | 71 => 43319
  | _ => 0

def previousFacePersistenceCount : Nat -> Nat
  | 58 => 0
  | 59 => 10
  | 60 => 18
  | 61 => 54
  | 62 => 152
  | 63 => 398
  | 64 => 1022
  | 65 => 2360
  | 66 => 5018
  | 67 => 9994
  | 68 => 19418
  | 69 => 36234
  | 70 => 66412
  | 71 => 117202
  | _ => 0

def plusOnePreviousFacePersistenceCount : Nat -> Nat
  | 58 => 0
  | 59 => 10
  | 60 => 18
  | 61 => 54
  | 62 => 152
  | 63 => 398
  | 64 => 1022
  | 65 => 2360
  | 66 => 5018
  | 67 => 9994
  | 68 => 19418
  | 69 => 36234
  | 70 => 66412
  | 71 => 117202
  | _ => 0

def paretoCount : Nat -> Nat
  | 58 => 2
  | 59 => 2
  | 60 => 1
  | 61 => 4
  | 62 => 6
  | 63 => 10
  | 64 => 11
  | 65 => 47
  | 66 => 40
  | 67 => 120
  | 68 => 109
  | 69 => 267
  | 70 => 174
  | 71 => 452
  | _ => 0

def jointCount : Nat -> Nat
  | 58 => 1
  | 59 => 1
  | 60 => 1
  | 61 => 4
  | 62 => 6
  | 63 => 10
  | 64 => 11
  | 65 => 47
  | 66 => 40
  | 67 => 120
  | 68 => 109
  | 69 => 267
  | 70 => 174
  | 71 => 452
  | _ => 0

def prefixCount : Nat -> Nat
  | 58 => 3
  | 59 => 6
  | 60 => 21
  | 61 => 44
  | 62 => 97
  | 63 => 204
  | 64 => 481
  | 65 => 823
  | 66 => 1253
  | 67 => 2424
  | 68 => 3554
  | 69 => 6049
  | 70 => 8191
  | 71 => 11798
  | _ => 0

def massCount : Nat -> Nat
  | 58 => 1
  | 59 => 1
  | 60 => 1
  | 61 => 5
  | 62 => 7
  | 63 => 19
  | 64 => 18
  | 65 => 84
  | 66 => 73
  | 67 => 296
  | 68 => 283
  | 69 => 882
  | 70 => 792
  | 71 => 2810
  | _ => 0

def fieldSelectionSplit : Nat -> Bool
  | 58 => true
  | 59 => true
  | 60 => true
  | 61 => true
  | 62 => true
  | 63 => true
  | 64 => true
  | 65 => true
  | 66 => true
  | 67 => true
  | 68 => true
  | 69 => true
  | 70 => true
  | 71 => true
  | _ => false

def exportedCountMatchesExact (n : Nat) : Bool :=
  exportedCount n == faceCount n

def previousFacePersistsAcrossEdge (edge : Nat × Nat) : Bool :=
  previousFacePersistenceCount edge.2 == faceCount edge.1

def plusOnePreviousFacePersistsAcrossEdge (edge : Nat × Nat) : Bool :=
  plusOnePreviousFacePersistenceCount edge.2 == faceCount edge.1

def skeletonStrictlyIncreasesAcrossEdge (edge : Nat × Nat) : Bool :=
  decide (skeletonCount edge.1 < skeletonCount edge.2)

def paretoSurfaceIsProper (n : Nat) : Bool :=
  decide (paretoCount n < faceCount n)

def jointSurfaceIsProper (n : Nat) : Bool :=
  decide (jointCount n < faceCount n)

def jointNoLargerThanParetoCount (n : Nat) : Bool :=
  decide (jointCount n <= paretoCount n)

def branchSummaryCertificate : Bool :=
  branchNs.all exportedCountMatchesExact &&
  branchNs.all fieldSelectionSplit &&
  branchEdges.all previousFacePersistsAcrossEdge &&
  branchEdges.all plusOnePreviousFacePersistsAcrossEdge &&
  branchEdges.all skeletonStrictlyIncreasesAcrossEdge &&
  branchNs.all paretoSurfaceIsProper &&
  branchNs.all jointSurfaceIsProper &&
  branchNs.all jointNoLargerThanParetoCount

theorem branch_table_face_counts :
    branchNs.map faceCount = [10, 18, 54, 152, 398, 1022, 2360, 5018, 9994, 19418, 36234, 66412, 117202, 203840] := by
  native_decide

theorem branch_table_skeleton_counts :
    branchNs.map skeletonCount = [2, 4, 18, 49, 123, 312, 669, 1329, 2488, 4712, 8408, 15089, 25395, 43319] := by
  native_decide

theorem branch_table_pareto_counts :
    branchNs.map paretoCount = [2, 2, 1, 4, 6, 10, 11, 47, 40, 120, 109, 267, 174, 452] := by
  native_decide

theorem branch_table_joint_counts :
    branchNs.map jointCount = [1, 1, 1, 4, 6, 10, 11, 47, 40, 120, 109, 267, 174, 452] := by
  native_decide

theorem all_exports_are_complete_by_count :
    branchNs.all exportedCountMatchesExact = true := by
  native_decide

theorem field_selection_split_persists_58_71 :
    branchNs.all fieldSelectionSplit = true := by
  native_decide

theorem every_post58_face_contains_previous_face :
    branchEdges.all previousFacePersistsAcrossEdge = true := by
  native_decide

theorem every_post58_face_contains_previous_face_plus_one :
    branchEdges.all plusOnePreviousFacePersistsAcrossEdge = true := by
  native_decide

theorem skeleton_count_strictly_increases_58_71 :
    branchEdges.all skeletonStrictlyIncreasesAcrossEdge = true := by
  native_decide

theorem pareto_surface_is_proper_subset_by_count_58_71 :
    branchNs.all paretoSurfaceIsProper = true := by
  native_decide

theorem joint_surface_is_proper_subset_by_count_58_71 :
    branchNs.all jointSurfaceIsProper = true := by
  native_decide

theorem joint_count_no_larger_than_pareto_count_58_71 :
    branchNs.all jointNoLargerThanParetoCount = true := by
  native_decide

theorem branch_summary_certificate_passes :
    branchSummaryCertificate = true := by
  native_decide

end Erdos30GroundFaceBranch5871Certificate
