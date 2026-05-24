import Lake
open Lake DSL

package erdos30_sidon where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

-- ══════════════════════════════════════════════════════════════
-- Core formalization (discussed in paper)
-- ══════════════════════════════════════════════════════════════

@[default_target]
lean_lib Erdos30_Sidon_Defs where
  srcDir := "lean"
  roots := #[`Erdos30_Sidon_Defs]

lean_lib Erdos30_Lindstrom where
  srcDir := "lean"
  roots := #[`Erdos30_Lindstrom]

lean_lib Erdos30_BFR where
  srcDir := "lean"
  roots := #[`Erdos30_BFR]

lean_lib Erdos30_Singer where
  srcDir := "lean"
  roots := #[`Erdos30_Singer]

lean_lib Erdos30_Complete where
  srcDir := "lean"
  roots := #[`Erdos30_Complete]

lean_lib Erdos30_OrderedElements where
  srcDir := "lean"
  roots := #[`Erdos30_OrderedElements]

lean_lib Erdos30_FaceField where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField]

lean_lib Erdos30_FaceField_ExactObservables where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_ExactObservables]

lean_lib Erdos30_FaceField_N30_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_N30_Certificate]

lean_lib Erdos30_FaceField_Window_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_Window_Certificate]

lean_lib Erdos30_FaceField_57_58_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_57_58_Certificate]

lean_lib Erdos30_FaceField_56_58_FullFace_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_56_58_FullFace_Certificate]

lean_lib Erdos30_FaceField_59_FullFace_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_59_FullFace_Certificate]

lean_lib Erdos30_FaceField_60_FullFace_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_60_FullFace_Certificate]

lean_lib Erdos30_FaceField_61_JointSurface_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_61_JointSurface_Certificate]

lean_lib Erdos30_FaceField_61_64_JointSurface_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_61_64_JointSurface_Certificate]

lean_lib Erdos30_FaceField_61_64_JointSurface_Transition_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_61_64_JointSurface_Transition_Certificate]

lean_lib Erdos30_FaceField_65_71_JointSurface_Transition_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_FaceField_65_71_JointSurface_Transition_Certificate]

lean_lib Erdos30_GroundFaceBranch_58_71_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_GroundFaceBranch_58_71_Certificate]

-- ══════════════════════════════════════════════════════════════
-- Scratch / supplementary files (not discussed in paper,
-- not imported by core files — kept for reference only)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos30_difference_counting where
  srcDir := "scratch"
  roots := #[`Erdos30_difference_counting]

lean_lib Sidon_SumCount_Fix where
  srcDir := "scratch"
  roots := #[`Sidon_SumCount_Fix]

lean_lib Ehp114LocalMixedRemainderScratch where
  srcDir := "scratch"
  roots := #[`Ehp114LocalMixedRemainderScratch]

-- ══════════════════════════════════════════════════════════════
-- Erdős #755 — B_h[g] Sequences (salvo attack 2026-04-19)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos755_BhG where
  srcDir := "lean"
  roots := #[`Erdos755_BhG]

lean_lib Erdos755_DifferenceCount where
  srcDir := "lean"
  roots := #[`Erdos755_DifferenceCount]

lean_lib Erdos755_Lindstrom where
  srcDir := "lean"
  roots := #[`Erdos755_Lindstrom]

lean_lib Erdos755_Singer_BhG where
  srcDir := "lean"
  roots := #[`Erdos755_Singer_BhG]

lean_lib Erdos755_Complete where
  srcDir := "lean"
  roots := #[`Erdos755_Complete]

lean_lib Erdos755_B3G where
  srcDir := "lean"
  roots := #[`Erdos755_B3G]

lean_lib Erdos755_BhG_General where
  srcDir := "lean"
  roots := #[`Erdos755_BhG_General]

lean_lib Erdos755_SymmetryQuotient where
  srcDir := "lean"
  roots := #[`Erdos755_SymmetryQuotient]

-- ══════════════════════════════════════════════════════════════
-- Erdős #1 — Distinct Subset Sums (salvo attack 2026-04-19)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos1_DistinctSubsetSums where
  srcDir := "lean"
  roots := #[`Erdos1_DistinctSubsetSums]

-- ══════════════════════════════════════════════════════════════
-- Erdős #755 — Higher-order specializations (h=4,5,6)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos755_HigherOrder where
  srcDir := "lean"
  roots := #[`Erdos755_HigherOrder]

-- ══════════════════════════════════════════════════════════════
-- Erdős #166 — Sum-Free Sets (salvo attack 2026-04-19)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos166_SumFree where
  srcDir := "lean"
  roots := #[`Erdos166_SumFree]

-- ══════════════════════════════════════════════════════════════
-- Erdős #30 — Sharp Sidon Difference Bound (PMF salvo 2026-04-19)
-- ══════════════════════════════════════════════════════════════

lean_lib Erdos30_SharpDiff where
  srcDir := "lean"
  roots := #[`Erdos30_SharpDiff]

lean_lib Erdos30_Singer57_Certificate where
  srcDir := "lean"
  roots := #[`Erdos30_Singer57_Certificate]

lean_lib Erdos30_CollisionChannel where
  srcDir := "lean"
  roots := #[`Erdos30_CollisionChannel]

-- ══════════════════════════════════════════════════════════════
-- Erdős #114 — EHP radial-direction Athena spike (2026-05-02)
-- Scaffold for n=14 closed-form theorem; namespace Erdos114.Radial
-- See: erdos-experiments/Erdos114/ATHENA_SPIKE_PROTOCOL_2026-05-02.md
-- ══════════════════════════════════════════════════════════════

lean_lib EhpRadialPuiseux where
  srcDir := "lean"
  roots := #[`EhpRadialPuiseux]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.27.0"
