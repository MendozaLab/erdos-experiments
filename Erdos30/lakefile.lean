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

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "f897ebcf72cd16f89ab4577d0c826cd14afaafc7"
