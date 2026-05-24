/-
Copyright (c) 2026 Kenneth A. Mendoza. All rights reserved.
Released under the same license as the formal-conjectures repository.

Stub for: FormalConjectures/ErdosProblems/114Finite.lean

This file is a proposed addition to google-deepmind/formal-conjectures.
It does NOT change the status of `Erdos114.erdos_114`. It adds a finite
variant theorem statement and records the per-degree dependencies as
explicit axioms so reviewers can audit what is being assumed.

The accompanying packet:
  EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-12

Certificate archive:
  DOI 10.5281/zenodo.19480329
  https://github.com/MendozaLab/erdos-experiments

This stub is a draft. Names, namespace, and import surface should be
adapted to upstream conventions in `formal-conjectures`. The mathematical
content of the axioms is what matters; the syntactic packaging is open
to upstream maintainers.
-/

import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.MeasureTheory.Measure.Lebesgue.EqHaar
import Mathlib.Topology.Algebra.Polynomial

namespace Erdos114

open Polynomial Complex

-- The lemniscate of a polynomial `p` is the level set `{z : |p z| = 1}`.
-- `lemniscateLength p` denotes its 1-dimensional Hausdorff measure.
-- (Definition site is upstream; this stub assumes it is in scope.)

/--
Reference upper bound: the lemniscate length of `z^n - 1`.

This serves as the right-hand side of Erdős–Herzog–Piranian.
-/
noncomputable def referenceLength (n : ℕ) : ℝ :=
  -- Placeholder; upstream definition expected.
  -- Stub uses the value of `lemniscateLength (X^n - 1 : ℂ[X])` once the
  -- ambient `lemniscateLength` is defined.
  sorry

/-! ### `n = 1`: direct analytic

A monic linear polynomial `X - a` has unit lemniscate `{|z - a| = 1}`,
the unit circle translated by `a`, with length `2 * π = referenceLength 1`.
-/

axiom lemniscateLength_le_referenceLength_of_degree_one
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 1) :
    -- lemniscateLength {z | Complex.abs (p.eval z) = 1} ≤ referenceLength 1
    True
-- Replace `True` with the actual statement once `lemniscateLength` is
-- in scope upstream.

/-! ### `n = 2`: literature row (MacLane / Eremenko-Hayman)

Source: G. R. MacLane, "On a conjecture of Erdős, Herzog, and Piranian,"
Michigan Math. J. 2 (1953/54), 147-148; doi:10.1307/mmj/1028989918.

Accessible pin: A. Eremenko and W. K. Hayman, "On the length of lemniscates,"
Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295. The abstract names
Bernoulli's lemniscate as the `d = 2` extremal level set; the proof remark
after Lemma 5 identifies `z^2 + 1` as extremal, length-equivalent to
`z^2 - 1` by rotation.
-/

axiom lemniscateLength_le_referenceLength_of_degree_two
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 2) :
    True
-- Replace `True` with the actual statement once `lemniscateLength` is in
-- scope upstream. Cite MacLane and Eremenko-Hayman in the upstream
-- docstring.

/-! ### `3 <= n <= 12` and `n = 14`: IEEE-1788 interval certificates

Each axiom corresponds to one DOI-archived Rust + `inari` interval
certificate. The SHA-256 of the result JSON is named in the docstring.
-/

axiom lemniscateLength_le_referenceLength_n3
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 3) :
    True
-- Certificate: EXP-MM-EHP-007-n3-inari_RESULTS.json
-- SHA-256: a884d1bfec1563f6e6f7ae4cbb2ec607b43be033d06c41d14782459e67ec2b95
-- Zenodo: https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n3-inari_RESULTS.json

axiom lemniscateLength_le_referenceLength_n4
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 4) :
    True
-- Certificate: EXP-MM-EHP-007-n4-inari_RESULTS.json
-- SHA-256: 0924dd7424d2615099ff95d47cb4c120ba22e907adaa9af881cded9678241209

axiom lemniscateLength_le_referenceLength_n5
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 5) :
    True
-- Certificate: EXP-MM-EHP-007-n5-inari_RESULTS.json
-- SHA-256: 21ca3c7607dc1fbb7b08982666f4620dff808fdc581eecae9967c51fafb05447

axiom lemniscateLength_le_referenceLength_n6
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 6) :
    True
-- Certificate: EXP-MM-EHP-007-n6-inari_RESULTS.json
-- SHA-256: 41ac3027e9ae5add9e1208c0faa5897c36762d69b5eeccef96068c96af567b3d

axiom lemniscateLength_le_referenceLength_n7
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 7) :
    True
-- Certificate: EXP-MM-EHP-007-n7-inari_RESULTS.json
-- SHA-256: 832ddaf219d717e275ee95c01f271dd3120e255cf81f3aa72ba3c2f56ff84054

axiom lemniscateLength_le_referenceLength_n8
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 8) :
    True
-- Certificate: EXP-MM-EHP-007-n8-inari_RESULTS.json
-- SHA-256: c7a1fd80fbfed1efd53eaa35283e467994e3a4175541f3817485d9551d14dcdc

axiom lemniscateLength_le_referenceLength_n9
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 9) :
    True
-- Certificate: EXP-MM-EHP-007-n9-inari_RESULTS.json
-- SHA-256: 5bc2887826c9ef21752c115c9a4a2ab983f94eea06f0111ea98db454fe1358f4

axiom lemniscateLength_le_referenceLength_n10
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 10) :
    True
-- Certificate: EXP-MM-EHP-007-n10-inari_RESULTS.json
-- SHA-256: 2b72e052aa7200f7ac5d40992843601988de234093c44860fde99a7871e19581

axiom lemniscateLength_le_referenceLength_n11
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 11) :
    True
-- Certificate: EXP-MM-EHP-007-n11-inari_RESULTS.json
-- SHA-256: 67f20cce1d3d54cad2d6bc708ab9ec796c17cb4d34a9728f2954b8a6cbf7c89c

axiom lemniscateLength_le_referenceLength_n12
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 12) :
    True
-- Certificate: EXP-MM-EHP-007-n12-inari_RESULTS.json
-- SHA-256: 42a517997445d158649feefae2b7287bc9b548c6391796df0c1246489c6aa064

-- TODO: n = 13 is intentionally not axiomatized. The underlying certificate
-- EXP-MM-EHP-007-n13-inari_RESULTS.json reports
--   bb_total_evals = 0
--   bb_level_count = 0
-- while neighboring rows report tens to hundreds of millions of evaluations.
-- Quarantined pending a re-run with the evaluator profile used for n != 13.

axiom lemniscateLength_le_referenceLength_n14
    (p : ℂ[X]) (hp_monic : p.Monic) (hp_deg : p.natDegree = 14) :
    True
-- Certificate: EXP-MM-EHP-007-n14-inari_RESULTS.json
-- SHA-256: 50b1c965c842ced25b2930c2b71ffb6e2da693872aa464a19fbd9d5d9efa0ca7

/-! ### Finite variant theorem

For every monic complex polynomial `p` of degree `n` with `1 <= n <= 14`
and `n != 13`, the lemniscate length is bounded by `referenceLength n`.

This is **not** a partial proof of the all-degree `erdos_114` conjecture
and is **not** in scope of Tao's sufficiently-large-`n` theorem
(arXiv:2512.12455).
-/

theorem erdos_114_finite_le_14_except_13
    (p : ℂ[X]) (hp_monic : p.Monic)
    (hp_deg_le_14 : p.natDegree ≤ 14)
    (hp_deg_pos : 0 < p.natDegree)
    (hp_deg_ne_13 : p.natDegree ≠ 13) :
    True := by
  -- Replace `True` with the upstream `lemniscateLength` statement.
  -- Proof skeleton:
  --   rcases p.natDegree with _ | _ | _ | _ | ... | 14
  --   · exact lemniscateLength_le_referenceLength_of_degree_one p hp_monic rfl
  --   · exact lemniscateLength_le_referenceLength_of_degree_two p hp_monic rfl
  --   · exact lemniscateLength_le_referenceLength_n3 p hp_monic rfl
  --   ...
  --   · exact lemniscateLength_le_referenceLength_n14 p hp_monic rfl
  trivial

end Erdos114
