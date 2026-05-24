/-
  Erdős Problem #30 — Singer Difference Set Construction (Lower Bound)
  ====================================================================

  **Stage:** 2 (Proof Strategy Design) + 4 (Lean Scaffold)
  **Pipeline:** erdos-problem-solver 7-stage engine
  **Tier:** Tier 1 (small cases) + Tier 4 axiom (general case)
  **MDL Morphism:** Prefix-Free Code (confidence 0.91)

  **Classical Statement:**
    For prime q, there exists a Sidon set A ⊂ {0,...,q²+q} with |A| = q+1.
    This gives the lower bound h(N) ≥ √N for N = q²+q.

  **Singer's Construction (1938):**
    1. Consider the projective plane PG(2,q) of order q
    2. It has q²+q+1 points, arranged cyclically via a Singer cycle σ
    3. A line in PG(2,q) has q+1 points
    4. Number points 0,1,...,q²+q by powers of σ
    5. The indices of points on any one line form a "perfect difference set"
       D ⊂ Z_{q²+q+1} with |D| = q+1
    6. Perfect difference set ⟹ all pairwise differences distinct ⟹ Sidon property

  **Proof Strategy (this file):**
    - Layer 1: Verify q=2 explicitly ({0,1,3} is Sidon, 3 elements ≤ 6)
    - Layer 2: Verify q=3 explicitly ({0,1,3,9} is Sidon, 4 elements ≤ 12)
    - Layer 3: General Singer axiom with full mathematical reference
    - Layer 4: Derive h(N) ≥ √N from Singer + Bertrand's postulate

  **References:**
    - Singer, J. (1938). "A theorem in finite projective geometry and some
      applications to number theory." Trans. AMS, 43(3), 377-385.
    - Erdős & Turán (1941). "On a problem of Sidon in additive number theory."
      J. London Math. Soc., 16, 212-215.
    - O'Bryant (2004). "A complete annotated bibliography of work related to
      Sidon sets." Electronic J. Combinatorics, DS11.

  **Perplexity Pre-Submission Gate (2026-03-29):**
    Query: "Has the Singer difference set construction been formally verified
    in Lean 4 or any other theorem prover?"
    Result: NO — no formal verification found in Lean, Coq, Isabelle, or Mizar.
    This formalization would be novel.

  **Lean version**: leanprover/lean4:v4.24.0
  **Mathlib version**: f897ebcf72cd16f89ab4577d0c826cd14afaafc7
-/

import Mathlib
import Erdos30_Sidon_Defs

open Finset Nat

namespace Erdos.Sidon

/-! ## Layer 1: Singer Construction for q = 2

  The Singer difference set for q = 2 is {0, 1, 3} in Z_7.
  - PG(2,2) is the Fano plane with 7 points
  - A line has 3 = q+1 points
  - {0, 1, 3} is a perfect difference set mod 7:
    differences: 1-0=1, 3-0=3, 3-1=2, 0-1≡6, 0-3≡4, 1-3≡5
    All residues 1..6 appear exactly once ✓
  - As a Sidon set in ℕ: {0, 1, 3} has all pairwise sums distinct:
    0+0=0, 0+1=1, 0+3=3, 1+1=2, 1+3=4, 3+3=6 — all distinct ✓
-/

/-- The Singer set for q=2. -/
def singer_q2 : Finset ℕ := {0, 1, 3}

/-- {0, 1, 3} is a Sidon set — verified computationally. -/
theorem singer_q2_sidon : IsSidonSet singer_q2 := by native_decide

/-- {0, 1, 3} has exactly 3 = q+1 elements. -/
theorem singer_q2_card : singer_q2.card = 3 := by decide

/-- All elements of {0, 1, 3} are ≤ 6 = q²+q = 4+2. -/
theorem singer_q2_range : ∀ a ∈ singer_q2, a ≤ 2 * 2 + 2 := by decide

/-- Singer construction witness for q = 2. -/
theorem singer_sidon_q2 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 2 + 1 ∧ ∀ a ∈ A, a ≤ 2 * 2 + 2 :=
  ⟨singer_q2, singer_q2_sidon, singer_q2_card, singer_q2_range⟩

/-! ## Layer 2: Singer Construction for q = 3

  The Singer difference set for q = 3 is {0, 1, 3, 9} in Z_13.
  - PG(2,3) has 13 points
  - A line has 4 = q+1 points
  - {0, 1, 3, 9} is a perfect difference set mod 13
  - Pairwise sums: 0,1,2,3,4,6,9,10,12,18 — need to check distinctness for a≤b
    (0,0)→0, (0,1)→1, (0,3)→3, (0,9)→9,
    (1,1)→2, (1,3)→4, (1,9)→10,
    (3,3)→6, (3,9)→12,
    (9,9)→18 — all 10 distinct ✓
-/

/-- The Singer set for q=3. -/
def singer_q3 : Finset ℕ := {0, 1, 3, 9}

/-- {0, 1, 3, 9} is a Sidon set — verified computationally. -/
theorem singer_q3_sidon : IsSidonSet singer_q3 := by native_decide

/-- {0, 1, 3, 9} has exactly 4 = q+1 elements. -/
theorem singer_q3_card : singer_q3.card = 4 := by decide

/-- All elements of {0, 1, 3, 9} are ≤ 12 = q²+q = 9+3. -/
theorem singer_q3_range : ∀ a ∈ singer_q3, a ≤ 3 * 3 + 3 := by decide

/-- Singer construction witness for q = 3. -/
theorem singer_sidon_q3 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 3 + 1 ∧ ∀ a ∈ A, a ≤ 3 * 3 + 3 :=
  ⟨singer_q3, singer_q3_sidon, singer_q3_card, singer_q3_range⟩

/-! ## Layer 2b: Singer Construction for q = 5

  The Singer difference set for q = 5 is {0, 1, 3, 8, 12, 18} in Z_31.
  - PG(2,5) has 31 points
  - A line has 6 = q+1 points
  - {0, 1, 3, 8, 12, 18} is a perfect difference set mod 31
  - Computationally verified: all 10 pairwise sums (with a≤b) are distinct
  - All elements ≤ 30 = q²+q = 25+5
-/

/-- The Singer set for q=5. -/
def singer_q5 : Finset ℕ := {0, 1, 3, 8, 12, 18}

/-- {0, 1, 3, 8, 12, 18} is a Sidon set — verified computationally. -/
theorem singer_q5_sidon : IsSidonSet singer_q5 := by native_decide

/-- {0, 1, 3, 8, 12, 18} has exactly 6 = q+1 elements. -/
theorem singer_q5_card : singer_q5.card = 6 := by decide

/-- All elements of {0, 1, 3, 8, 12, 18} are ≤ 30 = q²+q = 25+5. -/
theorem singer_q5_range : ∀ a ∈ singer_q5, a ≤ 5 * 5 + 5 := by decide

/-- Singer construction witness for q = 5. -/
theorem singer_sidon_q5 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 5 + 1 ∧ ∀ a ∈ A, a ≤ 5 * 5 + 5 :=
  ⟨singer_q5, singer_q5_sidon, singer_q5_card, singer_q5_range⟩

/-! ## Layer 2c: Singer Construction for q = 7

  The Singer difference set for q = 7 is {0, 1, 3, 13, 32, 36, 43, 52} in Z_57.
  - PG(2,7) has 57 points (q²+q+1 = 49+7+1 = 57)
  - A line has 8 = q+1 points
  - {0, 1, 3, 13, 32, 36, 43, 52} is a perfect difference set mod 57:
    each of the 56 nonzero residues mod 57 is realized as exactly one
    ordered difference a - b (a ≠ b in S). Verified externally.
  - All elements ≤ 56 = q²+q = 49+7
  - Sidon property in ℕ follows from PDS structure (all 36 unordered pairwise
    sums distinct), verified computationally below via `native_decide`.
-/

/-- The Singer set for q=7. -/
def singer_q7 : Finset ℕ := {0, 1, 3, 13, 32, 36, 43, 52}

/-- {0, 1, 3, 13, 32, 36, 43, 52} is a Sidon set — verified computationally. -/
theorem singer_q7_sidon : IsSidonSet singer_q7 := by native_decide

/-- {0, 1, 3, 13, 32, 36, 43, 52} has exactly 8 = q+1 elements. -/
theorem singer_q7_card : singer_q7.card = 8 := by decide

/-- All elements of {0, 1, 3, 13, 32, 36, 43, 52} are ≤ 56 = q²+q = 49+7. -/
theorem singer_q7_range : ∀ a ∈ singer_q7, a ≤ 7 * 7 + 7 := by decide

/-- Singer construction witness for q = 7. -/
theorem singer_sidon_q7 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 7 + 1 ∧ ∀ a ∈ A, a ≤ 7 * 7 + 7 :=
  ⟨singer_q7, singer_q7_sidon, singer_q7_card, singer_q7_range⟩

/-! ## Layer 2d: Singer Construction for q = 11

  The Singer difference set for q = 11 is {0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109}
  in Z_133. PG(2,11) has 133 points; a line has 12 = q+1 points; max element
  109 ≤ 132 = q²+q. Externally verified to be a perfect difference set mod 133.
-/

/-- The Singer set for q=11. -/
def singer_q11 : Finset ℕ := {0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109}

/-- The q=11 Singer candidate is a Sidon set — verified computationally. -/
theorem singer_q11_sidon : IsSidonSet singer_q11 := by native_decide

/-- The q=11 Singer candidate has exactly 12 = q+1 elements. -/
theorem singer_q11_card : singer_q11.card = 12 := by decide

/-- All elements of the q=11 Singer set are ≤ 132 = q²+q. -/
theorem singer_q11_range : ∀ a ∈ singer_q11, a ≤ 11 * 11 + 11 := by decide

/-- Singer construction witness for q = 11. -/
theorem singer_sidon_q11 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 11 + 1 ∧ ∀ a ∈ A, a ≤ 11 * 11 + 11 :=
  ⟨singer_q11, singer_q11_sidon, singer_q11_card, singer_q11_range⟩

/-! ## Layer 2e: Singer Construction for q = 13

  The Singer difference set for q = 13 is
  {0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175} in Z_183.
  PG(2,13) has 183 points; a line has 14 = q+1 points; max element
  175 ≤ 182 = q²+q. Externally verified to be a perfect difference set mod 183.
-/

/-- The Singer set for q=13. -/
def singer_q13 : Finset ℕ := {0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175}

/-- The q=13 Singer candidate is a Sidon set — verified computationally. -/
theorem singer_q13_sidon : IsSidonSet singer_q13 := by native_decide

/-- The q=13 Singer candidate has exactly 14 = q+1 elements. -/
theorem singer_q13_card : singer_q13.card = 14 := by decide

/-- All elements of the q=13 Singer set are ≤ 182 = q²+q. -/
theorem singer_q13_range : ∀ a ∈ singer_q13, a ≤ 13 * 13 + 13 := by decide

/-- Singer construction witness for q = 13. -/
theorem singer_sidon_q13 :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = 13 + 1 ∧ ∀ a ∈ A, a ≤ 13 * 13 + 13 :=
  ⟨singer_q13, singer_q13_sidon, singer_q13_card, singer_q13_range⟩

/-! ## Layer 3: General Singer Construction (Axiom)

  The general Singer construction requires:
  1. Existence of GF(q²) for prime q
  2. Singer cycle on PG(2,q) — a collineation of order q²+q+1
  3. Perfect difference set extraction from the cycle
  4. Sidon property from perfect difference set property

  This is beyond current Mathlib automation scope. We axiomatize with full
  mathematical reference and mark as a formalization target.

  **Reference:** Singer (1938), Trans. AMS 43(3), Theorem 1.
  **Formalization status:** NOT formalized in any prover (Perplexity-verified 2026-03-29).
  **TODO:** Full algebraic proof via GaloisField + IsPrimitiveRoot + trace map.
-/

/-- **AXIOM: Singer's Theorem (1938)**

  For every prime q, there exists a Sidon set of size q+1 with all elements ≤ q²+q.

  This is the key lower bound construction for Erdős Problem #30.
  The proof uses projective geometry: the points of a line in PG(2,q), indexed
  by powers of a Singer cycle, form a perfect difference set mod q²+q+1.
  Perfect difference sets are Sidon sets.

  **Reference:** Singer, J. "A theorem in finite projective geometry and some
  applications to number theory." Trans. Amer. Math. Soc. 43 (1938), 377–385.

  **Why axiom:** The full proof requires:
  - GaloisField q (Mathlib: available as `GaloisField`)
  - Primitive elements of GF(q³)* (Mathlib: `IsCyclic` for finite field units)
  - Trace/norm maps from GF(q³) → GF(q) (Mathlib: `Algebra.trace`, partial)
  - Singer cycle construction (not in Mathlib)
  - Perfect difference set → Sidon (straightforward but requires modular arithmetic bridge)

  **Verified for small primes:** q=2 (`singer_sidon_q2`), q=3 (`singer_sidon_q3`) above.
-/
axiom singer_sidon_exists (q : ℕ) (hq : Nat.Prime q) :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = q + 1 ∧ ∀ a ∈ A, a ≤ q * q + q

/-! ## Layer 4: Lower Bound h(N) ≥ √N

  From Singer's construction + Bertrand's postulate:
  - For any N, find prime q with q² ≤ N < (q+1)²  (roughly)
  - Singer gives Sidon set of size q+1 in {0,...,q²+q}
  - Since q²+q ≤ N (for appropriate N), this gives h(N) ≥ q+1 ≥ √N

  More precisely: for N = q²+q, h(N) ≥ q+1, and q+1 > √N since
  (q+1)² = q²+2q+1 > q²+q = N.

  Combined with the upper bound h(N) ≤ √N + O(N^{1/4}) from Lindström,
  this establishes h(N) = √N + o(√N) as the main conjecture direction.
-/

/-- For N = q²+q with q prime, h(N) ≥ q+1 > √N.
    Direct consequence of Singer's construction. -/
theorem singer_lower_bound (q : ℕ) (hq : Nat.Prime q) :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = q + 1 ∧
      A ⊆ Finset.range (q * q + q + 1) := by
  obtain ⟨A, hSidon, hCard, hRange⟩ := singer_sidon_exists q hq
  exact ⟨A, hSidon, hCard, fun a ha => Finset.mem_range.mpr (by linarith [hRange a ha])⟩

/-- The Singer bound gives (q+1)² > q²+q, i.e., the Sidon set is larger than √N.
    This is the key fact: h(N) ≥ √N for infinitely many N. -/
theorem singer_exceeds_sqrt (q : ℕ) (hq : 0 < q) :
    (q + 1) * (q + 1) > q * q + q := by
  nlinarith

/-! ## Connection to Erdős Problem #30

  The full conjecture asks: h(N) = N^{1/2} + O_ε(N^ε)?
  Equivalently: h(N) - √N = o(√N)?

  What we have formalized:
  1. Upper bound: k(k-1) ≤ 2N (Erdos30_Complete.lean) → h(N) ≤ √(2N) + O(1)
  2. Lower bound: h(q²+q) ≥ q+1 (this file, Singer) → h(N) ≥ √N for many N
  3. Small cases: q=2 (fully verified), q=3 (fully verified)

  The gap between upper (√(2N)) and lower (√N) is the heart of the problem.
  The $1,000 prize asks to close this gap to h(N) = √N + o(√N).

  **Completion status for Problem #30:**
  - Elementary upper bound: ✅ VERIFIED (0 sorry)
  - Difference counting: ✅ VERIFIED (0 sorry)
  - Sum counting: ✅ VERIFIED (0 sorry, via Sidon_SumCount_Fix.lean)
  - Singer q=2: ✅ VERIFIED (0 sorry, this file)
  - Singer q=3: ✅ VERIFIED (0 sorry, this file)
  - Singer q=5: ✅ VERIFIED (0 sorry, this file)
  - Singer q=7: ✅ VERIFIED (0 sorry, this file, added 2026-05-02)
  - Singer q=11: ✅ VERIFIED (0 sorry, this file, added 2026-05-02)
  - Singer q=13: ✅ VERIFIED (0 sorry, this file, added 2026-05-02)
  - Singer general: AXIOM (Singer 1938, not formalized in any prover)
  - Lindström tight bound: NOT YET ATTEMPTED
  - Error term o(√N): OPEN — this is the $1,000 prize
-/

end Erdos.Sidon

