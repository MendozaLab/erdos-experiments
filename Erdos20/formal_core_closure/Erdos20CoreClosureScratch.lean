import Mathlib

/-!
# Erdős #20 core-closure scratch

This file is a local structural scratch artifact only.  It encodes the
bookkeeping vocabulary for sunflower core/petal/core-closure work and proves
one lemma: a common pairwise-intersection core makes petals disjoint away from
that core.

It does not prove, improve, or advance the sunflower conjecture.
-/

namespace Erdos20FormalCoreClosure

open Finset

variable {α : Type*} [DecidableEq α]

/-- The petal of `A` away from a proposed core `C`. -/
def Petal (C A : Finset α) : Finset α :=
  A \ C

/-- `C` is the common pairwise intersection for every distinct pair in `F`. -/
def HasCore (C : Finset α) (F : Finset (Finset α)) : Prop :=
  ∀ A ∈ F, ∀ B ∈ F, A ≠ B → A ∩ B = C

/-- The petals of distinct family members are pairwise disjoint away from `C`. -/
def DisjointAwayFromCore (C : Finset α) (F : Finset (Finset α)) : Prop :=
  ∀ A ∈ F, ∀ B ∈ F, A ≠ B → Disjoint (Petal C A) (Petal C B)

/--
Lean-shaped carrier for the finite state needed by a core-closure channel.
`safeCandidates` is intended to be the subset of `candidateExtensions` that
does not close a forbidden sunflower triple through `core`.
-/
structure CoreClosureState where
  family : Finset (Finset α)
  core : Finset α
  candidateExtensions : Finset (Finset α)
  safeCandidates : Finset (Finset α)
  safe_subset_candidates : safeCandidates ⊆ candidateExtensions

/--
A finite bookkeeping cost: how many candidate extensions are ruled unsafe by
the current core-closure state.  This is not the analytic `-log₂ p_safe`; it is
the smallest Nat-valued structural proxy that compiles without importing
analysis/probability machinery.
-/
def CoreClosureCost (S : CoreClosureState (α := α)) : ℕ :=
  S.candidateExtensions.card - S.safeCandidates.card

/--
Smallest genuine bookkeeping lemma: if all distinct pairs in `F` have common
intersection `C`, then their petals are disjoint after removing `C`.
-/
lemma hasCore_disjointAwayFromCore {C : Finset α} {F : Finset (Finset α)}
    (hCore : HasCore C F) :
    DisjointAwayFromCore C F := by
  intro A hA B hB hAB
  rw [Finset.disjoint_iff_inter_eq_empty]
  apply Finset.eq_empty_iff_forall_notMem.mpr
  intro x hx
  have hxPetalA : x ∈ Petal C A := (Finset.mem_inter.mp hx).1
  have hxPetalB : x ∈ Petal C B := (Finset.mem_inter.mp hx).2
  have hxA : x ∈ A := (Finset.mem_sdiff.mp hxPetalA).1
  have hxNotC : x ∉ C := (Finset.mem_sdiff.mp hxPetalA).2
  have hxB : x ∈ B := (Finset.mem_sdiff.mp hxPetalB).1
  have hxInter : x ∈ A ∩ B := Finset.mem_inter.mpr ⟨hxA, hxB⟩
  have hxC : x ∈ C := by
    simpa [hCore A hA B hB hAB] using hxInter
  exact hxNotC hxC

/--
A `B` is a candidate extension of core `C` against family `F` if it has the
uniform width `w`, contains the core, and is not already in the family. This is
the structural shape only; it has no semantics about whether `B` would close a
sunflower triple.
-/
def CandidateExtension (w : ℕ) (C : Finset α) (F : Finset (Finset α))
    (B : Finset α) : Prop :=
  B.card = w ∧ C ⊆ B ∧ B ∉ F

/--
`B` closes a 3-sunflower through core `C` against family `F` if there exist
two distinct family members whose pairwise core is `C`, whose intersections
with `B` are also `C`, and whose petals (away from `C`) are pairwise disjoint
with `B`'s petal. Petal pairwise disjointness follows from
`hasCore_disjointAwayFromCore` reasoning: when `A1 ∩ A2 = C`, `A1 ∩ B = C`, and
`A2 ∩ B = C`, the three petals `A1 \ C`, `A2 \ C`, `B \ C` are pairwise
disjoint.
-/
def ClosesThreeSunflowerThroughCore (C : Finset α) (F : Finset (Finset α))
    (B : Finset α) : Prop :=
  ∃ A1 ∈ F, ∃ A2 ∈ F,
    A1 ≠ A2 ∧
    C ⊆ A1 ∧ C ⊆ A2 ∧ C ⊆ B ∧
    A1 ∩ A2 = C ∧
    A1 ∩ B = C ∧
    A2 ∩ B = C

/--
A safe extension is a candidate extension that does not close a 3-sunflower
through the core.
-/
def SafeExtension (w : ℕ) (C : Finset α) (F : Finset (Finset α))
    (B : Finset α) : Prop :=
  CandidateExtension w C F B ∧ ¬ ClosesThreeSunflowerThroughCore C F B

section SafeFilterTheorems

open scoped Classical

/--
First small bookkeeping theorem: if `safeCandidates` is the filter of
`candidateExtensions` by the `SafeExtension` predicate, then it is a subset of
`candidateExtensions`. Decidability of the filtering predicate is supplied via
`open scoped Classical`.
-/
theorem safeCandidates_eq_filter_and_subset
    (w : ℕ)
    (C : Finset α)
    (F : Finset (Finset α))
    (candidateExtensions : Finset (Finset α))
    (safeCandidates : Finset (Finset α))
    (hSafe :
      safeCandidates =
        candidateExtensions.filter (fun B => SafeExtension w C F B)) :
    safeCandidates ⊆ candidateExtensions := by
  intro B hB
  rw [hSafe] at hB
  exact (Finset.mem_filter.mp hB).1

/--
Second small bookkeeping theorem: assuming the filter equality on
`S.safeCandidates`, the structural cost is exactly the cardinality difference.
The proof is `rfl` because `CoreClosureCost` is already defined as that
difference; the filter hypothesis is recorded for downstream callers but not
needed by this proof.
-/
theorem coreClosureCost_from_filtered_safeCandidates
    (w : ℕ)
    (S : CoreClosureState (α := α))
    (_hSafe :
      S.safeCandidates =
        S.candidateExtensions.filter
          (fun B => SafeExtension w S.core S.family B)) :
    CoreClosureCost S =
      S.candidateExtensions.card - S.safeCandidates.card := by
  rfl

end SafeFilterTheorems

/-!
## W3N5M4 literal-sample bridge

This section encodes one concrete instance from Track A's
`EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01_RESULTS.json` as Lean
literals over `ℕ`, and discharges the filter-equality theorems on the
two qualitative core buckets:

* `C12 = {1,2}` lies fully inside `{1,2,3,4}`; its single candidate
  `{1,2,5}` is closed by the witness pair `A1 = {1,2,3}`, `A2 = {1,2,4}`,
  so `safe_C12 = ∅`.
* `C15 = {1,5}` includes the new element `5`; no member of `F_w3n5m4`
  contains `5`, so no 3-sunflower can be witnessed and all three
  candidates are safe.

The proofs use `decide` against decidability instances on
`CandidateExtension`, `ClosesThreeSunflowerThroughCore`, and
`SafeExtension`. This is a literal-sample bridge between Track A's
W3N5M4 enumerator output and the `SafeExtension` carrier in Lean; it
does not advance the sunflower conjecture.
-/

namespace W3N5M4Sample

open Finset

/-- Family `F` from W3N5M4: the four 3-element subsets of `{1,2,3,4}`. -/
def F_w3n5m4 : Finset (Finset ℕ) :=
  ({{1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}} : Finset (Finset ℕ))

/-- Core that lies entirely inside `{1,2,3,4}`. -/
def C12 : Finset ℕ := ({1,2} : Finset ℕ)

/-- Core that uses the new universe element `5`. -/
def C15 : Finset ℕ := ({1,5} : Finset ℕ)

/-- Track A's candidate set for `C12` at width 3: just `{1,2,5}`. -/
def cand_C12 : Finset (Finset ℕ) :=
  ({{1,2,5}} : Finset (Finset ℕ))

/-- Track A's candidate set for `C15` at width 3. -/
def cand_C15 : Finset (Finset ℕ) :=
  ({{1,2,5}, {1,3,5}, {1,4,5}} : Finset (Finset ℕ))

/-- Track A's safe set for `C12`: empty (the single candidate is closed). -/
def safe_C12 : Finset (Finset ℕ) := (∅ : Finset (Finset ℕ))

/-- Track A's safe set for `C15`: all three candidates. -/
def safe_C15 : Finset (Finset ℕ) :=
  ({{1,2,5}, {1,3,5}, {1,4,5}} : Finset (Finset ℕ))

-- Note: we deliberately do NOT register `Decidable` instances for
-- `CandidateExtension`, `ClosesThreeSunflowerThroughCore`, or
-- `SafeExtension` here. The bridge theorems below run under
-- `open scoped Classical` so that the `Finset.filter` predicate
-- elaborates with the same `Classical.dec` instance used by
-- `safeCandidates_eq_filter_and_subset` and
-- `coreClosureCost_from_filtered_safeCandidates`.

section ClassicalFilter

open scoped Classical

/-- Filter equality for the `C12` bucket: every candidate is closed.

The filter elaborates with the `Classical.dec` instance because
`safeCandidates_eq_filter_and_subset` was stated under
`open scoped Classical`. We discharge the equality at the membership
level via `Finset.ext` so we never need the two decidability instances
to be definitionally identical. -/
theorem hSafe_C12 :
    safe_C12 =
      cand_C12.filter (fun B => SafeExtension 3 C12 F_w3n5m4 B) := by
  apply Finset.ext
  intro B
  constructor
  · intro hB
    -- safe_C12 = ∅, so this branch is vacuous
    simp [safe_C12] at hB
  · intro hB
    rcases Finset.mem_filter.mp hB with ⟨hCand, hSafe⟩
    -- hCand : B ∈ cand_C12, so B = {1,2,5}; hSafe : SafeExtension ...
    -- We will derive a contradiction from hSafe by exhibiting the witness pair.
    exfalso
    have hB_eq : B = ({1,2,5} : Finset ℕ) := by
      simp [cand_C12] at hCand
      exact hCand
    -- hSafe.2 says no 3-sunflower closure through C12 exists, but it does:
    apply hSafe.2
    refine ⟨({1,2,3} : Finset ℕ), ?_, ({1,2,4} : Finset ℕ), ?_,
            ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
    · simp [F_w3n5m4]
    · simp [F_w3n5m4]
    · decide
    · decide
    · decide
    · subst hB_eq; decide
    · decide
    · subst hB_eq; decide
    · subst hB_eq; decide

/-- Filter equality for the `C15` bucket: every candidate is safe.

No element of `F_w3n5m4` contains `5`, so `A1 ∩ B = C15` forces
`5 ∈ A1`, contradiction. We discharge each candidate by `decide` after
substitution. -/
theorem hSafe_C15 :
    safe_C15 =
      cand_C15.filter (fun B => SafeExtension 3 C15 F_w3n5m4 B) := by
  apply Finset.ext
  intro B
  constructor
  · intro hB
    rw [Finset.mem_filter]
    simp [safe_C15] at hB
    refine ⟨?_, ?_, ?_⟩
    · -- B ∈ cand_C15
      simp [cand_C15]
      exact hB
    · -- CandidateExtension 3 C15 F_w3n5m4 B
      rcases hB with h | h | h
      all_goals (subst h; refine ⟨?_, ?_, ?_⟩ <;> decide)
    · -- ¬ ClosesThreeSunflowerThroughCore C15 F_w3n5m4 B
      -- Any A1 ∈ F_w3n5m4 has A1 ⊆ {1,2,3,4}, so 5 ∉ A1, but 5 ∈ C15 ⊆ A1
      -- would force 5 ∈ A1 — contradiction.
      rintro ⟨A1, hA1F, A2, hA2F, _, hCA1, _, _, _, _, _⟩
      have h5_in_C : (5 : ℕ) ∈ C15 := by decide
      have h5_in_A1 : (5 : ℕ) ∈ A1 := hCA1 h5_in_C
      have hA1_no5 : (5 : ℕ) ∉ A1 := by
        simp [F_w3n5m4] at hA1F
        rcases hA1F with h | h | h | h <;> (subst h; decide)
      exact hA1_no5 h5_in_A1
  · intro hB
    rcases Finset.mem_filter.mp hB with ⟨hCand, _⟩
    simp [safe_C15]
    simp [cand_C15] at hCand
    exact hCand

/-- Subset corollary for `C12`, via `safeCandidates_eq_filter_and_subset`. -/
theorem safe_C12_subset : safe_C12 ⊆ cand_C12 :=
  safeCandidates_eq_filter_and_subset
    3 C12 F_w3n5m4 cand_C12 safe_C12 hSafe_C12

/-- Subset corollary for `C15`. -/
theorem safe_C15_subset : safe_C15 ⊆ cand_C15 :=
  safeCandidates_eq_filter_and_subset
    3 C15 F_w3n5m4 cand_C15 safe_C15 hSafe_C15

/-- `CoreClosureState` instance for the `C12` bucket. -/
def S_C12 : CoreClosureState (α := ℕ) :=
  { family := F_w3n5m4
  , core := C12
  , candidateExtensions := cand_C12
  , safeCandidates := safe_C12
  , safe_subset_candidates := safe_C12_subset }

/-- `CoreClosureState` instance for the `C15` bucket. -/
def S_C15 : CoreClosureState (α := ℕ) :=
  { family := F_w3n5m4
  , core := C15
  , candidateExtensions := cand_C15
  , safeCandidates := safe_C15
  , safe_subset_candidates := safe_C15_subset }

/-- Cost lemma applied to `S_C12`. -/
theorem cost_C12_eq :
    CoreClosureCost S_C12 =
      S_C12.candidateExtensions.card - S_C12.safeCandidates.card :=
  coreClosureCost_from_filtered_safeCandidates 3 S_C12 hSafe_C12

/-- Cost lemma applied to `S_C15`. -/
theorem cost_C15_eq :
    CoreClosureCost S_C15 =
      S_C15.candidateExtensions.card - S_C15.safeCandidates.card :=
  coreClosureCost_from_filtered_safeCandidates 3 S_C15 hSafe_C15

/-- The `C12` bucket has cost `1 - 0 = 1`: one candidate, all unsafe. -/
example : CoreClosureCost S_C12 = 1 := by
  unfold CoreClosureCost S_C12 cand_C12 safe_C12
  decide

/-- The `C15` bucket has cost `3 - 3 = 0`: all candidates are safe. -/
example : CoreClosureCost S_C15 = 0 := by
  unfold CoreClosureCost S_C15 cand_C15 safe_C15
  decide

end ClassicalFilter

end W3N5M4Sample

end Erdos20FormalCoreClosure
