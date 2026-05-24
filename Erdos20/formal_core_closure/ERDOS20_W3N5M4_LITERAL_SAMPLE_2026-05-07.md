# Erdős #20 — W3N5M4 Literal-Sample Bridge (Worker D, 2026-05-07)

## Claim ceiling

A0. Required preserved phrase: **"shadow signature, not universal law"**.

Honest scope: this artifact is a compiled literal-sample bridge between
Track A's W3N5M4 enumerator output and the `SafeExtension` carrier in
Lean. It does not advance the sunflower conjecture. The encoded
instance is a shadow signature, not universal law.

## What was added

A new `namespace W3N5M4Sample` was appended inside
`Erdos20FormalCoreClosure` in
`erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean`.

New definitions:
- `F_w3n5m4 : Finset (Finset ℕ)` — `{{1,2,3}, {1,2,4}, {1,3,4}, {2,3,4}}`
- `C12 : Finset ℕ := {1,2}` — core inside `{1,2,3,4}`
- `C15 : Finset ℕ := {1,5}` — core touching the new universe element
- `cand_C12 : Finset (Finset ℕ) := {{1,2,5}}`
- `cand_C15 : Finset (Finset ℕ) := {{1,2,5}, {1,3,5}, {1,4,5}}`
- `safe_C12 : Finset (Finset ℕ) := ∅`
- `safe_C15 : Finset (Finset ℕ) := {{1,2,5}, {1,3,5}, {1,4,5}}`
- `S_C12`, `S_C15 : CoreClosureState (α := ℕ)` — full state instances

New theorems / examples:
- `hSafe_C12` and `hSafe_C15` — filter-equality theorems
- `safe_C12_subset`, `safe_C15_subset` — subset corollaries via
  `safeCandidates_eq_filter_and_subset`
- `cost_C12_eq`, `cost_C15_eq` — cost lemmas via
  `coreClosureCost_from_filtered_safeCandidates`
- Two `example` lemmas computing the cost values:
  - `CoreClosureCost S_C12 = 1`
  - `CoreClosureCost S_C15 = 0`

## Tactic used to discharge the filter equalities

Neither `decide` nor `native_decide` could close `hSafe_C12` /
`hSafe_C15` directly because the existing
`safeCandidates_eq_filter_and_subset` was elaborated under
`open scoped Classical`, so `Finset.filter (SafeExtension ...)`
elaborates with `Classical.dec`. Registering local `Decidable`
instances for `SafeExtension` produced a competing instance that
caused application-type-mismatch errors when feeding `hSafe_*` into
the upstream theorems.

Resolution: discharge the filter equality at the membership level via
`Finset.ext`, splitting on `B`'s membership in the candidate set:

- `hSafe_C12`: forward direction is vacuous (`safe_C12 = ∅`). Backward
  direction extracts `B = {1,2,5}` and exhibits the explicit witness
  pair `A1 = {1,2,3}`, `A2 = {1,2,4}` to contradict
  `¬ ClosesThreeSunflowerThroughCore`. Each membership and intersection
  fact closes by `decide`.

- `hSafe_C15`: backward direction is direct (`safe_C15 = cand_C15`).
  Forward direction shows safety: every `A1 ∈ F_w3n5m4` lies inside
  `{1,2,3,4}`, so `5 ∉ A1`, but `5 ∈ C15 ⊆ A1` would force `5 ∈ A1` —
  contradiction. Each `A1` case closes by `decide` after substitution.

The two `cost_*` examples close by `unfold` + `decide` on pure `Nat`
arithmetic over `Finset.card` of literal sets.

## Computed `CoreClosureCost` values

- `CoreClosureCost S_C12 = 1` (one candidate, zero safe → cost 1)
- `CoreClosureCost S_C15 = 0` (three candidates, all safe → cost 0)

Both verified by `decide` after unfolding.

## Build verification

Command:
```
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/Math-Problems/proofs/archive
lake env lean /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean
```

Exit code: **0** (clean, no warnings).

## sorry / axiom counts

```
$ grep -c "sorry" Erdos20CoreClosureScratch.lean
0
$ grep -c "^axiom" Erdos20CoreClosureScratch.lean
0
```

Token embargo (the 2026-04-24 embargoed token) check on this file: 0 occurrences after fix.

## Files changed

- `erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean`
  — appended `namespace W3N5M4Sample` block (~140 lines).
- `erdos-experiments/Erdos20/formal_core_closure/ERDOS20_W3N5M4_LITERAL_SAMPLE_2026-05-07.md`
  — new (this file).

No files outside `formal_core_closure/` were modified. No git commit
was made. No D1 / scorecard / public-doc updates.
