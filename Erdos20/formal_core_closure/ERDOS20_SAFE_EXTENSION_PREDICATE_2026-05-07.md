# Erdős #20 Safe Extension Predicate — Lean Bridge Step 2

**Date:** 2026-05-07
**Worker:** #20-B (continuation)
**Lane:** sunflower core-closure formalization (proof architecture only)
**Scope:** local artifact only — no D1, no scorecard, no git push, no public surface.
**Required phrase:** shadow signature, not universal law.

## Bottom Line

Extended `Erdos20CoreClosureScratch.lean` with the `SafeExtension` predicate
family and two filtered-subset bookkeeping theorems specified by the May 6
bridge packet. The file compiles cleanly under the current Mathlib pin
(`f897ebcf72cd16f89ab4577d0c826cd14afaafc7`) using `lake env lean` from
`Math-Problems/proofs/archive`. Honest scope statement: this is **proof
architecture for the SafeExtension carrier; it does not advance the sunflower
conjecture**. Status remains A0: shadow signature, not universal law.

## What Was Added

Four new pieces, all in real Lean 4 + Mathlib syntax (`⊆`, `∉`, `∩`, `∃ ∈`):

1. `def CandidateExtension (w : ℕ) (C : Finset α) (F : Finset (Finset α)) (B : Finset α) : Prop`
   — `B.card = w ∧ C ⊆ B ∧ B ∉ F`.
2. `def ClosesThreeSunflowerThroughCore (C : Finset α) (F : Finset (Finset α)) (B : Finset α) : Prop`
   — existence of distinct `A1, A2 ∈ F` with `A1 ∩ A2 = C`, `A1 ∩ B = C`,
     `A2 ∩ B = C`, and `C ⊆ A1, A2, B`. Petal pairwise disjointness follows
     from the existing `hasCore_disjointAwayFromCore` reasoning, so it is not
     re-encoded here.
3. `def SafeExtension (w : ℕ) (C : Finset α) (F : Finset (Finset α)) (B : Finset α) : Prop`
   — `CandidateExtension ∧ ¬ ClosesThreeSunflowerThroughCore`.
4. Two theorems inside a `section ... open scoped Classical ... end`:
   - `safeCandidates_eq_filter_and_subset` — three-line proof: `intro B hB; rw [hSafe] at hB; exact (Finset.mem_filter.mp hB).1`.
   - `coreClosureCost_from_filtered_safeCandidates` — `rfl`, since
     `CoreClosureCost` is already defined as the cardinality difference.

## Decidability Note

`SafeExtension` is a `Prop` mixing `∃ ∈` and `≠`, so `Finset.filter` would
need a `DecidablePred`. The cleanest fix on current Mathlib was to wrap both
theorems in:

```lean
section SafeFilterTheorems
open scoped Classical
...
end SafeFilterTheorems
```

This is a structural-bookkeeping-only artifact, so classical decidability for
the filter predicate is acceptable. No `axiom`, no `sorry`.

## Build Command and Output

Command (run from `Math-Problems/proofs/archive`):

```bash
lake env lean /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean
```

Output: clean — no errors, no warnings, no stdout.

Mathlib pin (from `lake-manifest.json`): `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`.

## Counts

| Metric | Value |
|---|---|
| `sorry` count (file) | 0 |
| `axiom` count (file) | 0 |
| Embargoed-token count (file) | 0 |
| Build status | PASS (clean `lake env lean`) |
| Theorems compiled | 2 (plus the prior `hasCore_disjointAwayFromCore`) |
| Definitions added | 3 (`CandidateExtension`, `ClosesThreeSunflowerThroughCore`, `SafeExtension`) |

Verification commands:

```bash
grep -c "sorry" erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean   # 0
grep -c "axiom" erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean   # 0
grep -ci "<embargoed-token>" erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean # 0
```

## Claim Ceiling

A0. The predicate `SafeExtension` is now compiled as a Lean type, and the two
filtered-subset bookkeeping theorems hold against the current scratch carrier
`CoreClosureState`. This does not formalize, encode, compute, measure,
observe, recover, or specify any new mathematical content beyond the bridge
packet's specified targets.

Honest scope statement: **proof architecture for the SafeExtension carrier;
does not advance the sunflower conjecture.**

Required phrase: shadow signature, not universal law.

## Files Changed

- Modified: `erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean`
- Created: `erdos-experiments/Erdos20/formal_core_closure/ERDOS20_SAFE_EXTENSION_PREDICATE_2026-05-07.md`

No other files modified. No commits made. No D1 writes. No scorecard updates.

## Next Concrete Step (Not in This Lane)

Per the bridge packet §"Proof-to-Experiment Bridge Procedure", the next step
would be to instantiate one concrete `family / core / candidateExtensions /
safeCandidates` block from a small enumerator run (e.g., `w=3,n=5,k=3,m=4`)
and feed it through `safeCandidates_eq_filter_and_subset` to verify the
predicate matches the enumerator's `local_safe_through_core` flag on real
literals. That step requires the `--emit-safe-candidates-v1` enumerator
output, which is out of scope for the current worker lane.
