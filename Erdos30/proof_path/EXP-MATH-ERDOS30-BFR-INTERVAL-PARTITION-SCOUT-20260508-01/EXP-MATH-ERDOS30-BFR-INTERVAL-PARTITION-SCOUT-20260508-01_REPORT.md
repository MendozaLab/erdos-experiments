# EXP-MATH-ERDOS30-BFR-INTERVAL-PARTITION-SCOUT-20260508-01

## Meaning

Eratosthenes took the interval-occupancy scaffold one step further: the BFR support range is now covered by a finite family of clipped intervals, and those intervals are pairwise disjoint. This is the floorplan needed before the discrepancy term can be written in Lean.

The strongest new statement is:

```lean
theorem bfrPartitionOccupancy_sum_eq_card_sq
    (A : Finset Nat) (N v : Nat)
    (hA : A ⊆ Finset.range (N + 1)) (hv : 0 < v) :
    ∑ i ∈ bfrPartitionIndexSet N v, bfrPartitionOccupancy A N v i =
      A.card * A.card
```

That says the interval rooms recover the whole ordered representation mass. This is still bookkeeping, not coefficient work.

## Verification

- Command: `lake env lean scratch/Erdos30_BFR_IntervalPartition_SCOUT.lean`
- Result: exit code `0`
- Escape-hatch scan: `rg -n "sorry|axiom|admit" scratch/Erdos30_BFR_IntervalPartition_SCOUT.lean`
- Matches: `0`
- Warnings: none
- Scout SHA-256: `4bea97b6ce68eeae52610c86e266160efcae793cc4d09a57c6edc4de35e946c4`

## What Was Added

- A finite interval family over `bfrSumRange N`.
- A quotient-index membership lemma for any support-range point.
- A cover theorem: the interval family covers the support range for positive width.
- A disjointness theorem: distinct intervals do not overlap.
- A mass theorem: summing interval occupancies recovers the support-range mass.
- A card-square theorem: under `A ⊆ {0,...,N}`, total partition occupancy is `|A|^2`.

## Boundary

This is a review-only scout. It does not alter D1, curated morphisms, public pages, or #114 artifacts. It does not change the public #30 bound.

## Next Move

Choose the exact interval-width convention from the BFR paper, then define the discrepancy variable on this partition and connect it to `bfr_cauchy_schwarz`.
