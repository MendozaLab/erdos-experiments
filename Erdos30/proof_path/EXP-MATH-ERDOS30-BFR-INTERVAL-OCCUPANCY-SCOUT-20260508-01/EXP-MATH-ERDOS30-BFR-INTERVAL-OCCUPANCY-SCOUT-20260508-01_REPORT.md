# EXP-MATH-ERDOS30-BFR-INTERVAL-OCCUPANCY-SCOUT-20260508-01

## Meaning

Eratosthenes used the newly landed representation-function bridge to begin the next BFR block: interval occupancy. The scout defines a clipped interval inside the natural support range, sums the representation function over that interval, and checks that one interval cannot exceed the total ordered representation mass.

This still does not touch the coefficient. It prepares the bookkeeping needed for the BFR discrepancy argument.

## Verification

- Module build: `lake build Erdos30_BFR`
- Module build result: exit code `0`
- Scout check: `lake env lean scratch/Erdos30_BFR_IntervalOccupancy_SCOUT.lean`
- Scout check result: exit code `0`
- Escape-hatch scan: `rg -n "sorry|axiom|admit" scratch/Erdos30_BFR_IntervalOccupancy_SCOUT.lean`
- Matches: `0`
- Scout SHA-256: `75a23b2e37833197d165692a2174988caaa2ee9d4f1e401fdd29dfd1ca845f6b`

## What Was Added

- `bfrIntervalInRange`: a width-`v` interval clipped to `bfrSumRange N`.
- `mem_bfrIntervalInRange_iff`: membership unpacking for that interval.
- `bfrIntervalOccupancy`: interval mass using `bfrRepFunction`.
- `bfrIntervalInRange_subset_sumRange`: the clipped interval stays inside the BFR range.
- `bfrIntervalOccupancy_le_total`: one interval is bounded by the total mass.
- `bfrIntervalOccupancy_le_card_sq`: one interval is bounded by `|A|^2`.

## Boundary

This is a review-only scout. It does not alter D1, curated morphisms, public pages, or #114 artifacts. It does not change the public #30 bound.

## Next Move

Build the actual finite interval partition: coverage, pairwise disjointness, and a theorem saying the sum of interval occupancies over the partition recovers the support-range mass.
