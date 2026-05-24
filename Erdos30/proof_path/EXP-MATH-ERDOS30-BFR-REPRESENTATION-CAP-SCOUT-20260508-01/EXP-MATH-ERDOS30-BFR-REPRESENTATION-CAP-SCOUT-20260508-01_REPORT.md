# EXP-MATH-ERDOS30-BFR-REPRESENTATION-CAP-SCOUT-20260508-01

## Meaning

Eratosthenes advanced the #30 BFR bridge by isolating the representation-function mechanics that must be true before any constant-bearing inequality is worth touching. In Hilbert terms, this is still the grammar of the operator: count the fibers correctly, separate the swap sector, and only then ask where the `0.98183` coefficient can be pressured.

The new Lean scout file now records the central cap:

```lean
theorem bfrRepFunction_sidon_le_two
    (A : Finset Nat) (s : Nat) (hS : IsSidonSet A) :
    bfrRepFunction A s <= 2
```

This does not improve Erdős #30. It says the local carrier and representation count are now disciplined enough to review as a bridge block.

## Verification

- Command: `lake env lean scratch/Erdos30_BFR_RepresentationBridge_SCOUT.lean`
- Result: exit code `0`
- Escape-hatch scan: `rg -n "sorry|axiom|admit" scratch/Erdos30_BFR_RepresentationBridge_SCOUT.lean`
- Matches: `0`
- Scout file SHA-256: `09e4f164329ab992c303cd17ee2acae2d327e8742fb35066c5229571875bfdcd`

## What Was Added

- `bfrRepFunction_le_half_card_le_one`: the `a <= b` half of the ordered representation fiber has cardinality at most one under the Sidon condition.
- `bfrRepFunction_gt_half_card_le_one`: the `b < a` half has cardinality at most one after reversing witnesses.
- `bfrRepFunction_sidon_le_two`: the two halves combine into the expected ordered cap.
- `bfr_sum_repFunction_eq_card_sq`: the support-range sum of the representation function is the ordered square cardinality.

## Boundary

This is a review-only scaffold. It does not change the public coefficient, update any registry, update D1, touch curated morphisms, or alter existing #30/#114 result artifacts.

## Next Move

Review the theorem names and local carrier choices against `lean/Erdos30_BFR.lean`. If the style is acceptable, the next phase is a narrow port of this bridge into the canonical BFR file, still without claiming any coefficient movement.
