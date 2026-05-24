# EXP-MATH-ERDOS30-BFR-REPRESENTATION-BRIDGE-SCOUT-20260508-01

Status: `REVIEW_ONLY / REPRESENTATION_BRIDGE_SCOUT_CHECK_PASSED`

## Decision

`A1_SUPPORT_RANGE_BLOCK_STARTED`

This packet executes the first concrete step from the #30 BFR Lean bridge queue. It adds an isolated scratch Lean file for the ordered representation function `r_A(s)` and the natural support range for sums when `A ⊆ {0,...,N}`.

## Lean Check

Command:

```bash
lake env lean scratch/Erdos30_BFR_RepresentationBridge_SCOUT.lean
```

Result: exit code `0`. This is a narrow single-file environment check, not a package-wide gate.

## New Formal Items

- `bfrRepFunction`
- `bfrRepFunction_le_card_product`
- `bfrRepFunction_le_card_sq`
- `bfrRepFunction_sum_mem_range_of_pair`
- `bfrRepFunction_eq_zero_of_not_mem_range`
- `bfrSumRange`
- `mem_bfrSumRange_iff`
- `bfrRepFunction_eq_zero_of_not_mem_bfrSumRange`

## Next Scout Targets

- `bfr_sum_repFunction_eq_card_sq`
- `bfrRepFunction_sidon_le_two`
- `bfr_addEnergy_sidon`

## Claim Ceiling

Review-only local formalization scout. This does not change the public #30 coefficient, does not close `bfr_core_bound`, and does not write to any registry or public surface.
