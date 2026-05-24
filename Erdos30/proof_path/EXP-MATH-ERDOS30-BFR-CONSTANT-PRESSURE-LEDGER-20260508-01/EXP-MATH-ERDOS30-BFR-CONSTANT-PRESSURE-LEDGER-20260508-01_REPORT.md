# EXP-MATH-ERDOS30-BFR-CONSTANT-PRESSURE-LEDGER-20260508-01

Status: REVIEW_ONLY / CONSTANT_PRESSURE_LEDGER_READY

This ledger separates coefficient pressure from rounding and finite-channel sanity checks.

## Decision

`NEEDS_EXTERNAL_THEOREM`

The next useful act is to unpack the constant-bearing BFR interface into named inequalities and Lean-sized bridge lemmas. No improvement to `0.98183` is claimed.

## Rows

- `bfr_core_bound`: `NEEDS_EXTERNAL_THEOREM` — constant-bearing step is external/axiom-facing until dependency theorem is fully unpacked
- `lindstrom_rounding`: `ROUNDING_ONLY` — Phase 2 flags Nat.sqrt/rounding as able to fake progress
- `collision_channel_sanity`: `PRESSURE_CONTEXT_ONLY` — R7 recovers known constructions and null separation but does not pressure asymptotic coefficient
