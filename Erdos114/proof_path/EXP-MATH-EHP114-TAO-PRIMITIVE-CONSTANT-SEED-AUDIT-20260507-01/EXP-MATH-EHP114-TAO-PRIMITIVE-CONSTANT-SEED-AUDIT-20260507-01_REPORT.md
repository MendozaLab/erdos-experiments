# EHP114 Tao Primitive Constant Seed Audit

Experiment: `EXP-MATH-EHP114-TAO-PRIMITIVE-CONSTANT-SEED-AUDIT-20260507-01`

## Verdict

Status: `PRIMITIVE_CONSTANT_SEED_AUDIT_PARTIAL`.

This packet separates primitive Tao-side constants into numeric seeds, symbolic seeds, and opaque asymptotic rows. It does not emit a global threshold.

## Summary

- Processed primitives: `12`
- Numeric seed rows: `6`
- Symbolic seed rows: `4`
- Opaque seed rows: `2`
- Source-label missing rows: `0`
- Tao source SHA status: `PASS_PINNED_SOURCE_SHA_MATCH`

## First Failed Condition

`poz and pocl exponential constants remain opaque`

## Claim Ceiling

Primitive constant seed audit only. This packet does not provide a numeric Tao threshold, does not authorize n=15 or higher computation, and does not prove the full all-degree statement.
