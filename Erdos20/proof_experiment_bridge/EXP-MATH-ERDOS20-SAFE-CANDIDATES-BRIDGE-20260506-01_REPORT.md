# EXP-MATH-ERDOS20-SAFE-CANDIDATES-BRIDGE-20260506-01 Report

## Status

- Status: `SAFE_CANDIDATES_BRIDGE_PASS_A0`
- Target: `w=3, n=7, k=3, core_size_s=2`
- Claim ceiling: A0, shadow signature, not universal law

## What This Binds

This bridge emits candidate-level rows for the Lean-shaped relation:

```text
safeCandidates = candidateExtensions.filter (SafeExtension core family)
```

Each unsafe candidate includes concrete witness pairs from the family when
available. This is still bookkeeping evidence, not theorem progress.

## Ratios

The run reports both controls requested after third-party review:

- `floor_ratio = I_core_local / I_floor`, with `I_floor = 1 bit`
- `ahs_ratio = I_core_local / log2(sqrt(10))`, using the Q1 gate's
  Abbott-Hansen-Sauer base `c_3 >= sqrt(10)` as a construction/literature
  control, not an exact local construction model

| m | core | candidates | safe | I_core_local | floor ratio | AHS ratio |
|---:|---|---:|---:|---:|---:|---:|
| 6 | [2, 3] | 4 | 4 | 0.0 | 0.0 | 0.0 |
| 7 | [2, 5] | 5 | 5 | 0.0 | 0.0 | 0.0 |
| 8 | [2, 5] | 4 | 4 | 0.0 | 0.0 | 0.0 |

## Interpretation

This turns the prior aggregate/per-core summaries into proof-facing rows. It
does not reproduce an Abbott-Hansen-Sauer construction, does not prove a
sunflower theorem, and does not upgrade the A-axis.
