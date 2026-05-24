# EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 59 through n = 64 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 10000 |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Reachability Pruning

Enabled with deficiency `0`. After each site, states are retained only if their current cardinality plus remaining sites can still reach `h(n)-0` using the reference `h(n)` for that row.

| n | min cardinality | pruned states | terminal retained states |
|---|---:|---:|---:|
| 59 | 10 | 33858072 | 18 |
| 60 | 10 | 41110095 | 54 |
| 61 | 10 | 49824011 | 152 |
| 62 | 10 | 60275116 | 398 |
| 63 | 10 | 72805828 | 1022 |
| 64 | 10 | 87775523 | 2360 |

## Parity Gate

Checked 6 n-values. h(n) matched in 6. Maximizer counts matched in 6. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 59 | 10 | 18 | 2.890372 | 0 | 0 | NA | 18 | 0 | 0.016531 |
| 60 | 10 | 54 | 3.988984 | 0 | 0 | NA | 54 | 0 | 0.000000 |
| 61 | 10 | 152 | 5.023881 | 0 | 0 | NA | 152 | 0 | 0.001754 |
| 62 | 10 | 398 | 5.986452 | 0 | 0 | NA | 398 | 0 | 0.000000 |
| 63 | 10 | 1022 | 6.929517 | 0 | 0 | NA | 1022 | 0 | 0.001678 |
| 64 | 10 | 2360 | 7.766417 | 0 | 0 | NA | 2360 | 0 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 59 | `[0, 3, 7, 19, 36, 37, 46, 51, 57, 59]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | `[0, 5, 14, 16, 29, 37, 49, 55, 56, 59]` | 0.016531 |
| 60 | `[0, 1, 9, 14, 24, 35, 41, 53, 57, 60]` | `[1, 6, 15, 17, 30, 38, 50, 56, 57, 60]` | `[1, 6, 15, 17, 30, 38, 50, 56, 57, 60]` | 0.000000 |
| 61 | `[0, 1, 7, 17, 29, 38, 42, 53, 56, 61]` | `[0, 4, 9, 24, 31, 43, 49, 57, 59, 60]` | `[0, 5, 12, 25, 28, 46, 47, 55, 57, 61]` | 0.001754 |
| 62 | `[0, 1, 7, 16, 27, 45, 48, 50, 58, 62]` | `[0, 1, 8, 26, 38, 42, 48, 57, 59, 62]` | `[0, 1, 8, 26, 38, 42, 48, 57, 59, 62]` | 0.000000 |
| 63 | `[0, 1, 8, 20, 26, 36, 49, 58, 60, 63]` | `[0, 1, 9, 27, 39, 42, 52, 56, 58, 63]` | `[0, 1, 9, 27, 39, 42, 52, 56, 58, 63]` | 0.001678 |
| 64 | `[0, 1, 8, 19, 34, 39, 51, 55, 61, 64]` | `[0, 1, 8, 28, 39, 45, 49, 58, 61, 63]` | `[0, 1, 15, 33, 36, 40, 45, 56, 62, 64]` | 0.000000 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | field split | pareto minima |
|---|---|---:|---:|---|---|
| 59 | EXPORTED_ALL | 18 | 18 | true | `[7, 13]` |
| 60 | EXPORTED_ALL | 54 | 54 | true | `[43]` |
| 61 | EXPORTED_ALL | 152 | 152 | true | `[85, 96, 113, 139]` |
| 62 | EXPORTED_ALL | 398 | 398 | true | `[50, 74, 202, 228, 363, 390]` |
| 63 | EXPORTED_ALL | 1022 | 1022 | true | `[138, 139, 371, 415, 477, 528, 564, 696, 889, 901]` |
| 64 | EXPORTED_ALL | 2360 | 2360 | true | `[375, 878, 1165, 1273, 1320, 1800, 1861, 1919, 2030, 2166, 2315]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01_RESULTS.sha256 | Integrity checksum |
