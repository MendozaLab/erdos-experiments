# EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 73 through n = 73 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 100 |
| Ground-face distance-edge cap | 100 |

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
| 73 | 11 | 371124979 | 8 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 73 | 11 | 8 | 2.079442 | 0 | 0 | NA | 8 | 0 | 0.057564 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 73 | `[1, 2, 5, 14, 29, 34, 48, 55, 65, 71, 73]` | `[1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73]` | `[1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73]` | 0.057564 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 73 | EXPORTED_ALL | 8 | 8 | EXPORTED_ALL | 28 | true | `[6]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_RESULTS.sha256 | Integrity checksum |
