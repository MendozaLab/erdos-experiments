# EXP-MM-030-PMF-GROUNDFACE-RESET-76-PAR8-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-RESET-76-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 76 through n = 76 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 300 |
| Ground-face distance-edge cap | 30000 |

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
| 76 | 11 | 622974564 | 214 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 76 | 11 | 214 | 5.365976 | 0 | 0 | NA | 214 | 0 | 0.002593 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 76 | `[0, 1, 8, 21, 33, 36, 47, 52, 70, 74, 76]` | `[0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76]` | `[0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76]` | 0.002593 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 76 | EXPORTED_ALL | 214 | 214 | EXPORTED_ALL | 22791 | true | `[40]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-RESET-76-PAR8-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-RESET-76-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-RESET-76-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
