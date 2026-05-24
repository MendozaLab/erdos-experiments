# EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 81 through n = 81 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 9000 |
| Ground-face distance-edge cap | 0 |

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
| 81 | 11 | 1440880837 | 8214 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 81 | 11 | 8214 | 9.013595 | 0 | 0 | NA | 8214 | 0 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 81 | `[0, 1, 9, 20, 35, 51, 56, 68, 74, 78, 81]` | `[0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81]` | `[0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81]` | 0.000000 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 81 | EXPORTED_ALL | 8214 | 8214 | SKIP_DISTANCE_EDGES_TOO_LARGE | 33730791 | true | `[795, 995, 1714, 1952, 2609, 2611, 2615, 3001, 3029, 3052, 3060, 3074, 3563, 3609, 3637, 3732, 3795, 3952, 3964, 4168, 4716, 5005, 5419, 5523, 5914, 6014, 6082, 6514, 6713, 6722, 7133, 7372, 7498, 7561, 7835]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
