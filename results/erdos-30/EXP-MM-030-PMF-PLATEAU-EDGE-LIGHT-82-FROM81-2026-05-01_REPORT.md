# EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-82-FROM81-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-82-FROM81-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 82 through n = 82 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | none |
| Ground-face distance-edge cap | none |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01_RESULTS.json |

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
| 82 | 11 | 1698049081 | 15958 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 82 | 11 | 15958 | 9.677716 | 0 | 0 | NA | 15958 | 0 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 82 | `[0, 1, 10, 21, 34, 37, 59, 63, 65, 77, 82]` | `[0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82]` | `[0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82]` | 0.000000 |

## Plateau Inheritance Probe

This light certificate checks whether the current exact face contains the previous exact face and the previous face shifted by `+1`, without requiring pairwise distance-edge export.

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 82 | 81 | 8214 | 12398 | 8214 | 8214 | 12398 | 3560 | true |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-82-FROM81-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-82-FROM81-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-82-FROM81-2026-05-01_RESULTS.sha256 | Integrity checksum |
