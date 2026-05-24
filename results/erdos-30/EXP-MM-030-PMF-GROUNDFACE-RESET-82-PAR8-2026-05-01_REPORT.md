# EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 82 through n = 82 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 17000 |
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

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 82 | EXPORTED_ALL | 15958 | 15958 | SKIP_DISTANCE_EDGES_TOO_LARGE | 127320903 | true | `[1696, 3600, 4475, 5452, 5687, 6154, 6205, 6247, 6298, 6300, 6328, 6725, 6729, 7268, 7514, 7607, 8538, 8734, 8767, 8771, 9462, 10663, 10782, 10831, 10836, 11100, 11398, 11473, 11497, 11557, 11679, 12492, 13221, 13240, 13463, 13670, 13751, 13764, 13774, 14738, 15070, 15601, 15702, 15812, 15844, 15850]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
