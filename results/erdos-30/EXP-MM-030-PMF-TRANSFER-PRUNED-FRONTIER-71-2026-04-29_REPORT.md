# EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 71 through n = 71 |
| Frontier k | 10 |
| Prune deficiency | 0 |

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
| 71 | 10 | 311297246 | 203840 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 71 | 10 | 203840 | 12.225091 | 0 | 0 | NA | 203840 | 0 | 0.001424 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 71 | `[0, 4, 13, 21, 40, 45, 56, 68, 70, 71]` | `[0, 1, 6, 31, 44, 53, 55, 63, 67, 70]` | `[0, 4, 13, 23, 34, 51, 63, 65, 66, 71]` | 0.001424 |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29_RESULTS.sha256 | Integrity checksum |
