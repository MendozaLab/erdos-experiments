# EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 56 through n = 59 |
| Frontier k | 5 |
| Prune deficiency | 2 |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Reachability Pruning

Enabled with deficiency `2`. After each site, states are retained only if their current cardinality plus remaining sites can still reach `h(n)-2` using the reference `h(n)` for that row.

| n | min cardinality | pruned states | terminal retained states |
|---|---:|---:|---:|
| 56 | 8 | 22899068 | 4882634 |
| 57 | 8 | 27146114 | 6631722 |
| 58 | 8 | 32095667 | 8885328 |
| 59 | 8 | 37797715 | 11830796 |

## Parity Gate

Checked 4 n-values. h(n) matched in 4. Maximizer counts matched in 4. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 56 | 10 | 4 | 1.386294 | 69564 | 4813066 | 1 | 4882634 | 0 | 0.011840 |
| 57 | 10 | 6 | 1.791759 | 120704 | 6511012 | 1 | 6631722 | 0 | 0.028889 |
| 58 | 10 | 10 | 2.302585 | 200946 | 8684372 | 1 | 8885328 | 0 | 0.036164 |
| 59 | 10 | 18 | 2.890372 | 330056 | 11500722 | 1 | 11830796 | 0 | 0.016531 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 56 | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | 0.011840 |
| 57 | `[2, 3, 8, 12, 25, 28, 36, 43, 55, 57]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | 0.028889 |
| 58 | `[0, 2, 15, 21, 22, 32, 46, 50, 55, 58]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | 0.036164 |
| 59 | `[0, 3, 7, 19, 36, 37, 46, 51, 57, 59]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | `[0, 5, 14, 16, 29, 37, 49, 55, 56, 59]` | 0.016531 |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30_RESULTS.sha256 | Integrity checksum |
