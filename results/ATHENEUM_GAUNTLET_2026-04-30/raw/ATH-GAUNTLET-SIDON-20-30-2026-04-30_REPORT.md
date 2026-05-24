# ATH-GAUNTLET-SIDON-20-30-2026-04-30 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | ATH-GAUNTLET-SIDON-20-30-2026-04-30 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 20 through n = 30 |
| Frontier k | 3 |
| Prune deficiency | none |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Reachability Pruning

Disabled. All exact transfer states are retained.

| n | min cardinality | pruned states | terminal retained states |
|---|---:|---:|---:|
| 20 | NA | 0 | 9188 |
| 21 | NA | 0 | 12366 |
| 22 | NA | 0 | 16417 |
| 23 | NA | 0 | 21787 |
| 24 | NA | 0 | 28708 |
| 25 | NA | 0 | 37722 |
| 26 | NA | 0 | 49083 |
| 27 | NA | 0 | 63921 |
| 28 | NA | 0 | 82640 |
| 29 | NA | 0 | 106722 |
| 30 | NA | 0 | 136675 |

## Parity Gate

Checked 11 n-values. h(n) matched in 11. Maximizer counts matched in 11. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 6 | 206 | 5.327876 | 3734 | 3786 | 1 | 9188 | 9188 | 0.000000 |
| 21 | 6 | 504 | 6.222576 | 5428 | 4750 | 1 | 12366 | 12366 | 0.007602 |
| 22 | 6 | 1004 | 6.911747 | 7612 | 5874 | 1 | 16417 | 16417 | 0.000000 |
| 23 | 6 | 1910 | 7.554859 | 10488 | 7196 | 1 | 21787 | 21787 | 0.006708 |
| 24 | 6 | 3380 | 8.125631 | 14126 | 8720 | 1 | 28708 | 28708 | 0.000000 |
| 25 | 7 | 10 | 2.302585 | 5688 | 18744 | 1 | 37722 | 37722 | 0.023926 |
| 26 | 7 | 34 | 3.526361 | 9036 | 24390 | 1 | 49083 | 49083 | 0.000000 |
| 27 | 7 | 98 | 4.584967 | 14106 | 31436 | 1 | 63921 | 63921 | 0.000000 |
| 28 | 7 | 282 | 5.641907 | 21190 | 39914 | 1 | 82640 | 82640 | 0.000000 |
| 29 | 7 | 760 | 6.633318 | 31158 | 50212 | 1 | 106722 | 106722 | 0.000000 |
| 30 | 7 | 1618 | 7.388946 | 44370 | 62390 | 1 | 136675 | 136675 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 20 | `[0, 3, 7, 12, 18, 20]` | `[0, 3, 13, 15, 19, 20]` | `[0, 3, 13, 15, 19, 20]` | 0.000000 |
| 21 | `[0, 3, 8, 15, 17, 21]` | `[0, 3, 14, 16, 20, 21]` | `[0, 3, 14, 16, 20, 21]` | 0.007602 |
| 22 | `[0, 4, 9, 19, 20, 22]` | `[0, 4, 14, 16, 21, 22]` | `[0, 4, 14, 16, 21, 22]` | 0.000000 |
| 23 | `[0, 4, 9, 15, 22, 23]` | `[0, 2, 14, 19, 22, 23]` | `[0, 4, 13, 20, 21, 23]` | 0.006708 |
| 24 | `[0, 5, 11, 15, 23, 24]` | `[0, 2, 15, 20, 23, 24]` | `[0, 5, 11, 21, 23, 24]` | 0.000000 |
| 25 | `[0, 1, 7, 11, 20, 23, 25]` | `[0, 4, 9, 15, 22, 23, 25]` | `[0, 4, 9, 15, 22, 23, 25]` | 0.023926 |
| 26 | `[0, 1, 6, 14, 17, 24, 26]` | `[0, 2, 12, 18, 21, 25, 26]` | `[0, 2, 12, 18, 21, 25, 26]` | 0.000000 |
| 27 | `[0, 2, 9, 15, 23, 26, 27]` | `[0, 2, 12, 18, 23, 26, 27]` | `[0, 2, 12, 18, 23, 26, 27]` | 0.000000 |
| 28 | `[0, 2, 7, 16, 24, 27, 28]` | `[0, 4, 15, 18, 20, 27, 28]` | `[0, 4, 15, 18, 20, 27, 28]` | 0.000000 |
| 29 | `[0, 3, 8, 15, 19, 28, 29]` | `[0, 3, 12, 22, 23, 27, 29]` | `[0, 3, 12, 22, 23, 27, 29]` | 0.000000 |
| 30 | `[0, 3, 9, 16, 20, 28, 30]` | `[0, 1, 16, 21, 24, 28, 30]` | `[0, 3, 17, 19, 25, 26, 30]` | 0.000000 |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| ATH-GAUNTLET-SIDON-20-30-2026-04-30_RESULTS.json | Structured transfer-operator results |
| ATH-GAUNTLET-SIDON-20-30-2026-04-30_REPORT.md | Human-readable report |
| ATH-GAUNTLET-SIDON-20-30-2026-04-30_RESULTS.sha256 | Integrity checksum |
