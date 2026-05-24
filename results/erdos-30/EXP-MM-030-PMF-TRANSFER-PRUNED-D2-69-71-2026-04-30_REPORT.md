# EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 69 through n = 71 |
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
| 69 | 8 | 165499468 | 139686424 |
| 70 | 8 | 189156721 | 173484140 |
| 71 | 8 | 215603601 | 214886212 |

## Parity Gate

Checked 3 n-values. h(n) matched in 3. Maximizer counts matched in 3. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 69 | 10 | 66412 | 11.103633 | 16623922 | 122996090 | 1 | 139686424 | 0 | 0.001481 |
| 70 | 10 | 117202 | 11.671654 | 22801688 | 150565250 | 1 | 173484140 | 0 | 0.000000 |
| 71 | 10 | 203840 | 12.225091 | 31058406 | 183623966 | 1 | 214886212 | 0 | 0.001424 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 69 | `[0, 3, 11, 20, 33, 43, 62, 64, 68, 69]` | `[0, 1, 5, 28, 43, 49, 57, 60, 67, 69]` | `[0, 3, 11, 21, 33, 49, 62, 64, 68, 69]` | 0.001481 |
| 70 | `[0, 4, 12, 21, 32, 39, 55, 65, 68, 70]` | `[0, 1, 5, 28, 44, 50, 58, 61, 68, 70]` | `[0, 4, 12, 31, 44, 46, 49, 60, 69, 70]` | 0.000000 |
| 71 | `[0, 4, 13, 21, 40, 45, 56, 68, 70, 71]` | `[0, 1, 6, 31, 44, 53, 55, 63, 67, 70]` | `[0, 4, 13, 23, 34, 51, 63, 65, 66, 71]` | 0.001424 |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30_RESULTS.sha256 | Integrity checksum |
