# EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 56 through n = 58 |
| Frontier k | 5 |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Parity Gate

Checked 0 n-values. h(n) matched in 0. Maximizer counts matched in 0. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 56 | 10 | 4 | 1.386294 | 69564 | 4813066 | 1 | 33911490 | 33911490 | 0.011840 |
| 57 | 10 | 6 | 1.791759 | 120704 | 6511012 | 1 | 40821696 | 40821696 | 0.009630 |
| 58 | 10 | 10 | 2.302585 | 200946 | 8684372 | 1 | 48993329 | 48993329 | 0.007522 |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29_RESULTS.sha256 | Integrity checksum |
