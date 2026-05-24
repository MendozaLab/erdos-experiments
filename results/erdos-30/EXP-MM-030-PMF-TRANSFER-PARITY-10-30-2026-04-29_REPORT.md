# EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 10 through n = 30 |
| Frontier k | 5 |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Parity Gate

Checked 21 n-values. h(n) matched in 11. Maximizer counts matched in 0. Mismatches: [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 4 | 110 | 4.700480 | 140 | 55 | 1 | 317 | 317 | 0.000000 |
| 11 | 5 | 4 | 1.386294 | 190 | 190 | 1 | 463 | 463 | 0.110970 |
| 12 | 5 | 22 | 3.091042 | 304 | 250 | 1 | 668 | 668 | 0.000000 |
| 13 | 5 | 68 | 4.219508 | 466 | 322 | 1 | 962 | 962 | 0.000000 |
| 14 | 5 | 156 | 5.049856 | 676 | 406 | 1 | 1359 | 1359 | 0.000000 |
| 15 | 5 | 320 | 5.768321 | 958 | 504 | 1 | 1919 | 1919 | 0.000000 |
| 16 | 5 | 584 | 6.369901 | 1312 | 616 | 1 | 2666 | 2666 | 0.022097 |
| 17 | 6 | 8 | 2.079442 | 1008 | 1762 | 1 | 3694 | 3694 | 0.030495 |
| 18 | 6 | 24 | 3.178054 | 1622 | 2310 | 1 | 5035 | 5035 | 0.018793 |
| 19 | 6 | 80 | 4.382027 | 2520 | 2984 | 1 | 6845 | 6845 | 0.008723 |
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

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked until the transfer operator also passes the 56-58 and 69-71 defect windows.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29_RESULTS.sha256 | Integrity checksum |
