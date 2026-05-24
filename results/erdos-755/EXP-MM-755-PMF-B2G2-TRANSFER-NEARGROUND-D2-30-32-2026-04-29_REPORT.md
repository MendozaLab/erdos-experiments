# EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-30-32-2026-04-29 — B_2[2] Transfer-State Scout

## Question

Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?

## Result

Exact finite maxima were computed for 3 rows. The exact h(n) range was 12..12. Field splits occurred in 3 rows: [30, 31, 32].

| n | h(n) | maximizers | h-1 | h-2 | nodes | pruned | prefix winner | mass winner | joint winner |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| 30 | 12 | 6 | 27856 | 1195604 | 17170555 | 4769079 | `[0, 1, 5, 6, 12, 17, 19, 20, 21, 27, 29, 30]` | `[0, 1, 3, 9, 10, 11, 13, 18, 24, 25, 29, 30]` | `[0, 1, 5, 6, 12, 17, 19, 20, 21, 27, 29, 30]` |
| 31 | 12 | 40 | 89612 | 2311480 | 23470108 | 6481565 | `[1, 2, 6, 7, 13, 18, 20, 21, 22, 28, 30, 31]` | `[0, 1, 4, 5, 8, 18, 19, 20, 25, 27, 29, 31]` | `[1, 2, 4, 10, 11, 12, 14, 19, 25, 26, 30, 31]` |
| 32 | 12 | 218 | 240288 | 4195380 | 33583962 | 9199017 | `[1, 2, 6, 7, 13, 18, 20, 21, 22, 28, 30, 31]` | `[0, 1, 2, 9, 10, 13, 14, 24, 27, 29, 30, 32]` | `[2, 3, 5, 11, 12, 13, 15, 20, 26, 27, 31, 32]` |

## Claim Boundary

This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet.
