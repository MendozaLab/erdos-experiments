# EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-31-35-2026-04-29 — B_2[2] Transfer-State Scout

## Question

Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?

## Result

Exact finite maxima were computed for 5 rows. The exact h(n) range was 12..12. Field splits occurred in 5 rows: [31, 32, 33, 34, 35].

| n | h(n) | maximizers | nodes | pruned | prefix winner | mass winner | joint winner |
|---|---:|---:|---:|---:|---|---|---|
| 31 | 12 | 40 | 23470108 | 6481565 | `[1, 2, 6, 7, 13, 18, 20, 21, 22, 28, 30, 31]` | `[0, 1, 4, 5, 8, 18, 19, 20, 25, 27, 29, 31]` | `[1, 2, 4, 10, 11, 12, 14, 19, 25, 26, 30, 31]` |
| 32 | 12 | 218 | 33583962 | 9199017 | `[1, 2, 6, 7, 13, 18, 20, 21, 22, 28, 30, 31]` | `[0, 1, 2, 9, 10, 13, 14, 24, 27, 29, 30, 32]` | `[2, 3, 5, 11, 12, 13, 15, 20, 26, 27, 31, 32]` |
| 33 | 12 | 1372 | 49376054 | 13396285 | `[0, 2, 5, 10, 12, 13, 17, 23, 26, 27, 32, 33]` | `[0, 1, 2, 3, 14, 15, 21, 24, 25, 29, 31, 33]` | `[0, 2, 5, 10, 12, 13, 17, 23, 26, 27, 32, 33]` |
| 34 | 12 | 6014 | 73817418 | 19817234 | `[0, 3, 7, 10, 11, 18, 20, 22, 23, 28, 32, 34]` | `[0, 1, 2, 3, 12, 19, 21, 25, 26, 29, 32, 34]` | `[0, 3, 7, 10, 11, 18, 20, 22, 23, 28, 32, 34]` |
| 35 | 12 | 24728 | 111385011 | 29581466 | `[1, 3, 7, 12, 13, 15, 17, 24, 25, 28, 32, 35]` | `[0, 1, 2, 3, 12, 20, 22, 25, 28, 29, 33, 35]` | `[1, 3, 7, 12, 13, 15, 17, 24, 25, 28, 32, 35]` |

## Claim Boundary

This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet.
