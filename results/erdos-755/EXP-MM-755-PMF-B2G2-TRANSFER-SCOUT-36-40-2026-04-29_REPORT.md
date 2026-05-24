# EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-36-40-2026-04-29 — B_2[2] Transfer-State Scout

## Question

Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?

## Result

Exact finite maxima were computed for 5 rows. The exact h(n) range was 13..13. Field splits occurred in 4 rows: [36, 38, 39, 40].

| n | h(n) | maximizers | nodes | pruned | prefix winner | mass winner | joint winner |
|---|---:|---:|---:|---:|---|---|---|
| 36 | 13 | 2 | 157102788 | 41458866 | `[0, 3, 7, 11, 12, 15, 16, 26, 28, 33, 34, 35, 36]` | `[0, 1, 2, 3, 8, 10, 20, 21, 24, 25, 29, 33, 36]` | `[0, 3, 7, 11, 12, 15, 16, 26, 28, 33, 34, 35, 36]` |
| 37 | 13 | 30 | 214410083 | 56210818 | `[0, 1, 5, 6, 10, 15, 22, 23, 24, 26, 34, 35, 37]` | `[0, 1, 5, 6, 10, 15, 22, 23, 24, 26, 34, 35, 37]` | `[0, 1, 5, 6, 10, 15, 22, 23, 24, 26, 34, 35, 37]` |
| 38 | 13 | 238 | 301220806 | 78338404 | `[1, 3, 4, 12, 14, 15, 16, 23, 28, 32, 33, 37, 38]` | `[0, 1, 3, 11, 12, 13, 15, 22, 30, 31, 35, 36, 38]` | `[1, 2, 6, 7, 11, 16, 23, 24, 25, 27, 35, 36, 38]` |
| 39 | 13 | 1570 | 431889545 | 111330539 | `[0, 1, 4, 10, 12, 13, 16, 24, 29, 32, 34, 37, 39]` | `[0, 1, 2, 3, 16, 17, 23, 26, 27, 30, 34, 35, 39]` | `[0, 2, 7, 10, 11, 16, 17, 21, 24, 33, 36, 37, 39]` |
| 40 | 13 | 7400 | 626135408 | 159867698 | `[0, 4, 5, 9, 12, 17, 23, 26, 27, 33, 36, 37, 39]` | `[0, 1, 2, 3, 8, 22, 23, 26, 27, 34, 36, 38, 40]` | `[0, 2, 5, 10, 13, 14, 18, 24, 27, 33, 34, 39, 40]` |

## Claim Boundary

This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet.
