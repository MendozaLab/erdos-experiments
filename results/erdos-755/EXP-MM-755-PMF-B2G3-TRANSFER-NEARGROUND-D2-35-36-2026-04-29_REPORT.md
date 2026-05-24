# EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-35-36-2026-04-29 — B_2[3] Transfer-State Scout

## Question

Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?

## Result

Exact finite maxima were computed for 2 rows. The exact h(n) range was 16..16. Field splits occurred in 2 rows: [35, 36].

| n | h(n) | maximizers | h-1 | h-2 | nodes | pruned | prefix winner | mass winner | joint winner |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| 35 | 16 | 14 | 300424 | 30868234 | 895841364 | 277856781 | `[0, 1, 4, 5, 8, 16, 20, 21, 23, 25, 26, 27, 32, 33, 34, 35]` | `[0, 1, 2, 3, 4, 14, 18, 20, 21, 23, 25, 26, 29, 30, 33, 35]` | `[0, 2, 5, 6, 9, 10, 12, 14, 15, 17, 21, 31, 32, 33, 34, 35]` |
| 36 | 16 | 154 | 1317756 | 75689678 | 1307823032 | 403253763 | `[0, 1, 3, 7, 10, 13, 14, 15, 17, 25, 26, 30, 31, 34, 35, 36]` | `[0, 1, 2, 3, 4, 16, 17, 20, 22, 23, 26, 27, 28, 31, 33, 35]` | `[0, 1, 3, 7, 10, 13, 14, 15, 17, 25, 26, 30, 31, 34, 35, 36]` |

## Claim Boundary

This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet.
