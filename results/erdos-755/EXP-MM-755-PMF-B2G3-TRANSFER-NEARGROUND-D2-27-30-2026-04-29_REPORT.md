# EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-27-30-2026-04-29 — B_2[3] Transfer-State Scout

## Question

Does the Sidon-adjacent bounded-sum deformation show the same field-sensitive ground-face behavior seen in #30?

## Result

Exact finite maxima were computed for 4 rows. The exact h(n) range was 14..14. Field splits occurred in 4 rows: [27, 28, 29, 30].

| n | h(n) | maximizers | h-1 | h-2 | nodes | pruned | prefix winner | mass winner | joint winner |
|---|---:|---:|---:|---:|---:|---:|---|---|---|
| 27 | 14 | 48 | 139090 | 2846606 | 20512275 | 6785688 | `[0, 2, 4, 5, 8, 14, 16, 18, 19, 20, 21, 25, 26, 27]` | `[0, 1, 2, 3, 4, 12, 13, 15, 17, 19, 20, 22, 25, 27]` | `[0, 1, 3, 4, 7, 14, 16, 17, 18, 19, 22, 24, 26, 27]` |
| 28 | 14 | 752 | 495394 | 6316830 | 31647045 | 10408049 | `[0, 2, 4, 6, 8, 10, 15, 16, 18, 19, 22, 25, 27, 28]` | `[0, 1, 2, 3, 4, 13, 16, 17, 18, 21, 22, 24, 27, 28]` | `[0, 2, 4, 6, 8, 10, 15, 16, 18, 19, 22, 25, 27, 28]` |
| 29 | 14 | 8652 | 1614908 | 13413866 | 52680018 | 17164381 | `[0, 2, 4, 6, 9, 11, 15, 17, 18, 21, 22, 25, 28, 29]` | `[0, 1, 2, 3, 4, 14, 16, 17, 19, 22, 23, 25, 28, 29]` | `[0, 2, 4, 6, 8, 12, 14, 17, 18, 19, 21, 26, 27, 28]` |
| 30 | 14 | 52958 | 4417694 | 26461952 | 91765283 | 29546540 | `[1, 3, 5, 7, 9, 11, 16, 17, 19, 20, 23, 26, 28, 29]` | `[0, 1, 2, 3, 4, 10, 18, 19, 21, 23, 24, 26, 29, 30]` | `[0, 3, 5, 7, 10, 11, 13, 17, 18, 21, 22, 26, 27, 30]` |

## Claim Boundary

This is an exact finite computational scout for the B_2[g] transfer state, not an independent theorem or parity check against a prior #755 reference packet.
