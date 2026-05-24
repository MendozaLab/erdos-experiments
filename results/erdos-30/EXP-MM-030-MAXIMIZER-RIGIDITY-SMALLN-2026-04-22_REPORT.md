# EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Scan window | n = 10 through n = 50 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Probe Question

> On the actual h(n) sets, not just the floor-sqrt(n) corridor, do the new general prefix and mass templates still organize the small-n geometry?

## Answer

Yes in the limited but exact sense that matters here. Across 76,368 exact maximizers in the scan window, the general prefix wrapper and the explicit mass center still track the data at a bounded small-n scale.

The worst observed prefix residual after subtracting the new general drift `max(| |A| - sqrt(n) |, 1) * sqrt(n)` stayed below 0.5303 in n^(7/8) units, and the worst observed mass deviation stayed below 0.9505 in n^(11/8) units. That is not a proof of sharpness, but it does show that the new theorem is pointed at the right regime rather than only at the floor-sqrt(n) corridor.

## Per-n Summary

| n | h(n) | maximizers | gap = ||A|-sqrt(n)| | prefix drift | best prefix residual | worst prefix residual | best mass dev | worst mass dev |
|---|---|---|---|---|---|---|---|---|
| 10 | 4 | 110 | 0.8377 | 3.1623 | 0.0000 | 3.4868 | 2.6228 | 20.6228 |
| 11 | 5 | 4 | 1.6834 | 5.5831 | 0.0000 | 0.3668 | 19.7494 | 24.7494 |
| 12 | 5 | 22 | 1.5359 | 5.3205 | 0.0000 | 2.0718 | 14.9615 | 28.9615 |
| 13 | 5 | 68 | 1.3944 | 5.0278 | 0.0000 | 3.3944 | 12.0833 | 31.0833 |
| 14 | 5 | 156 | 1.2583 | 4.7083 | 0.0000 | 4.2583 | 9.1249 | 33.1249 |
| 15 | 5 | 320 | 1.1270 | 4.3649 | 0.0000 | 5.1270 | 6.0948 | 35.0948 |
| 16 | 5 | 584 | 1.0000 | 4.0000 | 0.0000 | 6.0000 | 3.0000 | 37.0000 |
| 17 | 6 | 8 | 1.8769 | 7.7386 | 0.0000 | 0.8769 | 28.5852 | 42.5852 |
| 18 | 6 | 24 | 1.7574 | 7.4558 | 0.0000 | 2.5147 | 23.0955 | 47.0955 |
| 19 | 6 | 80 | 1.6411 | 7.1534 | 0.0000 | 3.2822 | 19.5369 | 49.5369 |
| 20 | 6 | 206 | 1.5279 | 6.8328 | 0.0000 | 4.0557 | 15.9149 | 51.9149 |
| 21 | 6 | 504 | 1.4174 | 6.4955 | 0.0000 | 5.8348 | 12.2341 | 54.2341 |
| 22 | 6 | 1,004 | 1.3096 | 6.1425 | 0.0000 | 6.6192 | 8.4987 | 56.4987 |
| 23 | 6 | 1,910 | 1.2042 | 5.7750 | 0.0000 | 7.4083 | 4.7125 | 58.7125 |
| 24 | 6 | 3,380 | 1.1010 | 5.3939 | 0.0000 | 8.2020 | 0.8786 | 60.8786 |
| 25 | 7 | 10 | 2.0000 | 10.0000 | 0.0000 | 2.0000 | 42.0000 | 63.0000 |
| 26 | 7 | 34 | 1.9010 | 9.6931 | 0.0000 | 3.8020 | 37.7725 | 65.7725 |
| 27 | 7 | 98 | 1.8038 | 9.3731 | 0.0000 | 4.6077 | 29.4923 | 72.4923 |
| 28 | 7 | 282 | 1.7085 | 9.0405 | 0.0000 | 5.4170 | 25.1621 | 75.1621 |
| 29 | 7 | 760 | 1.6148 | 8.6962 | 0.0000 | 7.2297 | 20.7846 | 77.7846 |
| 30 | 7 | 1,618 | 1.5228 | 8.3406 | 0.0000 | 8.0455 | 16.3623 | 80.3623 |
| 31 | 7 | 3,334 | 1.4322 | 7.9744 | 0.0000 | 8.8645 | 11.8974 | 82.8974 |
| 32 | 7 | 6,360 | 1.3431 | 7.5980 | 0.0000 | 9.6863 | 7.3919 | 85.3919 |
| 33 | 7 | 11,482 | 1.2554 | 7.2119 | 0.0000 | 10.5109 | 2.8478 | 87.8478 |
| 34 | 8 | 2 | 2.1690 | 12.6476 | 0.0000 | 1.6762 | 54.9143 | 92.9143 |
| 35 | 8 | 22 | 2.0839 | 12.3286 | 0.0000 | 2.3357 | 49.9789 | 95.9789 |
| 36 | 8 | 70 | 2.0000 | 12.0000 | 0.0000 | 4.0000 | 45.0000 | 99.0000 |
| 37 | 8 | 214 | 1.9172 | 11.6621 | 0.0000 | 5.7517 | 39.9795 | 101.9795 |
| 38 | 8 | 540 | 1.8356 | 11.3153 | 0.0000 | 6.5068 | 34.9189 | 104.9189 |
| 39 | 8 | 1,250 | 1.7550 | 10.9600 | 0.0000 | 8.5100 | 25.8199 | 111.8199 |
| 40 | 8 | 2,718 | 1.6754 | 10.5964 | 0.0000 | 9.3509 | 20.6840 | 114.6840 |
| 41 | 8 | 5,712 | 1.5969 | 10.2250 | 0.0000 | 10.1938 | 15.5125 | 117.5125 |
| 42 | 8 | 10,910 | 1.5193 | 9.8459 | 0.0000 | 12.0385 | 10.3067 | 120.3067 |
| 43 | 8 | 20,418 | 1.4426 | 9.4595 | 0.0000 | 12.8851 | 5.0678 | 123.0678 |
| 44 | 9 | 2 | 2.3668 | 15.6992 | 0.0000 | 0.0000 | 92.4962 | 108.4962 |
| 45 | 9 | 12 | 2.2918 | 15.3738 | 0.0000 | 3.8754 | 80.8692 | 117.8692 |
| 46 | 9 | 30 | 2.2177 | 15.0410 | 0.0000 | 4.6530 | 68.2048 | 128.2048 |
| 47 | 9 | 90 | 2.1443 | 14.7009 | 0.0000 | 5.4330 | 61.5045 | 132.5045 |
| 48 | 9 | 230 | 2.0718 | 14.3538 | 0.0000 | 6.2872 | 55.7691 | 135.7691 |
| 49 | 9 | 562 | 2.0000 | 14.0000 | 0.0000 | 8.0000 | 50.0000 | 139.0000 |
| 50 | 9 | 1,228 | 1.9289 | 13.6396 | 0.0000 | 8.7868 | 43.1981 | 143.1981 |

## Strongest Small-n Witnesses

For n = 46, the mass-best maximizer is [0, 5, 13, 25, 27, 36, 42, 43, 46] with deviation 68.2048, and the prefix-best maximizer is [2, 5, 11, 19, 21, 34, 41, 45, 46] with residual 0.0000.
For n = 47, the mass-best maximizer is [0, 4, 12, 26, 31, 37, 44, 46, 47] with deviation 61.5045, and the prefix-best maximizer is [3, 6, 12, 20, 22, 35, 42, 46, 47] with residual 0.0000.
For n = 48, the mass-best maximizer is [1, 5, 13, 27, 32, 38, 45, 47, 48] with deviation 55.7691, and the prefix-best maximizer is [4, 7, 13, 21, 23, 36, 43, 47, 48] with residual 0.0000.
For n = 49, the mass-best maximizer is [2, 6, 14, 28, 33, 39, 46, 48, 49] with deviation 50.0000, and the prefix-best maximizer is [5, 8, 14, 22, 24, 37, 44, 48, 49] with residual 0.0000.
For n = 50, the mass-best maximizer is [0, 4, 19, 29, 36, 41, 47, 49, 50] with deviation 43.1981, and the prefix-best maximizer is [6, 9, 15, 23, 25, 38, 45, 49, 50] with residual 0.0000.

## Interpretation

This does not mean the general theorem is sharp. It means the new regime correction was the right one: once the theorem is re-parameterized to see true maximizers, its prefix and mass centers stop missing the exact data for purely scope reasons.

The next mathematical question is no longer whether the theorem is pointed at the right class of sets. It is whether the `max(gap, 1) * sqrt(n)` drift can be compressed further on true maximizers, or whether the observed residual scale is already close to the right endpoint bookkeeping.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22_REPORT.md | Human-readable report |
| EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22_RESULTS.sha256 | Integrity checksum |
