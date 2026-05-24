# EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23 |
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

## Prefix Drift Calibration

| Drift / recentering ansatz | Worst residual in n^(7/8) units | Worst n | Worst t |
|---|---|---|---|
| 0 | 0.8927 | 13 | 6 |
| sqrt(n) | 0.5875 | 43 | 17 |
| abs(card-sqrt(n)) * sqrt(n) | 0.5334 | 10 | 6 |
| max(abs(card-sqrt(n)), 1) * sqrt(n) | 0.5303 | 16 | 6 |
| affine endpoint recentered | 0.8856 | 50 | 32 |
| density-adjusted affine | 0.6199 | 24 | 6 |

Neither the endpoint-bridge recentering nor the density-adjusted affine center beats the current theorem drift on this window. That means the remaining slack is not explained by a simple affine reanchoring alone, and any future no-drift program has to be more structural than just changing slope or endpoint.

## Endpoint vs Interior

| n | h(n) | endpoint bridge | density slope | raw max location | recentered max location | density-adjusted max location | worst density-adjusted ratio | raw witness set | density-adjusted witness set |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 4 | -2.6491 | 2.5000 | 6 | 6 | 3 | 0.6001 | [0, 2, 5, 6] | [0, 2, 3, 10] |
| 11 | 5 | -5.5831 | 2.2000 | 4 | 7 | 1 | 0.4171 | [0, 3, 4, 9, 11] | [0, 1, 4, 9, 11] |
| 12 | 5 | -5.3205 | 2.4000 | 3 | 8 | 3 | 0.4775 | [0, 2, 3, 8, 12] | [0, 2, 3, 8, 12] |
| 13 | 5 | -5.0278 | 2.6000 | 6 | 9 | 3 | 0.5088 | [0, 2, 5, 6, 13] | [0, 2, 3, 9, 13] |
| 14 | 5 | -4.7083 | 2.8000 | 6 | 10 | 3 | 0.5365 | [0, 2, 5, 6, 14] | [0, 2, 3, 10, 14] |
| 15 | 5 | -4.3649 | 3.0000 | 6 | 11 | 6 | 0.5611 | [0, 2, 5, 6, 15] | [0, 2, 5, 6, 15] |
| 16 | 5 | -4.0000 | 3.2000 | 6 | 12 | 6 | 0.6010 | [0, 2, 5, 6, 16] | [0, 2, 5, 6, 16] |
| 17 | 6 | -7.7386 | 2.8333 | 12 | 10 | 1 | 0.3912 | [0, 1, 4, 10, 12, 17] | [0, 1, 8, 12, 14, 17] |
| 18 | 6 | -7.4558 | 3.0000 | 7 | 11 | 3 | 0.4784 | [0, 2, 6, 7, 15, 18] | [0, 1, 3, 8, 14, 18] |
| 19 | 6 | -7.1534 | 3.1667 | 7 | 12 | 3 | 0.4943 | [0, 2, 6, 7, 16, 19] | [0, 2, 3, 10, 15, 19] |
| 20 | 6 | -6.8328 | 3.3333 | 7 | 13 | 3 | 0.5090 | [0, 4, 6, 7, 15, 20] | [0, 2, 3, 10, 15, 19] |
| 21 | 6 | -6.4955 | 3.5000 | 6 | 14 | 6 | 0.5574 | [0, 2, 5, 6, 14, 21] | [0, 2, 5, 6, 14, 21] |
| 22 | 6 | -6.1425 | 3.6667 | 6 | 15 | 6 | 0.5797 | [0, 2, 5, 6, 15, 22] | [0, 2, 5, 6, 15, 22] |
| 23 | 6 | -5.7750 | 3.8333 | 6 | 16 | 6 | 0.6005 | [0, 2, 5, 6, 16, 23] | [0, 2, 5, 6, 16, 23] |
| 24 | 6 | -5.3939 | 4.0000 | 6 | 17 | 6 | 0.6199 | [0, 2, 5, 6, 17, 24] | [0, 2, 5, 6, 17, 24] |
| 25 | 7 | -10.0000 | 3.5714 | 3 | 18 | 3 | 0.4614 | [0, 2, 3, 10, 16, 21, 25] | [0, 2, 3, 10, 16, 21, 25] |
| 26 | 7 | -9.6931 | 3.7143 | 12 | 19 | 3 | 0.4706 | [0, 1, 7, 9, 12, 22, 26] | [0, 2, 3, 10, 16, 21, 25] |
| 27 | 7 | -9.3731 | 3.8571 | 12 | 20 | 3 | 0.4793 | [0, 1, 7, 9, 12, 22, 26] | [0, 2, 3, 10, 16, 21, 25] |
| 28 | 7 | -9.0405 | 4.0000 | 12 | 21 | 7 | 0.4875 | [0, 3, 4, 10, 12, 23, 28] | [0, 2, 6, 7, 16, 24, 27] |
| 29 | 7 | -8.6962 | 4.1429 | 11 | 22 | 11 | 0.5103 | [0, 2, 7, 10, 11, 23, 29] | [0, 2, 7, 10, 11, 23, 29] |
| 30 | 7 | -8.3406 | 4.2857 | 11 | 23 | 6 | 0.5682 | [0, 2, 7, 10, 11, 24, 30] | [0, 1, 4, 6, 14, 23, 30] |
| 31 | 7 | -7.9744 | 4.4286 | 11 | 24 | 6 | 0.5805 | [0, 2, 7, 10, 11, 25, 31] | [0, 2, 5, 6, 16, 24, 31] |
| 32 | 7 | -7.5980 | 4.5714 | 11 | 25 | 6 | 0.5921 | [0, 2, 7, 10, 11, 26, 32] | [0, 2, 5, 6, 16, 24, 31] |
| 33 | 7 | -7.2119 | 4.7143 | 11 | 21 | 6 | 0.6032 | [0, 3, 4, 9, 11, 23, 33] | [0, 2, 5, 6, 18, 26, 33] |
| 34 | 8 | -12.6476 | 4.2500 | 9 | 22 | 4 | 0.3999 | [0, 1, 4, 9, 15, 22, 32, 34] | [0, 1, 4, 9, 15, 22, 32, 34] |
| 35 | 8 | -12.3286 | 4.3750 | 9 | 23 | 4 | 0.4066 | [0, 2, 5, 9, 15, 23, 34, 35] | [0, 3, 4, 12, 22, 28, 33, 35] |
| 36 | 8 | -12.0000 | 4.5000 | 14 | 24 | 3 | 0.4565 | [0, 3, 5, 13, 14, 20, 32, 36] | [0, 2, 3, 14, 21, 27, 31, 36] |
| 37 | 8 | -11.6621 | 4.6250 | 13 | 25 | 7 | 0.4881 | [0, 1, 6, 10, 13, 21, 35, 37] | [0, 1, 5, 7, 16, 24, 34, 37] |
| 38 | 8 | -11.3153 | 4.7500 | 13 | 26 | 7 | 0.4976 | [0, 1, 6, 10, 13, 21, 35, 37] | [0, 1, 5, 7, 16, 24, 34, 37] |
| 39 | 8 | -10.9600 | 4.8750 | 18 | 27 | 7 | 0.5067 | [0, 1, 3, 8, 14, 18, 30, 39] | [0, 2, 6, 7, 16, 28, 36, 39] |
| 40 | 8 | -10.5964 | 5.0000 | 18 | 28 | 6 | 0.5550 | [0, 1, 3, 8, 14, 18, 30, 39] | [0, 2, 5, 6, 16, 25, 33, 40] |
| 41 | 8 | -10.2250 | 5.1250 | 18 | 29 | 6 | 0.5626 | [0, 1, 3, 8, 14, 18, 30, 39] | [0, 2, 5, 6, 16, 25, 33, 40] |
| 42 | 8 | -9.8459 | 5.2500 | 17 | 30 | 11 | 0.5793 | [0, 1, 8, 12, 14, 17, 32, 42] | [0, 2, 7, 10, 11, 24, 36, 42] |
| 43 | 8 | -9.4595 | 5.3750 | 17 | 31 | 11 | 0.5908 | [0, 2, 7, 13, 16, 17, 35, 43] | [0, 2, 7, 10, 11, 24, 36, 42] |
| 44 | 9 | -15.6992 | 4.8889 | 44 | 32 | 5 | 0.3526 | [0, 3, 9, 17, 19, 32, 39, 43, 44] | [0, 1, 5, 12, 25, 27, 35, 41, 44] |
| 45 | 9 | -15.3738 | 5.0000 | 21 | 33 | 9 | 0.3934 | [0, 3, 9, 16, 20, 21, 35, 43, 45] | [0, 3, 7, 9, 21, 29, 34, 44, 45] |
| 46 | 9 | -15.0410 | 5.1111 | 21 | 34 | 9 | 0.4015 | [0, 3, 9, 16, 20, 21, 35, 43, 45] | [0, 3, 7, 9, 21, 29, 34, 44, 45] |
| 47 | 9 | -14.7009 | 5.2222 | 21 | 35 | 3 | 0.4361 | [0, 3, 9, 16, 20, 21, 35, 43, 45] | [0, 1, 3, 10, 16, 21, 35, 43, 47] |
| 48 | 9 | -14.3538 | 5.3333 | 14 | 30 | 8 | 0.4507 | [0, 4, 6, 13, 14, 25, 30, 45, 48] | [0, 3, 7, 8, 21, 31, 37, 46, 48] |
| 49 | 9 | -14.0000 | 5.4444 | 13 | 31 | 7 | 0.4906 | [0, 1, 8, 10, 13, 24, 30, 45, 49] | [0, 1, 5, 7, 19, 29, 32, 40, 49] |
| 50 | 9 | -13.6396 | 5.5556 | 20 | 32 | 7 | 0.4965 | [0, 1, 3, 11, 15, 20, 36, 43, 49] | [0, 2, 6, 7, 16, 24, 35, 47, 50] |

## Mass Center Calibration

| Center ansatz | Worst residual in n^(11/8) units | Worst n | Witness set |
|---|---|---|---|
| sqrt(n)-centered mass template | 0.9505 | 12 | [0, 1, 3, 7, 12] |
| density-adjusted mass template | 0.5904 | 10 | [0, 1, 4, 6] |

The density-adjusted mass center beats the old sqrt(n)-centered mass template on this window. That means the mass slack really does look like a centering problem, not just a coarse theorem constant.

## Observable Split

> Do the same exact maximizers optimize both the best prefix observable and the best density-adjusted mass observable, or do the two observables genuinely pull toward different witness sets?

Mostly they split. In only 8 of the 41 scanned values of n does the prefix-best maximizer coincide with the density-adjusted-mass-best maximizer.

The strongest split in this window occurs at n = 24: the prefix-best witness is [6, 10, 16, 21, 23, 24], the mass-best witness is [6, 8, 12, 13, 21, 24], and the combined normalized split score is 0.2770.

The best joint compromise witness appears at n = 10, where [3, 4, 8, 10] minimizes the summed normalized prefix-plus-mass score at 0.0000.

So the exact data does not force a total decoupling, but it does say that one affine coordinate is not obviously organizing both observables at once.

## Per-n Summary

| n | h(n) | maximizers | gap = ||A|-sqrt(n)| | endpoint bridge | density slope | prefix drift | best prefix residual | worst prefix residual | best density-adjusted residual | worst density-adjusted residual | best mass dev | worst mass dev | best density-adjusted mass dev | worst density-adjusted mass dev |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 4 | 110 | 0.8377 | -2.6491 | 2.5000 | 3.1623 | 0.0000 | 3.4868 | 1.5000 | 4.5000 | 2.6228 | 20.6228 | 0.0000 | 14.0000 |
| 11 | 5 | 4 | 1.6834 | -5.5831 | 2.2000 | 5.5831 | 0.0000 | 0.3668 | 2.4000 | 3.4000 | 19.7494 | 24.7494 | 3.0000 | 8.0000 |
| 12 | 5 | 22 | 1.5359 | -5.3205 | 2.4000 | 5.3205 | 0.0000 | 2.0718 | 2.2000 | 4.2000 | 14.9615 | 28.9615 | 0.0000 | 13.0000 |
| 13 | 5 | 68 | 1.3944 | -5.0278 | 2.6000 | 5.0278 | 0.0000 | 3.3944 | 2.2000 | 4.8000 | 12.0833 | 31.0833 | 0.0000 | 16.0000 |
| 14 | 5 | 156 | 1.2583 | -4.7083 | 2.8000 | 4.7083 | 0.0000 | 4.2583 | 2.4000 | 5.4000 | 9.1249 | 33.1249 | 0.0000 | 19.0000 |
| 15 | 5 | 320 | 1.1270 | -4.3649 | 3.0000 | 4.3649 | 0.0000 | 5.1270 | 2.0000 | 6.0000 | 6.0948 | 35.0948 | 0.0000 | 22.0000 |
| 16 | 5 | 584 | 1.0000 | -4.0000 | 3.2000 | 4.0000 | 0.0000 | 6.0000 | 2.4000 | 6.8000 | 3.0000 | 37.0000 | 0.0000 | 25.0000 |
| 17 | 6 | 8 | 1.8769 | -7.7386 | 2.8333 | 7.7386 | 0.0000 | 0.8769 | 3.6667 | 4.6667 | 28.5852 | 42.5852 | 1.5000 | 15.5000 |
| 18 | 6 | 24 | 1.7574 | -7.4558 | 3.0000 | 7.4558 | 0.0000 | 2.5147 | 3.0000 | 6.0000 | 23.0955 | 47.0955 | 1.0000 | 21.0000 |
| 19 | 6 | 80 | 1.6411 | -7.1534 | 3.1667 | 7.1534 | 0.0000 | 3.2822 | 3.1667 | 6.5000 | 19.5369 | 49.5369 | 0.5000 | 24.5000 |
| 20 | 6 | 206 | 1.5279 | -6.8328 | 3.3333 | 6.8328 | 0.0000 | 4.0557 | 3.0000 | 7.0000 | 15.9149 | 51.9149 | 0.0000 | 28.0000 |
| 21 | 6 | 504 | 1.4174 | -6.4955 | 3.5000 | 6.4955 | 0.0000 | 5.8348 | 3.0000 | 8.0000 | 12.2341 | 54.2341 | 0.5000 | 31.5000 |
| 22 | 6 | 1,004 | 1.3096 | -6.1425 | 3.6667 | 6.1425 | 0.0000 | 6.6192 | 3.0000 | 8.6667 | 8.4987 | 56.4987 | 0.0000 | 35.0000 |
| 23 | 6 | 1,910 | 1.2042 | -5.7750 | 3.8333 | 5.7750 | 0.0000 | 7.4083 | 3.1667 | 9.3333 | 4.7125 | 58.7125 | 0.5000 | 38.5000 |
| 24 | 6 | 3,380 | 1.1010 | -5.3939 | 4.0000 | 5.3939 | 0.0000 | 8.2020 | 3.0000 | 10.0000 | 0.8786 | 60.8786 | 0.0000 | 42.0000 |
| 25 | 7 | 10 | 2.0000 | -10.0000 | 3.5714 | 10.0000 | 0.0000 | 2.0000 | 5.2857 | 7.7143 | 42.0000 | 63.0000 | 2.0000 | 23.0000 |
| 26 | 7 | 34 | 1.9010 | -9.6931 | 3.7143 | 9.6931 | 0.0000 | 3.8020 | 5.4286 | 8.1429 | 37.7725 | 65.7725 | 0.0000 | 27.0000 |
| 27 | 7 | 98 | 1.8038 | -9.3731 | 3.8571 | 9.3731 | 0.0000 | 4.6077 | 4.7143 | 8.5714 | 29.4923 | 72.4923 | 0.0000 | 35.0000 |
| 28 | 7 | 282 | 1.7085 | -9.0405 | 4.0000 | 9.0405 | 0.0000 | 5.4170 | 5.0000 | 9.0000 | 25.1621 | 75.1621 | 0.0000 | 39.0000 |
| 29 | 7 | 760 | 1.6148 | -8.6962 | 4.1429 | 8.6962 | 0.0000 | 7.2297 | 4.4286 | 9.7143 | 20.7846 | 77.7846 | 0.0000 | 43.0000 |
| 30 | 7 | 1,618 | 1.5228 | -8.3406 | 4.2857 | 8.3406 | 0.0000 | 8.0455 | 4.4286 | 11.1429 | 16.3623 | 80.3623 | 0.0000 | 47.0000 |
| 31 | 7 | 3,334 | 1.4322 | -7.9744 | 4.4286 | 7.9744 | 0.0000 | 8.8645 | 4.4286 | 11.7143 | 11.8974 | 82.8974 | 0.0000 | 51.0000 |
| 32 | 7 | 6,360 | 1.3431 | -7.5980 | 4.5714 | 7.5980 | 0.0000 | 9.6863 | 4.1429 | 12.2857 | 7.3919 | 85.3919 | 0.0000 | 55.0000 |
| 33 | 7 | 11,482 | 1.2554 | -7.2119 | 4.7143 | 7.2119 | 0.0000 | 10.5109 | 4.2857 | 12.8571 | 2.8478 | 87.8478 | 0.0000 | 59.0000 |
| 34 | 8 | 2 | 2.1690 | -12.6476 | 4.2500 | 12.6476 | 0.0000 | 1.6762 | 7.7500 | 8.7500 | 54.9143 | 92.9143 | 2.0000 | 36.0000 |
| 35 | 8 | 22 | 2.0839 | -12.3286 | 4.3750 | 12.3286 | 0.0000 | 2.3357 | 6.1250 | 9.1250 | 49.9789 | 95.9789 | 0.5000 | 40.5000 |
| 36 | 8 | 70 | 2.0000 | -12.0000 | 4.5000 | 12.0000 | 0.0000 | 4.0000 | 6.0000 | 10.5000 | 45.0000 | 99.0000 | 1.0000 | 45.0000 |
| 37 | 8 | 214 | 1.9172 | -11.6621 | 4.6250 | 11.6621 | 0.0000 | 5.7517 | 5.8750 | 11.5000 | 39.9795 | 101.9795 | 0.5000 | 49.5000 |
| 38 | 8 | 540 | 1.8356 | -11.3153 | 4.7500 | 11.3153 | 0.0000 | 6.5068 | 5.5000 | 12.0000 | 34.9189 | 104.9189 | 0.0000 | 54.0000 |
| 39 | 8 | 1,250 | 1.7550 | -10.9600 | 4.8750 | 10.9600 | 0.0000 | 8.5100 | 5.5000 | 12.5000 | 25.8199 | 111.8199 | 0.5000 | 62.5000 |
| 40 | 8 | 2,718 | 1.6754 | -10.5964 | 5.0000 | 10.5964 | 0.0000 | 9.3509 | 5.0000 | 14.0000 | 20.6840 | 114.6840 | 0.0000 | 67.0000 |
| 41 | 8 | 5,712 | 1.5969 | -10.2250 | 5.1250 | 10.2250 | 0.0000 | 10.1938 | 5.2500 | 14.5000 | 15.5125 | 117.5125 | 0.5000 | 71.5000 |
| 42 | 8 | 10,910 | 1.5193 | -9.8459 | 5.2500 | 9.8459 | 0.0000 | 12.0385 | 5.2500 | 15.2500 | 10.3067 | 120.3067 | 0.0000 | 76.0000 |
| 43 | 8 | 20,418 | 1.4426 | -9.4595 | 5.3750 | 9.4595 | 0.0000 | 12.8851 | 4.7500 | 15.8750 | 5.0678 | 123.0678 | 0.5000 | 80.5000 |
| 44 | 9 | 2 | 2.3668 | -15.6992 | 4.8889 | 15.6992 | 0.0000 | 0.0000 | 8.6667 | 9.6667 | 92.4962 | 108.4962 | 14.0000 | 30.0000 |
| 45 | 9 | 12 | 2.2918 | -15.3738 | 5.0000 | 15.3738 | 0.0000 | 3.8754 | 8.0000 | 11.0000 | 80.8692 | 117.8692 | 4.0000 | 41.0000 |
| 46 | 9 | 30 | 2.2177 | -15.0410 | 5.1111 | 15.0410 | 0.0000 | 4.6530 | 7.3333 | 11.4444 | 68.2048 | 128.2048 | 0.0000 | 53.0000 |
| 47 | 9 | 90 | 2.1443 | -14.7009 | 5.2222 | 14.7009 | 0.0000 | 5.4330 | 7.4444 | 12.6667 | 61.5045 | 132.5045 | 2.0000 | 59.0000 |
| 48 | 9 | 230 | 2.0718 | -14.3538 | 5.3333 | 14.3538 | 0.0000 | 6.2872 | 7.0000 | 13.3333 | 55.7691 | 135.7691 | 0.0000 | 64.0000 |
| 49 | 9 | 562 | 2.0000 | -14.0000 | 5.4444 | 14.0000 | 0.0000 | 8.0000 | 7.2222 | 14.7778 | 50.0000 | 139.0000 | 0.0000 | 69.0000 |
| 50 | 9 | 1,228 | 1.9289 | -13.6396 | 5.5556 | 13.6396 | 0.0000 | 8.7868 | 6.7778 | 15.2222 | 43.1981 | 143.1981 | 0.0000 | 75.0000 |

## Strongest Small-n Witnesses

For n = 46, the mass-best maximizer is [0, 5, 13, 25, 27, 36, 42, 43, 46] with deviation 68.2048, and the prefix-best maximizer is [2, 5, 11, 19, 21, 34, 41, 45, 46] with residual 0.0000.
For n = 47, the mass-best maximizer is [0, 4, 12, 26, 31, 37, 44, 46, 47] with deviation 61.5045, and the prefix-best maximizer is [3, 6, 12, 20, 22, 35, 42, 46, 47] with residual 0.0000.
For n = 48, the mass-best maximizer is [1, 5, 13, 27, 32, 38, 45, 47, 48] with deviation 55.7691, and the prefix-best maximizer is [4, 7, 13, 21, 23, 36, 43, 47, 48] with residual 0.0000.
For n = 49, the mass-best maximizer is [2, 6, 14, 28, 33, 39, 46, 48, 49] with deviation 50.0000, and the prefix-best maximizer is [5, 8, 14, 22, 24, 37, 44, 48, 49] with residual 0.0000.
For n = 50, the mass-best maximizer is [0, 4, 19, 29, 36, 41, 47, 49, 50] with deviation 43.1981, and the prefix-best maximizer is [6, 9, 15, 23, 25, 38, 45, 49, 50] with residual 0.0000.

## Interpretation

This does not mean the general theorem is sharp. It means the new regime correction was the right one: once the theorem is re-parameterized to see true maximizers, its prefix and mass centers stop missing the exact data for purely scope reasons.

The next mathematical question is no longer whether the theorem is pointed at the right class of sets. It is whether the remaining slack is best understood as theorem-constant waste or as a need for a more structural center than either the current drift wrapper or the first affine recentering candidates.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23_REPORT.md | Human-readable report |
| EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23_RESULTS.sha256 | Integrity checksum |
