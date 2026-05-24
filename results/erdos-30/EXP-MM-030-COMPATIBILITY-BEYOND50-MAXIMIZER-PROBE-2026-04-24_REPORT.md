# EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Scan window | n = 51 through n = 55 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Probe Question

> On the actual h(n) sets, not just the floor-sqrt(n) corridor, do the new general prefix and mass templates still organize the small-n geometry?

## Answer

Yes in the limited but exact sense that matters here. Across 40,666 exact maximizers in the scan window, the general prefix wrapper and the explicit mass center still track the data at a bounded small-n scale.

The worst observed prefix residual after subtracting the new general drift `max(| |A| - sqrt(n) |, 1) * sqrt(n)` stayed below 0.4056 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6614 in n^(11/8) units. That is not a proof of sharpness, but it does show that the new theorem is pointed at the right regime rather than only at the floor-sqrt(n) corridor.

## Prefix Drift Calibration

| Drift / recentering ansatz | Worst residual in n^(7/8) units | Worst n | Worst t |
|---|---|---|---|
| 0 | 0.7756 | 54 | 26 |
| sqrt(n) | 0.5516 | 54 | 26 |
| abs(card-sqrt(n)) * sqrt(n) | 0.4056 | 54 | 26 |
| max(abs(card-sqrt(n)), 1) * sqrt(n) | 0.4056 | 54 | 26 |
| affine endpoint recentered | 0.8982 | 55 | 37 |
| density-adjusted affine | 0.5488 | 54 | 6 |

Neither the endpoint-bridge recentering nor the density-adjusted affine center beats the current theorem drift on this window. That means the remaining slack is not explained by a simple affine reanchoring alone, and any future no-drift program has to be more structural than just changing slope or endpoint.

## Endpoint vs Interior

| n | h(n) | endpoint bridge | density slope | raw max location | recentered max location | density-adjusted max location | worst density-adjusted ratio | raw witness set | density-adjusted witness set |
|---|---|---|---|---|---|---|---|---|---|
| 51 | 9 | -13.2729 | 5.6667 | 20 | 33 | 7 | 0.5022 | [0, 4, 7, 9, 19, 20, 37, 45, 51] | [0, 4, 6, 7, 17, 22, 31, 43, 51] |
| 52 | 9 | -12.8999 | 5.7778 | 12 | 34 | 6 | 0.5392 | [0, 1, 4, 10, 12, 27, 34, 47, 52] | [0, 2, 5, 6, 18, 28, 37, 45, 52] |
| 53 | 9 | -12.5210 | 5.8889 | 26 | 35 | 6 | 0.5441 | [0, 1, 5, 15, 18, 24, 26, 46, 53] | [0, 2, 5, 6, 19, 28, 35, 43, 53] |
| 54 | 9 | -12.1362 | 6.0000 | 26 | 36 | 6 | 0.5488 | [0, 1, 5, 15, 18, 24, 26, 46, 53] | [0, 2, 5, 6, 19, 28, 35, 43, 53] |
| 55 | 10 | -19.1620 | 5.5000 | 10 | 37 | 10 | 0.3601 | [0, 1, 6, 10, 23, 26, 34, 41, 53, 55] | [0, 1, 6, 10, 23, 26, 34, 41, 53, 55] |

## Mass Center Calibration

| Center ansatz | Worst residual in n^(11/8) units | Worst n | Witness set |
|---|---|---|---|
| sqrt(n)-centered mass template | 0.6614 | 51 | [0, 2, 5, 9, 15, 23, 34, 35, 51] |
| density-adjusted mass template | 0.3983 | 54 | [0, 2, 5, 9, 15, 23, 34, 35, 51] |

The density-adjusted mass center beats the old sqrt(n)-centered mass template on this window. That means the mass slack really does look like a centering problem, not just a coarse theorem constant.

## Observable Split

> Do the same exact maximizers optimize both the best prefix observable and the best density-adjusted mass observable, or do the two observables genuinely pull toward different witness sets?

Mostly they split. In only 1 of the 5 scanned values of n does the prefix-best maximizer coincide with the density-adjusted-mass-best maximizer.

The strongest split in this window occurs at n = 54: the prefix-best witness is [10, 13, 19, 27, 29, 42, 49, 53, 54], the mass-best witness is [7, 9, 13, 22, 27, 43, 44, 51, 54], and the combined normalized split score is 0.1079.

The best joint compromise witness appears at n = 52, where [6, 10, 12, 19, 33, 36, 41, 51, 52] minimizes the summed normalized prefix-plus-mass score at 0.0000.

So the exact data does not force a total decoupling, but it does say that one affine coordinate is not obviously organizing both observables at once.

## Finite Compatibility Candidate

The first compatibility signal is one-sided. Optimizing density-adjusted mass usually keeps the prefix observable near its best face, while optimizing prefix often leaves visible density-adjusted mass cost.

Using numerical zero tolerance `1e-09`, the mass-best witness has zero prefix residual in 4 of 5 values of n. The prefix-best witness has zero density-adjusted mass deviation in 0 of 5 values.

In normalized units, the mass-best witness has prefix cost at most 0.0321 in n^(7/8) units, with mean 0.0064. The prefix-best witness has density-adjusted mass cost as high as 0.1079 in n^(11/8) units, with mean 0.0698.

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 5 of 5 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 4 of 5 values and also equals the prefix-best witness in 1 of 5 values; these counts can overlap when the same witness optimizes both observables. A third witness is joint-best in 1 of 5 values.

## Per-n Summary

| n | h(n) | maximizers | gap = ||A|-sqrt(n)| | endpoint bridge | density slope | prefix drift | best prefix residual | worst prefix residual | best density-adjusted residual | worst density-adjusted residual | best mass dev | worst mass dev | best density-adjusted mass dev | worst density-adjusted mass dev |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 51 | 9 | 2,690 | 1.8586 | -13.2729 | 5.6667 | 13.2729 | 0.0000 | 9.5757 | 7.0000 | 15.6667 | 36.3643 | 147.3643 | 0.0000 | 81.0000 |
| 52 | 9 | 5,550 | 1.7889 | -12.8999 | 5.7778 | 12.8999 | 0.0000 | 11.1556 | 6.6667 | 17.1111 | 30.4996 | 150.4996 | 0.0000 | 86.0000 |
| 53 | 9 | 11,260 | 1.7199 | -12.5210 | 5.8889 | 12.5210 | 0.0000 | 12.4398 | 6.6667 | 17.5556 | 24.6049 | 153.6049 | 0.0000 | 91.0000 |
| 54 | 9 | 21,164 | 1.6515 | -12.1362 | 6.0000 | 12.1362 | 0.0000 | 13.3031 | 7.0000 | 18.0000 | 18.6811 | 156.6811 | 0.0000 | 96.0000 |
| 55 | 10 | 2 | 2.5838 | -19.1620 | 5.5000 | 19.1620 | 0.0000 | 0.5028 | 11.0000 | 12.0000 | 106.8909 | 158.8909 | 1.5000 | 53.5000 |

## Strongest Small-n Witnesses

For n = 51, the mass-best maximizer is [0, 16, 17, 28, 36, 42, 46, 49, 51] with deviation 36.3643, and the prefix-best maximizer is [7, 10, 16, 24, 26, 39, 46, 50, 51] with residual 0.0000.
For n = 52, the mass-best maximizer is [1, 17, 18, 29, 37, 43, 47, 50, 52] with deviation 30.4996, and the prefix-best maximizer is [8, 11, 17, 25, 27, 40, 47, 51, 52] with residual 0.0000.
For n = 53, the mass-best maximizer is [2, 18, 19, 30, 38, 44, 48, 51, 53] with deviation 24.6049, and the prefix-best maximizer is [9, 12, 18, 26, 28, 41, 48, 52, 53] with residual 0.0000.
For n = 54, the mass-best maximizer is [3, 19, 20, 31, 39, 45, 49, 52, 54] with deviation 18.6811, and the prefix-best maximizer is [10, 13, 19, 27, 29, 42, 49, 53, 54] with residual 0.0000.
For n = 55, the mass-best maximizer is [0, 2, 14, 21, 29, 32, 45, 49, 54, 55] with deviation 106.8909, and the prefix-best maximizer is [0, 2, 14, 21, 29, 32, 45, 49, 54, 55] with residual 0.0000.

## Interpretation

This does not mean the general theorem is sharp. It means the new regime correction was the right one: once the theorem is re-parameterized to see true maximizers, its prefix and mass centers stop missing the exact data for purely scope reasons.

The next mathematical question is no longer whether the theorem is pointed at the right class of sets. It is whether the remaining slack is best understood as theorem-constant waste or as a need for a more structural center than either the current drift wrapper or the first affine recentering candidates.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24_RESULTS.sha256 | Integrity checksum |
