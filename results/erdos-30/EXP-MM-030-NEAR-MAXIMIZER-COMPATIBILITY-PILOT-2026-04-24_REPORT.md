# EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24 — Near-Maximizer Sidon Compatibility Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Scan window | n = 10 through n = 30 |
| Regime | exact size layers |A| = h(n) and |A| = h(n)-1 |
| Comparison target | whether the one-sided prefix/mass compatibility signal survives one layer below exact maximizers |

## Probe Question

> Does the one-sided compatibility pattern persist when the scan includes near-maximizers, or is it an exact-optimizer artifact?

## Answer

On the `h(n)` layer, the prefix-best and density-adjusted-mass-best witnesses coincide in 5 of 21 values of n.

The mass-best witness has numerically zero prefix residual in 16 of 21 values, while the prefix-best witness has zero density-adjusted mass deviation in 1 of 21 values.

In normalized units, the mass-best witness has prefix cost at most 0.0884, with mean 0.0141. The prefix-best witness has density-adjusted mass cost as high as 0.2025, with mean 0.0835.

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 20 of 21 values. The exceptional n-values are [14].

The best joint witness equals mass-best in 16 of 21 values, equals prefix-best in 5 of 21 values, and is a third witness in 5 of 21 values.

On the `h(n)-1` layer, the prefix-best and density-adjusted-mass-best witnesses coincide in 0 of 21 values of n.

The mass-best witness has numerically zero prefix residual in 2 of 21 values, while the prefix-best witness has zero density-adjusted mass deviation in 2 of 21 values.

In normalized units, the mass-best witness has prefix cost at most 0.2062, with mean 0.1058. The prefix-best witness has density-adjusted mass cost as high as 0.1495, with mean 0.0663.

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 9 of 21 values. The exceptional n-values are [10, 14, 15, 16, 20, 21, 22, 23, 24, 28, 29, 30].

The best joint witness equals mass-best in 2 of 21 values, equals prefix-best in 2 of 21 values, and is a third witness in 17 of 21 values.

## Per-n Layer Summary

| n | layer | size | sets | same best | mass-best prefix ratio | prefix-best mass ratio | joint score |
|---|---|---|---|---|---|---|---|
| 10 | h(n) | 4 | 110 | False | 0.0649 | 0.1687 | 0.0000 |
| 10 | h(n)-1 | 3 | 140 | False | 0.1117 | 0.0000 | 0.0000 |
| 11 | h(n) | 5 | 4 | True | 0.0000 | 0.1110 | 0.1110 |
| 11 | h(n)-1 | 4 | 190 | False | 0.0000 | 0.1295 | 0.0185 |
| 12 | h(n) | 5 | 22 | False | 0.0000 | 0.0328 | 0.0000 |
| 12 | h(n)-1 | 4 | 304 | False | 0.0000 | 0.1313 | 0.0000 |
| 13 | h(n) | 5 | 68 | False | 0.0000 | 0.0588 | 0.0000 |
| 13 | h(n)-1 | 4 | 466 | False | 0.0418 | 0.1029 | 0.0147 |
| 14 | h(n) | 5 | 156 | True | 0.0000 | 0.0000 | 0.0000 |
| 14 | h(n)-1 | 4 | 676 | False | 0.1250 | 0.0266 | 0.0000 |
| 15 | h(n) | 5 | 320 | False | 0.0000 | 0.0483 | 0.0000 |
| 15 | h(n)-1 | 4 | 958 | False | 0.1054 | 0.0362 | 0.0121 |
| 16 | h(n) | 5 | 584 | False | 0.0884 | 0.0884 | 0.0000 |
| 16 | h(n)-1 | 4 | 1,312 | False | 0.1768 | 0.0000 | 0.0000 |
| 17 | h(n) | 6 | 8 | True | 0.0000 | 0.0305 | 0.0305 |
| 17 | h(n)-1 | 5 | 1,008 | False | 0.0413 | 0.1016 | 0.0000 |
| 18 | h(n) | 6 | 24 | True | 0.0000 | 0.0188 | 0.0188 |
| 18 | h(n)-1 | 5 | 1,622 | False | 0.0604 | 0.1316 | 0.0000 |
| 19 | h(n) | 6 | 80 | False | 0.0000 | 0.0611 | 0.0087 |
| 19 | h(n)-1 | 5 | 2,520 | False | 0.0488 | 0.0698 | 0.0000 |
| 20 | h(n) | 6 | 206 | False | 0.0000 | 0.0976 | 0.0000 |
| 20 | h(n)-1 | 5 | 3,734 | False | 0.1111 | 0.0650 | 0.0000 |
| 21 | h(n) | 6 | 504 | False | 0.0000 | 0.1292 | 0.0076 |
| 21 | h(n)-1 | 5 | 5,428 | False | 0.1684 | 0.0608 | 0.0000 |
| 22 | h(n) | 6 | 1,004 | False | 0.0000 | 0.1569 | 0.0000 |
| 22 | h(n)-1 | 5 | 7,612 | False | 0.1545 | 0.0143 | 0.0000 |
| 23 | h(n) | 6 | 1,910 | False | 0.0131 | 0.1811 | 0.0067 |
| 23 | h(n)-1 | 5 | 10,488 | False | 0.2062 | 0.0268 | 0.0000 |
| 24 | h(n) | 6 | 3,380 | False | 0.0745 | 0.2025 | 0.0000 |
| 24 | h(n)-1 | 5 | 14,126 | False | 0.1922 | 0.0127 | 0.0000 |
| 25 | h(n) | 7 | 10 | True | 0.0000 | 0.0239 | 0.0239 |
| 25 | h(n)-1 | 6 | 5,688 | False | 0.0598 | 0.1495 | 0.0060 |
| 26 | h(n) | 7 | 34 | False | 0.0000 | 0.0113 | 0.0000 |
| 26 | h(n)-1 | 6 | 9,036 | False | 0.0521 | 0.0793 | 0.0000 |
| 27 | h(n) | 7 | 98 | False | 0.0000 | 0.0430 | 0.0000 |
| 27 | h(n)-1 | 6 | 14,106 | False | 0.1009 | 0.1022 | 0.0054 |
| 28 | h(n) | 7 | 282 | False | 0.0542 | 0.0717 | 0.0000 |
| 28 | h(n)-1 | 6 | 21,190 | False | 0.1467 | 0.0614 | 0.0000 |
| 29 | h(n) | 7 | 760 | False | 0.0000 | 0.0975 | 0.0000 |
| 29 | h(n)-1 | 6 | 31,158 | False | 0.1899 | 0.0634 | 0.0049 |
| 30 | h(n) | 7 | 1,618 | False | 0.0000 | 0.1210 | 0.0000 |
| 30 | h(n)-1 | 6 | 44,370 | False | 0.1286 | 0.0279 | 0.0000 |

## Interpretation

This packet is a robustness check for the finite Pareto story. If the `h(n)-1` layer keeps the same one-sided direction, the compatibility signal is less likely to be a pure exact-optimizer artifact. If it weakens or reverses there, the current Pareto note should stay explicitly finite and maximizer-scoped.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24_RESULTS.sha256 | Integrity checksum |
