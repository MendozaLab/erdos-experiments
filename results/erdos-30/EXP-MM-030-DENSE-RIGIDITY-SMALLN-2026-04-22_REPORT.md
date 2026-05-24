# EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22 — Exact Dense Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Scan window | n = 10 through n = 50 |
| Dense regime | all Sidon sets A ⊆ [0,n] with |A| = floor(sqrt(n)) |
| Comparison target | ordered profile, prefix discrepancy, and mass center from the current Lean external-interface package |

## Probe Question

> In the exact small-n regime where we can enumerate every floor-sqrt(n) Sidon set, do the theorem-aligned profile, prefix, and mass quantities already look geometry-constrained?

## Answer

The exact scan shows two things at once. First, the Lean dense-Sidon regime is genuinely narrower than the true extremal regime in this window: for 41 of 41 scanned values of n, the exact maximum h(n) sits above floor(sqrt(n)). Second, inside the floor-sqrt(n) regime the profile, prefix, and mass quantities are already tightly organized enough to look like a real rigidity signal rather than noise.

Across the full scan we enumerated 25,189,976 exact dense Sidon sets. The largest observed prefix residual beyond the deterministic sqrt(n) step was 0.7798 in n^(7/8) units, and the largest observed mass deviation was 0.6408 in n^(11/8) units. That does not prove the literature-scale theorem locally, but it is exactly the kind of bounded small-n behavior you would want to see before investing in a new discrepancy argument.

## Per-n Summary

| n | floor(sqrt n) | h(n) | gap | dense sets | best mass dev | worst mass dev | best prefix residual | worst prefix residual |
|---|---|---|---|---|---|---|---|---|
| 10 | 3 | 4 | 1 | 140 | 0.0263 | 14.9737 | 0.0000 | 3.3246 |
| 11 | 3 | 5 | 2 | 190 | 0.1003 | 15.8997 | 0.0000 | 3.6834 |
| 12 | 3 | 5 | 2 | 250 | 0.2154 | 16.7846 | 0.0000 | 4.5359 |
| 13 | 3 | 5 | 2 | 322 | 0.3667 | 17.6333 | 0.0000 | 5.3944 |
| 14 | 3 | 5 | 2 | 406 | 0.4499 | 18.4499 | 0.0000 | 6.2583 |
| 15 | 3 | 5 | 2 | 504 | 0.2379 | 19.2379 | 0.0000 | 7.1270 |
| 16 | 4 | 5 | 1 | 1,312 | 0.0000 | 29.0000 | 0.0000 | 6.0000 |
| 17 | 4 | 6 | 2 | 1,762 | 0.2311 | 30.2311 | 0.0000 | 6.3693 |
| 18 | 4 | 6 | 2 | 2,310 | 0.4264 | 31.4264 | 0.0000 | 6.7574 |
| 19 | 4 | 6 | 2 | 2,984 | 0.4110 | 32.5890 | 0.0000 | 7.6411 |
| 20 | 4 | 6 | 2 | 3,786 | 0.2786 | 33.7214 | 0.0000 | 8.5279 |
| 21 | 4 | 6 | 2 | 4,750 | 0.1742 | 34.8258 | 0.0000 | 9.4174 |
| 22 | 4 | 6 | 2 | 5,874 | 0.0958 | 35.9042 | 0.0000 | 10.3096 |
| 23 | 4 | 6 | 2 | 7,196 | 0.0417 | 36.9583 | 0.0000 | 11.2042 |
| 24 | 4 | 6 | 2 | 8,720 | 0.0102 | 37.9898 | 0.0000 | 12.1010 |
| 25 | 5 | 7 | 2 | 18,744 | 0.0000 | 52.0000 | 0.0000 | 9.0000 |
| 26 | 5 | 7 | 2 | 24,390 | 0.4853 | 53.4853 | 0.0000 | 9.3961 |
| 27 | 5 | 7 | 2 | 31,436 | 0.0577 | 54.9423 | 0.0000 | 9.8038 |
| 28 | 5 | 7 | 2 | 39,914 | 0.3725 | 56.3725 | 0.0000 | 10.7085 |
| 29 | 5 | 7 | 2 | 50,212 | 0.2225 | 57.7775 | 0.0000 | 11.6148 |
| 30 | 5 | 7 | 2 | 62,390 | 0.1584 | 59.1584 | 0.0000 | 12.5228 |
| 31 | 5 | 7 | 2 | 76,932 | 0.4835 | 60.5165 | 0.0000 | 13.4322 |
| 32 | 5 | 7 | 2 | 93,918 | 0.1472 | 61.8528 | 0.0000 | 14.3431 |
| 33 | 5 | 7 | 2 | 113,960 | 0.1684 | 63.1684 | 0.0000 | 15.2554 |
| 34 | 5 | 8 | 3 | 137,058 | 0.4643 | 64.4643 | 0.0000 | 16.1690 |
| 35 | 5 | 8 | 3 | 163,896 | 0.2588 | 65.7412 | 0.0000 | 17.0839 |
| 36 | 6 | 8 | 2 | 262,156 | 0.0000 | 84.0000 | 0.0000 | 13.0000 |
| 37 | 6 | 8 | 2 | 336,960 | 0.2620 | 85.7380 | 0.0000 | 13.4138 |
| 38 | 6 | 8 | 2 | 427,488 | 0.4527 | 87.4527 | 0.0000 | 13.8356 |
| 39 | 6 | 8 | 2 | 538,690 | 0.1450 | 89.1450 | 0.0000 | 14.7550 |
| 40 | 6 | 8 | 2 | 671,604 | 0.1843 | 90.8157 | 0.0000 | 15.6754 |
| 41 | 6 | 8 | 2 | 831,926 | 0.4656 | 92.4656 | 0.0000 | 16.5969 |
| 42 | 6 | 8 | 2 | 1,021,238 | 0.0956 | 94.0956 | 0.0000 | 17.5193 |
| 43 | 6 | 8 | 2 | 1,246,604 | 0.2938 | 95.7062 | 0.0000 | 18.4426 |
| 44 | 6 | 9 | 3 | 1,510,056 | 0.2982 | 97.2982 | 0.0000 | 19.3668 |
| 45 | 6 | 9 | 3 | 1,820,580 | 0.1277 | 98.8723 | 0.0000 | 20.2918 |
| 46 | 6 | 9 | 3 | 2,179,480 | 0.4289 | 100.4289 | 0.0000 | 21.2177 |
| 47 | 6 | 9 | 3 | 2,597,758 | 0.0313 | 101.9687 | 0.0000 | 22.1443 |
| 48 | 6 | 9 | 3 | 3,078,578 | 0.4923 | 103.4923 | 0.0000 | 23.0718 |
| 49 | 7 | 9 | 2 | 3,444,944 | 0.0000 | 123.0000 | 0.0000 | 18.0000 |
| 50 | 7 | 9 | 2 | 4,368,558 | 0.0101 | 124.9899 | 0.0000 | 18.3553 |

## Strongest Small-n Witnesses

For n = 46, the mass-best dense set is [16, 18, 19, 23, 29, 37] with deviation 0.4289, and the prefix-best dense set is [7, 14, 20, 28, 32, 37] with residual 0.0000.
For n = 47, the mass-best dense set is [17, 18, 20, 25, 29, 35] with deviation 0.0313, and the prefix-best dense set is [7, 14, 20, 28, 32, 37] with residual 0.0000.
For n = 48, the mass-best dense set is [17, 18, 20, 24, 29, 37] with deviation 0.4923, and the prefix-best dense set is [7, 14, 20, 28, 32, 37] with residual 0.0000.
For n = 49, the mass-best dense set is [17, 19, 23, 24, 33, 36, 44] with deviation 0.0000, and the prefix-best dense set is [8, 15, 21, 29, 33, 38, 49] with residual 0.0000.
For n = 50, the mass-best dense set is [17, 19, 22, 23, 31, 38, 48] with deviation 0.0101, and the prefix-best dense set is [8, 15, 21, 29, 33, 38, 49] with residual 0.0000.

## Interpretation

The right reading is not that small-n exact data proves the Balasubramanian-Dutta scale. It does not. The right reading is that once we restrict to the same floor-sqrt(n) dense regime used by the current Lean package, the explicit affine and mass centers are not fighting the data. They are organizing it.

The other important outcome is methodological. If the extremal regime usually sits above floor(sqrt(n)), then future local theorem targets need to stay explicit about whether they are studying true maximizers or the dense-below-sqrt(n) corridor. That scope distinction is mathematically load-bearing and the exact scan makes it visible.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22_REPORT.md | Human-readable report |
| EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22_RESULTS.sha256 | Integrity checksum |
