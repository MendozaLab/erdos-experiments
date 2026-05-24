# EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 71 through n = 71 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 203840 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.4732 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6221 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 71 | 10 | 203840 | 0 | 9 | 2185808767 | 314565085 | 0.4732 | 0.6221 | 0.0000 | 0.0527 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 71 | prefix | 1 | 0.0000 | 0.0071 | 0.0071 | [0, 4, 13, 21, 40, 45, 56, 68, 70, 71] |
| 71 | prefix | 2 | 0.0000 | 0.0328 | 0.0328 | [0, 4, 13, 21, 41, 44, 59, 60, 66, 71] |
| 71 | prefix | 3 | 0.0000 | 0.0242 | 0.0242 | [0, 4, 13, 23, 29, 49, 63, 64, 66, 71] |
| 71 | mass | 1 | 0.1444 | 0.0014 | 0.1458 | [0, 1, 6, 31, 44, 53, 55, 63, 67, 70] |
| 71 | mass | 2 | 0.1444 | 0.0014 | 0.1458 | [0, 1, 6, 33, 42, 53, 55, 63, 67, 70] |
| 71 | mass | 3 | 0.1444 | 0.0014 | 0.1458 | [0, 1, 6, 33, 45, 47, 55, 64, 68, 71] |
| 71 | joint | 1 | 0.0000 | 0.0014 | 0.0014 | [0, 4, 13, 23, 34, 51, 63, 65, 66, 71] |
| 71 | joint | 2 | 0.0000 | 0.0014 | 0.0014 | [0, 4, 13, 28, 44, 46, 54, 65, 66, 71] |
| 71 | joint | 3 | 0.0000 | 0.0014 | 0.0014 | [0, 4, 13, 33, 43, 45, 48, 64, 70, 71] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29_RESULTS.sha256 | Integrity checksum |
