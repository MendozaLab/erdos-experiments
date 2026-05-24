# EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 57 through n = 57 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 6 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0582 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6403 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 0 of 1 values of n. The exceptional n-values are [57].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 1 of 1 values, and is a third witness in 0 of 1 values.

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 57 | 10 | 6 | 0 | 8 | 167555344 | 26505952 | 0.0582 | 0.6403 | 0.0291 | 0.0289 |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29_RESULTS.sha256 | Integrity checksum |
