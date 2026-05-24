# EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 70 through n = 70 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 117202 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.4592 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6248 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 0 of 1 values of n. The exceptional n-values are [70].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 1 of 1 values.

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 70 | 10 | 117202 | 0 | 9 | 1823945271 | 264202675 | 0.4592 | 0.6248 | 0.0616 | 0.0407 |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24_RESULTS.sha256 | Integrity checksum |
