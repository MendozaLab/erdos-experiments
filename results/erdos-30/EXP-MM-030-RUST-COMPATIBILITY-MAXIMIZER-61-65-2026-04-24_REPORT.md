# EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 61 through n = 65 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 8950 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.3581 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6372 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 5 of 5 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 3 of 5 values, equals the prefix-best witness in 0 of 5 values, and is a third witness in 2 of 5 values.

## Per-n Summary

| n | h(n) | maximizers | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|
| 61 | 10 | 152 | 0.1904 | 0.6336 | 0.0000 | 0.0895 |
| 62 | 10 | 398 | 0.2568 | 0.6316 | 0.0000 | 0.1029 |
| 63 | 10 | 1022 | 0.2981 | 0.6295 | 0.0266 | 0.1158 |
| 64 | 10 | 2360 | 0.3416 | 0.6372 | 0.0000 | 0.1281 |
| 65 | 10 | 5018 | 0.3581 | 0.6348 | 0.0259 | 0.1399 |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24_RESULTS.sha256 | Integrity checksum |
