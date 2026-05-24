# EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 66 through n = 68 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 65646 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.4302 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6329 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 3 of 3 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 2 of 3 values, equals the prefix-best witness in 0 of 3 values, and is a third witness in 1 of 3 values.

## Per-n Summary

| n | h(n) | maximizers | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|
| 66 | 10 | 9994 | 0.3998 | 0.6323 | 0.0256 | 0.1511 |
| 67 | 10 | 19418 | 0.4151 | 0.6329 | 0.0000 | 0.1619 |
| 68 | 10 | 36234 | 0.4302 | 0.6302 | 0.0000 | 0.1723 |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24_RESULTS.sha256 | Integrity checksum |
