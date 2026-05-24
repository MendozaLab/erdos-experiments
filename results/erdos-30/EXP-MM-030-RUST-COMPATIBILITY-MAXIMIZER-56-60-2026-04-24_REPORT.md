# EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 56 through n = 60 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 92 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1743 in n^(7/8) units, and the worst observed mass deviation stayed below 0.6417 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 4 of 5 values of n. The exceptional n-values are [57].

The best joint witness equals the mass-best witness in 3 of 5 values, equals the prefix-best witness in 2 of 5 values, and is a third witness in 1 of 5 values.

## Per-n Summary

| n | h(n) | maximizers | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|
| 56 | 10 | 4 | 0.0325 | 0.6417 | 0.0000 | 0.0118 |
| 57 | 10 | 6 | 0.0582 | 0.6403 | 0.0291 | 0.0289 |
| 58 | 10 | 10 | 0.0859 | 0.6388 | 0.0286 | 0.0451 |
| 59 | 10 | 18 | 0.1489 | 0.6372 | 0.0564 | 0.0606 |
| 60 | 10 | 54 | 0.1743 | 0.6355 | 0.0000 | 0.0754 |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24_RESULTS.sha256 | Integrity checksum |
