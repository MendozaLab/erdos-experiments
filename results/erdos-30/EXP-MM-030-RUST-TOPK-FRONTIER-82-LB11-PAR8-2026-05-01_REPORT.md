# EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 82 through n = 82 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 15958 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.3761 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5996 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 82, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 82 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 82 | 11 | 15958 | 0 | 11 | 12518977464 | 1698049081 | 0.3761 | 0.5996 | 0.0000 | 0.0023 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 82 | prefix | 1 | 0.0000 | 0.1005 | 0.1005 | [0, 1, 10, 21, 34, 37, 59, 63, 65, 77, 82] |
| 82 | prefix | 2 | 0.0000 | 0.0421 | 0.0421 | [0, 1, 10, 21, 45, 48, 53, 60, 76, 78, 82] |
| 82 | prefix | 3 | 0.0000 | 0.0561 | 0.0561 | [0, 1, 10, 23, 35, 50, 52, 68, 71, 76, 82] |
| 82 | mass | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82] |
| 82 | mass | 2 | 0.0212 | 0.0000 | 0.0212 | [0, 1, 14, 22, 37, 48, 65, 72, 75, 77, 81] |
| 82 | mass | 3 | 0.0212 | 0.0000 | 0.0212 | [0, 1, 14, 25, 40, 45, 62, 72, 74, 78, 81] |
| 82 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 11, 27, 45, 49, 62, 68, 70, 77, 82] |
| 82 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 16, 24, 36, 39, 64, 71, 77, 81, 82] |
| 82 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 3, 10, 31, 47, 48, 53, 67, 71, 80, 82] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
