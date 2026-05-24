# EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 83 through n = 83 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 30510 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.3885 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5980 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 1 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 83, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 2, 14, 21, 29, 57, 60, 73, 77, 82, 83].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 83 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 2, 14, 21, 29, 57, 60, 73, 77, 82, 83] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 83 | 11 | 30510 | 0 | 11 | 14819152913 | 1999097010 | 0.3885 | 0.5980 | 0.0000 | 0.0023 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 83 | prefix | 1 | 0.0000 | 0.0712 | 0.0712 | [0, 2, 11, 21, 29, 45, 57, 60, 77, 82, 83] |
| 83 | prefix | 2 | 0.0000 | 0.0253 | 0.0253 | [0, 2, 11, 21, 33, 46, 60, 75, 76, 80, 83] |
| 83 | prefix | 3 | 0.0000 | 0.0345 | 0.0345 | [0, 2, 11, 21, 37, 51, 55, 68, 75, 80, 83] |
| 83 | mass | 1 | 0.1071 | 0.0000 | 0.1071 | [0, 1, 5, 27, 41, 52, 65, 71, 73, 80, 83] |
| 83 | mass | 2 | 0.0862 | 0.0000 | 0.0862 | [0, 1, 6, 32, 45, 54, 62, 65, 69, 81, 83] |
| 83 | mass | 3 | 0.0652 | 0.0000 | 0.0652 | [0, 1, 7, 26, 36, 50, 66, 70, 78, 81, 83] |
| 83 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 14, 21, 29, 57, 60, 73, 77, 82, 83] |
| 83 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 14, 31, 34, 47, 55, 73, 77, 82, 83] |
| 83 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 15, 22, 38, 48, 66, 69, 77, 78, 83] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
