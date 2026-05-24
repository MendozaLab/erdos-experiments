# EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 79 through n = 79 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 1974 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.2745 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5917 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 79, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 5, 20, 24, 38, 45, 51, 67, 68, 77, 79].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 79 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 5, 20, 24, 38, 45, 51, 67, 68, 77, 79] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 79 | 11 | 1974 | 0 | 11 | 7494653059 | 1034009248 | 0.2745 | 0.5917 | 0.0000 | 0.0000 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 79 | prefix | 1 | 0.0000 | 0.0738 | 0.0738 | [0, 1, 8, 24, 27, 45, 57, 59, 70, 74, 79] |
| 79 | prefix | 2 | 0.0000 | 0.0074 | 0.0074 | [0, 1, 8, 25, 35, 46, 61, 65, 74, 77, 79] |
| 79 | prefix | 3 | 0.0000 | 0.0836 | 0.0836 | [0, 1, 8, 25, 36, 45, 48, 58, 63, 77, 79] |
| 79 | mass | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 5, 20, 24, 38, 45, 51, 67, 68, 77, 79] |
| 79 | mass | 2 | 0.0219 | 0.0000 | 0.0219 | [0, 5, 21, 25, 36, 43, 49, 66, 75, 76, 78] |
| 79 | mass | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 7, 11, 26, 36, 38, 56, 70, 73, 78, 79] |
| 79 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 5, 20, 24, 38, 45, 51, 67, 68, 77, 79] |
| 79 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 7, 11, 26, 36, 38, 56, 70, 73, 78, 79] |
| 79 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 10, 13, 17, 41, 43, 59, 64, 70, 78, 79] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
