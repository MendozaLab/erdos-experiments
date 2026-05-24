# EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 84 through n = 84 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 56110 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.4040 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5964 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 1 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 84, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 2, 14, 27, 35, 58, 64, 65, 75, 80, 84].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 84 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 2, 14, 27, 35, 58, 64, 65, 75, 80, 84] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 84 | 11 | 56110 | 0 | 11 | 17521913289 | 2350900638 | 0.4040 | 0.5964 | 0.0000 | 0.0000 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 84 | prefix | 1 | 0.0000 | 0.0520 | 0.0520 | [0, 2, 11, 21, 34, 48, 54, 72, 76, 79, 84] |
| 84 | prefix | 2 | 0.0000 | 0.0271 | 0.0271 | [0, 2, 11, 21, 34, 48, 54, 76, 79, 83, 84] |
| 84 | prefix | 3 | 0.0000 | 0.0249 | 0.0249 | [0, 2, 11, 21, 37, 51, 59, 66, 79, 83, 84] |
| 84 | mass | 1 | 0.0348 | 0.0000 | 0.0348 | [0, 1, 9, 24, 35, 52, 66, 72, 79, 82, 84] |
| 84 | mass | 2 | 0.0348 | 0.0000 | 0.0348 | [0, 1, 9, 34, 39, 55, 57, 70, 74, 81, 84] |
| 84 | mass | 3 | 0.0414 | 0.0000 | 0.0414 | [0, 1, 10, 31, 50, 53, 58, 65, 76, 78, 82] |
| 84 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 14, 27, 35, 58, 64, 65, 75, 80, 84] |
| 84 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 17, 24, 37, 40, 65, 73, 79, 83, 84] |
| 84 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 3, 11, 27, 31, 52, 69, 70, 75, 82, 84] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-84-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
