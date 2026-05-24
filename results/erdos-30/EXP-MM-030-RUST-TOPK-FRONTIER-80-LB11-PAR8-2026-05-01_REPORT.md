# EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 80 through n = 80 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 4030 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.3279 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5905 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 1 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 80, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 3, 11, 20, 36, 46, 59, 73, 74, 78, 80].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 80 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 3, 11, 20, 36, 46, 59, 73, 74, 78, 80] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 80 | 11 | 4030 | 0 | 11 | 8903300745 | 1221326452 | 0.3279 | 0.5905 | 0.0000 | 0.0024 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 80 | prefix | 1 | 0.0000 | 0.0411 | 0.0411 | [0, 1, 9, 21, 45, 48, 52, 62, 67, 78, 80] |
| 80 | prefix | 2 | 0.0000 | 0.0628 | 0.0628 | [0, 1, 9, 23, 30, 41, 61, 65, 67, 77, 80] |
| 80 | prefix | 3 | 0.0000 | 0.0145 | 0.0145 | [0, 1, 9, 23, 43, 47, 62, 68, 75, 78, 80] |
| 80 | mass | 1 | 0.0649 | 0.0000 | 0.0649 | [0, 2, 11, 24, 41, 51, 57, 69, 72, 76, 77] |
| 80 | mass | 2 | 0.0745 | 0.0000 | 0.0745 | [0, 3, 5, 24, 39, 50, 62, 66, 72, 79, 80] |
| 80 | mass | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 3, 11, 20, 36, 46, 59, 73, 74, 78, 80] |
| 80 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 3, 11, 20, 36, 46, 59, 73, 74, 78, 80] |
| 80 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 4, 11, 26, 43, 49, 51, 67, 70, 79, 80] |
| 80 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [1, 5, 11, 33, 36, 44, 59, 60, 73, 78, 80] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
