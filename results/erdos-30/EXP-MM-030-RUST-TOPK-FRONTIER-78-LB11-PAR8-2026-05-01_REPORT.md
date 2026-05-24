# EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 78 through n = 78 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 970 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.2580 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5878 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 0 of 1 values of n. The exceptional n-values are [78].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 0 of 1 values of n. The first-hit failures are [78].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: [78]. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 78, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 1, 15, 18, 41, 46, 57, 66, 70, 76, 78].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 78 | no | 0.000000 | 0.000000 | 0.000000 | [0, 1, 15, 18, 41, 46, 57, 66, 70, 76, 78] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 78 | 11 | 970 | 0 | 11 | 6300979290 | 874346701 | 0.2580 | 0.5878 | 0.0000 | 0.0000 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 78 | prefix | 1 | 0.0000 | 0.1076 | 0.1076 | [0, 1, 8, 18, 30, 41, 44, 57, 72, 76, 78] |
| 78 | prefix | 2 | 0.0000 | 0.0400 | 0.0400 | [0, 1, 8, 23, 29, 40, 60, 64, 73, 76, 78] |
| 78 | prefix | 3 | 0.0000 | 0.0350 | 0.0350 | [0, 1, 11, 18, 41, 50, 54, 56, 70, 75, 78] |
| 78 | mass | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 15, 18, 41, 46, 57, 66, 70, 76, 78] |
| 78 | mass | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 11, 23, 43, 49, 56, 59, 73, 74, 78] |
| 78 | mass | 3 | 0.0000 | 0.0000 | 0.0000 | [1, 3, 20, 25, 35, 41, 48, 66, 74, 77, 78] |
| 78 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 15, 18, 41, 46, 57, 66, 70, 76, 78] |
| 78 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 11, 23, 43, 49, 56, 59, 73, 74, 78] |
| 78 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [1, 3, 20, 25, 35, 41, 48, 66, 74, 77, 78] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
