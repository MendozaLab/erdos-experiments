# EXP-MM-030-RUST-TOPK-FRONTIER-74-LB11-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-74-LB11-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 74 through n = 74 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 34 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1294 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5831 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [74].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 74, score = 0.013452, prefix ratio = 0.000000, mass ratio = 0.013452, witness = [0, 2, 8, 25, 40, 43, 53, 62, 69, 73, 74].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 74 | yes | 0.013452 | 0.000000 | 0.013452 | [0, 2, 8, 25, 40, 43, 53, 62, 69, 73, 74] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 74 | 11 | 34 | 0 | 11 | 3107632106 | 441618147 | 0.1294 | 0.5831 | 0.0000 | 0.0430 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 74 | prefix | 1 | 0.0000 | 0.0565 | 0.0565 | [0, 1, 9, 24, 30, 35, 55, 57, 67, 71, 74] |
| 74 | prefix | 2 | 0.0000 | 0.0377 | 0.0377 | [0, 1, 11, 17, 29, 42, 50, 65, 69, 72, 74] |
| 74 | prefix | 3 | 0.0000 | 0.0161 | 0.0161 | [0, 1, 14, 19, 35, 45, 57, 65, 68, 72, 74] |
| 74 | mass | 1 | 0.0000 | 0.0135 | 0.0135 | [0, 2, 8, 25, 40, 43, 53, 62, 69, 73, 74] |
| 74 | mass | 2 | 0.0000 | 0.0161 | 0.0161 | [0, 1, 14, 19, 35, 45, 57, 65, 68, 72, 74] |
| 74 | mass | 3 | 0.0000 | 0.0377 | 0.0377 | [0, 1, 11, 17, 29, 42, 50, 65, 69, 72, 74] |
| 74 | joint | 1 | 0.0000 | 0.0135 | 0.0135 | [0, 2, 8, 25, 40, 43, 53, 62, 69, 73, 74] |
| 74 | joint | 2 | 0.0000 | 0.0161 | 0.0161 | [0, 1, 14, 19, 35, 45, 57, 65, 68, 72, 74] |
| 74 | joint | 3 | 0.0000 | 0.0377 | 0.0377 | [0, 1, 11, 17, 29, 42, 50, 65, 69, 72, 74] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-74-LB11-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-74-LB11-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-74-LB11-2026-05-01_RESULTS.sha256 | Integrity checksum |
