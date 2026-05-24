# EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 81 through n = 81 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 8214 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.3421 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5964 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 81, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 81 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 81 | 11 | 8214 | 0 | 11 | 10563744181 | 1440880837 | 0.3421 | 0.5964 | 0.0000 | 0.0048 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 81 | prefix | 1 | 0.0000 | 0.0309 | 0.0309 | [0, 1, 9, 20, 35, 51, 56, 68, 74, 78, 81] |
| 81 | prefix | 2 | 0.0000 | 0.1022 | 0.1022 | [0, 1, 9, 21, 34, 36, 58, 62, 65, 76, 81] |
| 81 | prefix | 3 | 0.0000 | 0.0333 | 0.0333 | [0, 1, 9, 21, 35, 52, 59, 62, 75, 77, 81] |
| 81 | mass | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81] |
| 81 | mass | 2 | 0.0214 | 0.0000 | 0.0214 | [0, 1, 9, 23, 43, 47, 62, 68, 75, 78, 80] |
| 81 | mass | 3 | 0.0214 | 0.0000 | 0.0214 | [0, 1, 11, 26, 43, 46, 62, 67, 74, 76, 80] |
| 81 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 9, 21, 46, 53, 57, 63, 76, 79, 81] |
| 81 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 12, 27, 36, 56, 58, 64, 74, 77, 81] |
| 81 | joint | 3 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 9, 26, 36, 48, 61, 66, 77, 80, 81] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
