# EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 77 through n = 77 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 482 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.2263 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5812 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 0 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 1 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 77, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 2, 13, 16, 37, 44, 59, 67, 71, 76, 77].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 77 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 2, 13, 16, 37, 44, 59, 67, 71, 76, 77] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 77 | 11 | 482 | 0 | 11 | 5290742483 | 738457565 | 0.2263 | 0.5812 | 0.0000 | 0.0025 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 77 | prefix | 1 | 0.0000 | 0.0713 | 0.0713 | [0, 1, 12, 16, 36, 46, 49, 54, 68, 75, 77] |
| 77 | prefix | 2 | 0.0000 | 0.0178 | 0.0178 | [0, 1, 12, 18, 40, 44, 53, 67, 69, 74, 77] |
| 77 | prefix | 3 | 0.0000 | 0.0357 | 0.0357 | [0, 1, 12, 19, 35, 44, 50, 64, 72, 74, 77] |
| 77 | mass | 1 | 0.0447 | 0.0000 | 0.0447 | [0, 2, 11, 24, 40, 45, 57, 65, 71, 72, 75] |
| 77 | mass | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 13, 16, 37, 44, 59, 67, 71, 76, 77] |
| 77 | mass | 3 | 0.0000 | 0.0000 | 0.0000 | [1, 6, 9, 18, 37, 50, 52, 66, 70, 76, 77] |
| 77 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 2, 13, 16, 37, 44, 59, 67, 71, 76, 77] |
| 77 | joint | 2 | 0.0000 | 0.0000 | 0.0000 | [1, 6, 9, 18, 37, 50, 52, 66, 70, 76, 77] |
| 77 | joint | 3 | 0.0000 | 0.0025 | 0.0025 | [1, 4, 11, 26, 35, 37, 58, 64, 72, 76, 77] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
