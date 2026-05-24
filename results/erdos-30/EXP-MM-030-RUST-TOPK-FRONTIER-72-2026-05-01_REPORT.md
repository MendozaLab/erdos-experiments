# EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 72 through n = 72 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 4 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0258 in n^(7/8) units, and the worst observed mass deviation stayed below 0.4862 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 1 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [72].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 72, score = 0.072634, prefix ratio = 0.000000, mass ratio = 0.072634, witness = [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 72 | yes | 0.072634 | 0.000000 | 0.072634 | [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 72 | 11 | 4 | 0 | 9 | 2597444754 | 371655085 | 0.0258 | 0.4862 | 0.0000 | 0.0726 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 72 | prefix | 1 | 0.0000 | 0.1145 | 0.1145 | [0, 1, 9, 19, 24, 31, 52, 56, 58, 69, 72] |
| 72 | prefix | 2 | 0.0000 | 0.0726 | 0.0726 | [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72] |
| 72 | prefix | 3 | 0.0028 | 0.1285 | 0.1313 | [0, 1, 4, 13, 28, 33, 47, 54, 64, 70, 72] |
| 72 | mass | 1 | 0.0000 | 0.0726 | 0.0726 | [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72] |
| 72 | mass | 2 | 0.0258 | 0.0866 | 0.1124 | [0, 3, 14, 16, 20, 41, 48, 53, 63, 71, 72] |
| 72 | mass | 3 | 0.0000 | 0.1145 | 0.1145 | [0, 1, 9, 19, 24, 31, 52, 56, 58, 69, 72] |
| 72 | joint | 1 | 0.0000 | 0.0726 | 0.0726 | [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72] |
| 72 | joint | 2 | 0.0258 | 0.0866 | 0.1124 | [0, 3, 14, 16, 20, 41, 48, 53, 63, 71, 72] |
| 72 | joint | 3 | 0.0000 | 0.1145 | 0.1145 | [0, 1, 9, 19, 24, 31, 52, 56, 58, 69, 72] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01_RESULTS.sha256 | Integrity checksum |
