# EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 75 through n = 75 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 84 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1454 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5826 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 75, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 1, 15, 25, 33, 46, 52, 63, 68, 72, 75].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 75 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 1, 15, 25, 33, 46, 52, 63, 68, 72, 75] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 75 | 11 | 84 | 0 | 11 | 3715690427 | 524838192 | 0.1454 | 0.5826 | 0.0000 | 0.0291 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 75 | prefix | 1 | 0.0000 | 0.0766 | 0.0766 | [0, 1, 9, 22, 27, 33, 56, 58, 68, 72, 75] |
| 75 | prefix | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 15, 25, 33, 46, 52, 63, 68, 72, 75] |
| 75 | prefix | 3 | 0.0000 | 0.0898 | 0.0898 | [0, 1, 16, 18, 28, 37, 41, 61, 67, 72, 75] |
| 75 | mass | 1 | 0.0229 | 0.0000 | 0.0229 | [0, 1, 14, 19, 35, 45, 57, 65, 68, 72, 74] |
| 75 | mass | 2 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 15, 25, 33, 46, 52, 63, 68, 72, 75] |
| 75 | mass | 3 | 0.0229 | 0.0026 | 0.0255 | [0, 2, 8, 25, 40, 43, 53, 62, 69, 73, 74] |
| 75 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 1, 15, 25, 33, 46, 52, 63, 68, 72, 75] |
| 75 | joint | 2 | 0.0229 | 0.0000 | 0.0229 | [0, 1, 14, 19, 35, 45, 57, 65, 68, 72, 74] |
| 75 | joint | 3 | 0.0000 | 0.0238 | 0.0238 | [1, 2, 12, 18, 30, 43, 51, 66, 70, 73, 75] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01_RESULTS.sha256 | Integrity checksum |
