# EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 76 through n = 76 |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 214 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1803 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5819 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [76].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 76, score = 0.002593, prefix ratio = 0.000000, mass ratio = 0.002593, witness = [0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 76 | yes | 0.002593 | 0.000000 | 0.002593 | [0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 76 | 11 | 214 | 0 | 11 | 4436801435 | 622974564 | 0.1803 | 0.5819 | 0.0000 | 0.0026 |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 76 | prefix | 1 | 0.0000 | 0.0986 | 0.0986 | [0, 1, 8, 21, 33, 36, 47, 52, 70, 74, 76] |
| 76 | prefix | 2 | 0.0000 | 0.0804 | 0.0804 | [0, 1, 9, 21, 28, 39, 45, 61, 71, 74, 76] |
| 76 | prefix | 3 | 0.0000 | 0.0337 | 0.0337 | [0, 1, 14, 19, 30, 42, 52, 67, 69, 73, 76] |
| 76 | mass | 1 | 0.0000 | 0.0026 | 0.0026 | [0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76] |
| 76 | mass | 2 | 0.0000 | 0.0052 | 0.0052 | [0, 5, 17, 23, 32, 36, 52, 66, 73, 74, 76] |
| 76 | mass | 3 | 0.0000 | 0.0078 | 0.0078 | [0, 4, 9, 24, 32, 46, 57, 63, 73, 75, 76] |
| 76 | joint | 1 | 0.0000 | 0.0026 | 0.0026 | [0, 1, 16, 24, 34, 45, 51, 64, 71, 73, 76] |
| 76 | joint | 2 | 0.0000 | 0.0052 | 0.0052 | [0, 5, 17, 23, 32, 36, 52, 66, 73, 74, 76] |
| 76 | joint | 3 | 0.0000 | 0.0078 | 0.0078 | [0, 4, 9, 24, 32, 46, 57, 63, 73, 75, 76] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01_RESULTS.sha256 | Integrity checksum |
