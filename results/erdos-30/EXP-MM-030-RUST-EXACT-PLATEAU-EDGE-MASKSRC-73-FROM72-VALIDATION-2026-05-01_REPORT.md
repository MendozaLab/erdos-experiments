# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-73-FROM72-VALIDATION-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-73-FROM72-VALIDATION-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 73 through n = 73 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 8 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0407 in n^(7/8) units, and the worst observed mass deviation stayed below 0.4877 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [73].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 73, score = 0.057564, prefix ratio = 0.000000, mass ratio = 0.057564, witness = [1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 73 | yes | 0.057564 | 0.000000 | 0.057564 | [1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 73 | 11 | 8 | 0 | 11 | 2595520615 | 371124979 | 0.0407 | 0.4877 | 0.0000 | 0.0576 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 73 | 72 | 4 | 8 | 4 | 4 | 8 | 0 | true |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 73 | prefix | 1 | 0.0000 | 0.1124 | 0.1124 | [1, 2, 5, 14, 29, 34, 48, 55, 65, 71, 73] |
| 73 | prefix | 2 | 0.0000 | 0.0987 | 0.0987 | [1, 2, 10, 20, 25, 32, 53, 57, 59, 70, 73] |
| 73 | prefix | 3 | 0.0000 | 0.0576 | 0.0576 | [1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73] |
| 73 | mass | 1 | 0.0000 | 0.0576 | 0.0576 | [1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73] |
| 73 | mass | 2 | 0.0172 | 0.0713 | 0.0885 | [1, 4, 15, 17, 21, 42, 49, 54, 64, 72, 73] |
| 73 | mass | 3 | 0.0234 | 0.0877 | 0.1111 | [0, 2, 8, 18, 25, 39, 44, 59, 68, 71, 72] |
| 73 | joint | 1 | 0.0000 | 0.0576 | 0.0576 | [1, 3, 9, 19, 26, 40, 45, 60, 69, 72, 73] |
| 73 | joint | 2 | 0.0172 | 0.0713 | 0.0885 | [1, 4, 15, 17, 21, 42, 49, 54, 64, 72, 73] |
| 73 | joint | 3 | 0.0000 | 0.0987 | 0.0987 | [1, 2, 10, 20, 25, 32, 53, 57, 59, 70, 73] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-73-FROM72-VALIDATION-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-73-FROM72-VALIDATION-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-73-FROM72-VALIDATION-2026-05-01_RESULTS.sha256 | Integrity checksum |
