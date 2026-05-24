# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 87 through n = 87 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-86-FROM85-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 6 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0402 in n^(7/8) units, and the worst observed mass deviation stayed below 0.4836 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [87].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 87, score = 0.052763, prefix ratio = 0.000000, mass ratio = 0.052763, witness = [2, 11, 12, 19, 32, 44, 47, 58, 63, 81, 85, 87].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 87 | yes | 0.052763 | 0.000000 | 0.052763 | [2, 11, 12, 19, 32, 44, 47, 58, 63, 81, 85, 87] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 87 | 12 | 6 | 0 | 12 | 24414456532 | 3240856425 | 0.0402 | 0.4836 | 0.0000 | 0.0528 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 87 | 86 | 4 | 6 | 4 | 4 | 6 | 0 | true |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 87 | prefix | 1 | 0.0000 | 0.0829 | 0.0829 | [2, 4, 8, 26, 31, 42, 45, 57, 70, 77, 78, 87] |
| 87 | prefix | 2 | 0.0000 | 0.0528 | 0.0528 | [2, 11, 12, 19, 32, 44, 47, 58, 63, 81, 85, 87] |
| 87 | prefix | 3 | 0.0201 | 0.1088 | 0.1288 | [1, 3, 7, 25, 30, 41, 44, 56, 69, 76, 77, 86] |
| 87 | mass | 1 | 0.0000 | 0.0528 | 0.0528 | [2, 11, 12, 19, 32, 44, 47, 58, 63, 81, 85, 87] |
| 87 | mass | 2 | 0.0201 | 0.0786 | 0.0987 | [1, 10, 11, 18, 31, 43, 46, 57, 62, 80, 84, 86] |
| 87 | mass | 3 | 0.0000 | 0.0829 | 0.0829 | [2, 4, 8, 26, 31, 42, 45, 57, 70, 77, 78, 87] |
| 87 | joint | 1 | 0.0000 | 0.0528 | 0.0528 | [2, 11, 12, 19, 32, 44, 47, 58, 63, 81, 85, 87] |
| 87 | joint | 2 | 0.0000 | 0.0829 | 0.0829 | [2, 4, 8, 26, 31, 42, 45, 57, 70, 77, 78, 87] |
| 87 | joint | 3 | 0.0201 | 0.0786 | 0.0987 | [1, 10, 11, 18, 31, 43, 46, 57, 62, 80, 84, 86] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01_RESULTS.sha256 | Integrity checksum |
