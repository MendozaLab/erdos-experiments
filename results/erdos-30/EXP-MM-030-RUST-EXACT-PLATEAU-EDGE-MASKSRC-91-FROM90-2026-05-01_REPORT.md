# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 91 through n = 91 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-90-FROM89-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 28 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1291 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5022 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [91].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 91, score = 0.005061, prefix ratio = 0.000000, mass ratio = 0.005061, witness = [6, 15, 16, 23, 36, 48, 51, 62, 67, 85, 89, 91].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 91 | yes | 0.005061 | 0.000000 | 0.005061 | [6, 15, 16, 23, 36, 48, 51, 62, 67, 85, 89, 91] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 91 | 12 | 28 | 0 | 12 | 46797043762 | 6084574527 | 0.1291 | 0.5022 | 0.0000 | 0.0051 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 91 | 90 | 14 | 18 | 14 | 14 | 18 | 10 | true |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 91 | prefix | 1 | 0.0000 | 0.0476 | 0.0476 | [0, 2, 12, 23, 30, 36, 50, 65, 82, 87, 90, 91] |
| 91 | prefix | 2 | 0.0000 | 0.0658 | 0.0658 | [0, 4, 6, 22, 27, 46, 53, 61, 78, 81, 90, 91] |
| 91 | prefix | 3 | 0.0000 | 0.0881 | 0.0881 | [0, 6, 10, 15, 27, 35, 51, 58, 77, 88, 90, 91] |
| 91 | mass | 1 | 0.0000 | 0.0051 | 0.0051 | [6, 15, 16, 23, 36, 48, 51, 62, 67, 85, 89, 91] |
| 91 | mass | 2 | 0.0000 | 0.0091 | 0.0091 | [0, 7, 16, 27, 31, 41, 53, 70, 83, 88, 89, 91] |
| 91 | mass | 3 | 0.0193 | 0.0294 | 0.0487 | [5, 14, 15, 22, 35, 47, 50, 61, 66, 84, 88, 90] |
| 91 | joint | 1 | 0.0000 | 0.0051 | 0.0051 | [6, 15, 16, 23, 36, 48, 51, 62, 67, 85, 89, 91] |
| 91 | joint | 2 | 0.0000 | 0.0091 | 0.0091 | [0, 7, 16, 27, 31, 41, 53, 70, 83, 88, 89, 91] |
| 91 | joint | 3 | 0.0000 | 0.0334 | 0.0334 | [6, 8, 12, 30, 35, 46, 49, 61, 74, 81, 82, 91] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01_RESULTS.sha256 | Integrity checksum |
