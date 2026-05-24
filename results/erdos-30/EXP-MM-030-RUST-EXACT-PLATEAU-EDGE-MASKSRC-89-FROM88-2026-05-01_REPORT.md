# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 89 through n = 89 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-88-FROM87-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 10 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0788 in n^(7/8) units, and the worst observed mass deviation stayed below 0.4860 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [89].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 89, score = 0.028179, prefix ratio = 0.000000, mass ratio = 0.028179, witness = [4, 13, 14, 21, 34, 46, 49, 60, 65, 83, 87, 89].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 89 | yes | 0.028179 | 0.000000 | 0.028179 | [4, 13, 14, 21, 34, 46, 49, 60, 65, 83, 87, 89] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 89 | 12 | 10 | 0 | 12 | 33871634479 | 4449308156 | 0.0788 | 0.4860 | 0.0000 | 0.0282 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 89 | 88 | 8 | 10 | 8 | 8 | 10 | 0 | true |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 89 | prefix | 1 | 0.0000 | 0.0574 | 0.0574 | [4, 6, 10, 28, 33, 44, 47, 59, 72, 79, 80, 89] |
| 89 | prefix | 2 | 0.0000 | 0.0282 | 0.0282 | [4, 13, 14, 21, 34, 46, 49, 60, 65, 83, 87, 89] |
| 89 | prefix | 3 | 0.0197 | 0.0824 | 0.1021 | [3, 5, 9, 27, 32, 43, 46, 58, 71, 78, 79, 88] |
| 89 | mass | 1 | 0.0000 | 0.0282 | 0.0282 | [4, 13, 14, 21, 34, 46, 49, 60, 65, 83, 87, 89] |
| 89 | mass | 2 | 0.0197 | 0.0532 | 0.0729 | [3, 12, 13, 20, 33, 45, 48, 59, 64, 82, 86, 88] |
| 89 | mass | 3 | 0.0000 | 0.0574 | 0.0574 | [4, 6, 10, 28, 33, 44, 47, 59, 72, 79, 80, 89] |
| 89 | joint | 1 | 0.0000 | 0.0282 | 0.0282 | [4, 13, 14, 21, 34, 46, 49, 60, 65, 83, 87, 89] |
| 89 | joint | 2 | 0.0000 | 0.0574 | 0.0574 | [4, 6, 10, 28, 33, 44, 47, 59, 72, 79, 80, 89] |
| 89 | joint | 3 | 0.0197 | 0.0532 | 0.0729 | [3, 12, 13, 20, 33, 45, 48, 59, 64, 82, 86, 88] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01_RESULTS.sha256 | Integrity checksum |
