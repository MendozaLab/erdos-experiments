# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 85 through n = 85 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-PLATEAU-EDGE-MASKSRC-84-FROM83-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 2 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.0000 in n^(7/8) units, and the worst observed mass deviation stayed below 0.4806 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 0 of 1 values. Nonzero joint-frontier n-values are [85].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 85, score = 0.078936, prefix ratio = 0.000000, mass ratio = 0.078936, witness = [0, 9, 10, 17, 30, 42, 45, 56, 61, 79, 83, 85].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 85 | yes | 0.078936 | 0.000000 | 0.078936 | [0, 9, 10, 17, 30, 42, 45, 56, 61, 79, 83, 85] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 85 | 12 | 2 | 0 | 11 | 20667370552 | 2758279975 | 0.0000 | 0.4806 | 0.0000 | 0.0789 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 85 | 84 | 56110 | 81710 | 0 | 0 | 0 | 2 | false |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 85 | prefix | 1 | 0.0000 | 0.1101 | 0.1101 | [0, 2, 6, 24, 29, 40, 43, 55, 68, 75, 76, 85] |
| 85 | prefix | 2 | 0.0000 | 0.0789 | 0.0789 | [0, 9, 10, 17, 30, 42, 45, 56, 61, 79, 83, 85] |
| 85 | mass | 1 | 0.0000 | 0.0789 | 0.0789 | [0, 9, 10, 17, 30, 42, 45, 56, 61, 79, 83, 85] |
| 85 | mass | 2 | 0.0000 | 0.1101 | 0.1101 | [0, 2, 6, 24, 29, 40, 43, 55, 68, 75, 76, 85] |
| 85 | joint | 1 | 0.0000 | 0.0789 | 0.0789 | [0, 9, 10, 17, 30, 42, 45, 56, 61, 79, 83, 85] |
| 85 | joint | 2 | 0.0000 | 0.1101 | 0.1101 | [0, 2, 6, 24, 29, 40, 43, 55, 68, 75, 76, 85] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01_RESULTS.sha256 | Integrity checksum |
