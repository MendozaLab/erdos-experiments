# EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-92-FROM91-2026-05-01 — Exact Maximizer Sidon Rigidity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-92-FROM91-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact enumeration, no sampling |
| Implementation | Rust exact maximizer scanner |
| Scan window | n = 92 through n = 92 |
| Inheritance source | /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01_RESULTS.json |
| Ground-face mask export | true |
| Regime | all exact maximizers A ⊆ [0,n] with |A| = h(n) |
| Comparison target | maximizer-friendly prefix and mass theorems using max(| |A| - sqrt(n) |, 1) * sqrt(n) drift |

## Answer

The Rust scan preserves the Python packet contract while moving the bottleneck into a compiled exact enumerator. Across 60 exact maximizers in this window, the general prefix wrapper and explicit mass center still track the data at bounded small-n scale.

The worst observed prefix residual after subtracting the general drift stayed below 0.1886 in n^(7/8) units, and the worst observed mass deviation stayed below 0.5188 in n^(11/8) units.

## Finite Compatibility Candidate

The mass-best prefix penalty is no larger than the prefix-best density-adjusted mass penalty in 1 of 1 values of n. The exceptional n-values are [].

The best joint witness equals the mass-best witness in 1 of 1 values, equals the prefix-best witness in 0 of 1 values, and is a third witness in 0 of 1 values.

## First-hit vs Face-aware Summary

First-hit compatibility holds in 1 of 1 values of n. The first-hit failures are [].

The top-k joint frontier has near-zero joint score in 1 of 1 values. Nonzero joint-frontier n-values are [].

First-hit failures recovered by the face-aware frontier: []. First-hit failures persisting as nonzero face-level handoff cases: [].

Largest top-k joint score: n = 92, score = 0.000000, prefix ratio = 0.000000, mass ratio = 0.000000, witness = [0, 7, 13, 28, 30, 40, 54, 72, 83, 88, 91, 92].

| n | first-hit compatible? | top-k joint score | top-k prefix ratio | top-k mass ratio | top-k joint witness |
|---|---|---:|---:|---:|---|
| 92 | yes | 0.000000 | 0.000000 | 0.000000 | [0, 7, 13, 28, 30, 40, 54, 72, 83, 88, 91, 92] |

## Per-n Summary

| n | h(n) | maximizers | seed depth | initial lower bound | nodes | prunes | prefix residual max ratio | mass dev max ratio | mass-best prefix ratio | prefix-best mass ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 92 | 12 | 60 | 0 | 12 | 54923840771 | 7105685738 | 0.1886 | 0.5188 | 0.0000 | 0.0000 |

## Plateau Inheritance Probe

The inheritance probe checks whether the current exact maximizer face contains the previous exact face and the previous face shifted by `+1` during the same exact traversal used for h(n).

| n | source n | previous face | inherited union | prev exact present | prev +1 present | inherited present | new face | union contained |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 92 | 91 | 28 | 42 | 28 | 28 | 42 | 18 | true |

## Top-k Frontier Witnesses

| n | frontier | rank | prefix ratio | mass ratio | joint score | witness |
|---|---|---:|---:|---:|---:|---|
| 92 | prefix | 1 | 0.0000 | 0.0578 | 0.0578 | [0, 1, 7, 19, 35, 39, 56, 65, 79, 87, 89, 92] |
| 92 | prefix | 2 | 0.0000 | 0.1097 | 0.1097 | [0, 1, 12, 16, 29, 38, 48, 62, 69, 87, 89, 92] |
| 92 | prefix | 3 | 0.0000 | 0.0219 | 0.0219 | [0, 1, 18, 25, 38, 47, 61, 77, 80, 82, 88, 92] |
| 92 | mass | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 7, 13, 28, 30, 40, 54, 72, 83, 88, 91, 92] |
| 92 | mass | 2 | 0.0191 | 0.0040 | 0.0231 | [0, 7, 16, 27, 31, 41, 53, 70, 83, 88, 89, 91] |
| 92 | mass | 3 | 0.0000 | 0.0060 | 0.0060 | [7, 16, 17, 24, 37, 49, 52, 63, 68, 86, 90, 92] |
| 92 | joint | 1 | 0.0000 | 0.0000 | 0.0000 | [0, 7, 13, 28, 30, 40, 54, 72, 83, 88, 91, 92] |
| 92 | joint | 2 | 0.0000 | 0.0060 | 0.0060 | [7, 16, 17, 24, 37, 49, 52, 63, 68, 86, 90, 92] |
| 92 | joint | 3 | 0.0000 | 0.0199 | 0.0199 | [1, 8, 17, 28, 32, 42, 54, 71, 84, 89, 90, 92] |

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-92-FROM91-2026-05-01_RESULTS.json | Structured exact-enumeration results |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-92-FROM91-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-92-FROM91-2026-05-01_RESULTS.sha256 | Integrity checksum |
