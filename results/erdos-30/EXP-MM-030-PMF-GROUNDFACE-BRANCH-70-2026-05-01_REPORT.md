# EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 70 through n = 70 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 125000 |
| Ground-face distance-edge cap | 0 |

## State Model

- `State = { occupied_mask: u128, used_differences_mask: u128, cardinality: u8 }`
- Occupied mask: bit i is 1 iff lattice site i is occupied
- Difference memory: bit d is 1 iff a positive difference d has already been realized
- Occupied suffix: last 8 occupied sites are serialized for representative ground states.
- Transition: skip x always; occupy x iff every new difference |x-a| is absent from used_differences_mask

## Reachability Pruning

Enabled with deficiency `0`. After each site, states are retained only if their current cardinality plus remaining sites can still reach `h(n)-0` using the reference `h(n)` for that row.

| n | min cardinality | pruned states | terminal retained states |
|---|---:|---:|---:|
| 70 | 10 | 260934836 | 117202 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 70 | 10 | 117202 | 11.671654 | 0 | 0 | NA | 117202 | 0 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 70 | `[0, 4, 12, 21, 32, 39, 55, 65, 68, 70]` | `[0, 1, 5, 28, 44, 50, 58, 61, 68, 70]` | `[0, 4, 12, 31, 44, 46, 49, 60, 69, 70]` | 0.000000 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 70 | EXPORTED_ALL | 117202 | 117202 | SKIP_DISTANCE_EDGES_TOO_LARGE | 6868095801 | true | `[33850, 34327, 34396, 34403, 34414, 34560, 34806, 35346, 38522, 38571, 38807, 38977, 38998, 39262, 39307, 39418, 39551, 39564, 39574, 42063, 42293, 42454, 42653, 44438, 44502, 44527, 44711, 45055, 46298, 46331, 46625, 46686, 47877, 47881, 47902, 47986, 48517, 48731, 48735, 48772, 48834, 49215, 49313, 49783, 49897, 50218, 50271, 50353, 50587, 67621, 68005, 68239, 68267, 68450, 68457, 68554, 70931, 72004, 73659, 73665, 73801, 73841, 73851, 74026, 74244, 74263, 74425, 74484, 74524, 74616, 75930, 75931, 75965, 76220, 76648, 77384, 77405, 77541, 77584, 77897, 78521, 78539, 78540, 78825, 79155, 80090, 80263, 80425, 80584, 80911, 87870, 87976, 87991, 88069, 90646, 90804, 90811, 91117, 92333, 92390, 92408, 92984, 92988, 93079, 93166, 94119, 94200, 94279, 94354, 94364, 95371, 95443, 95889, 96502, 96960, 97235, 97328, 97459, 97546, 99565, 99797, 99933, 99985, 100105, 101419, 101921, 102068, 103005, 103130, 104539, 104687, 104713, 105257, 105291, 105855, 106019, 106224, 106254, 106475, 106633, 106973, 106996, 106998, 108070, 109152, 109448, 109472, 110224, 110266, 110633, 110655, 111081, 111403, 111888, 112078, 112673, 112767, 112825, 112836, 113307, 113346, 114134, 114674, 115147, 115194, 115707, 115811, 116279, 116316, 116386, 116405, 116828, 116948, 117074]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01_RESULTS.sha256 | Integrity checksum |
