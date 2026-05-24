# EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 69 through n = 69 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 75000 |
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
| 69 | 10 | 218433851 | 66412 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 69 | 10 | 66412 | 11.103633 | 0 | 0 | NA | 66412 | 0 | 0.001481 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 69 | `[0, 3, 11, 20, 33, 43, 62, 64, 68, 69]` | `[0, 1, 5, 28, 43, 49, 57, 60, 67, 69]` | `[0, 3, 11, 21, 33, 49, 62, 64, 68, 69]` | 0.001481 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 69 | EXPORTED_ALL | 66412 | 66412 | SKIP_DISTANCE_EDGES_TOO_LARGE | 2205243666 | true | `[16569, 16849, 17019, 17021, 17228, 17230, 17338, 17357, 17435, 17461, 17472, 17569, 17624, 17639, 17641, 17688, 17696, 17778, 20154, 20381, 20599, 20615, 20738, 20775, 20784, 20994, 20995, 21009, 21319, 22868, 22879, 22882, 22883, 22892, 23052, 23067, 23090, 23203, 23229, 23238, 23257, 23352, 23491, 23568, 23574, 25305, 25335, 25345, 25352, 25504, 25556, 25653, 25886, 25887, 25893, 26628, 26655, 26810, 27137, 27626, 27739, 27990, 28030, 28036, 28039, 28049, 28114, 28444, 28470, 28524, 28645, 28681, 28697, 28699, 28704, 28908, 28983, 29281, 29283, 29360, 29372, 29406, 29418, 29660, 29668, 29780, 29813, 29847, 36915, 37067, 37224, 37231, 37258, 37338, 37477, 37582, 37664, 37717, 37781, 37782, 39590, 39752, 39834, 39844, 39866, 39869, 39872, 39990, 40015, 40026, 40305, 40321, 40337, 40347, 41583, 41601, 41611, 41722, 41753, 41755, 42034, 42039, 42079, 42219, 42220, 42233, 42268, 42321, 42355, 43213, 43222, 43230, 43473, 43534, 43599, 43708, 43735, 44318, 44352, 44551, 44663, 45418, 45545, 45711, 45805, 45825, 45945, 45982, 46137, 46235, 46458, 46476, 46568, 46625, 46877, 48899, 48909, 49002, 49212, 49302, 49312, 49322, 49323, 49328, 49347, 50804, 50818, 50839, 50851, 51032, 51040, 51088, 51146, 51173, 51177, 52238, 52242, 52247, 52436, 52450, 52461, 52646, 52767, 52781, 53345, 53353, 53644, 53703, 53883, 54329, 54596, 54646, 55049, 55073, 55121, 55157, 55598, 55612, 55613, 55686, 55780, 55868, 55994, 56093, 56218, 57409, 57411, 57421, 57497, 57540, 57589, 57606, 57748, 58331, 58541, 58571, 58618, 58758, 59203, 59208, 59214, 59261, 59413, 59517, 60076, 60418, 60741, 60830, 60865, 61010, 61118, 61246, 61262, 62020, 62394, 62686, 62893, 63047, 63226, 63263, 63346, 63557, 63592, 63680, 63694, 63761, 63792, 63836, 63839, 63915, 64430, 64575, 64672, 64919, 65180, 65275, 65313, 65318, 65465, 65479, 65552, 65788, 65869, 65962, 66050, 66122, 66171]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01_RESULTS.sha256 | Integrity checksum |
