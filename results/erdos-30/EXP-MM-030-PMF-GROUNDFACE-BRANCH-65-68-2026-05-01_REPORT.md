# EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 65 through n = 68 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 50000 |
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
| 65 | 10 | 105668869 | 5018 |
| 66 | 10 | 127010338 | 9994 |
| 67 | 10 | 152419320 | 19418 |
| 68 | 10 | 182586477 | 36234 |

## Parity Gate

Checked 4 n-values. h(n) matched in 4. Maximizer counts matched in 4. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 65 | 10 | 5018 | 8.520787 | 0 | 0 | NA | 5018 | 0 | 0.001608 |
| 66 | 10 | 9994 | 9.209740 | 0 | 0 | NA | 9994 | 0 | 0.000000 |
| 67 | 10 | 19418 | 9.873956 | 0 | 0 | NA | 19418 | 0 | 0.001542 |
| 68 | 10 | 36234 | 10.497753 | 0 | 0 | NA | 36234 | 0 | 0.000000 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 65 | `[0, 1, 9, 21, 25, 35, 52, 58, 63, 65]` | `[0, 1, 8, 20, 36, 50, 54, 60, 63, 65]` | `[0, 1, 15, 25, 28, 47, 54, 59, 63, 65]` | 0.001608 |
| 66 | `[0, 2, 10, 19, 35, 42, 48, 62, 63, 66]` | `[0, 1, 10, 34, 39, 45, 47, 59, 62, 66]` | `[0, 2, 15, 22, 27, 48, 56, 62, 65, 66]` | 0.000000 |
| 67 | `[0, 2, 10, 19, 39, 40, 45, 52, 63, 67]` | `[0, 1, 8, 29, 42, 46, 52, 61, 64, 66]` | `[0, 2, 13, 23, 29, 47, 59, 62, 66, 67]` | 0.001542 |
| 68 | `[0, 3, 11, 20, 27, 49, 50, 62, 64, 68]` | `[0, 1, 8, 26, 36, 48, 59, 63, 65, 68]` | `[0, 3, 12, 26, 28, 50, 57, 63, 67, 68]` | 0.000000 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 65 | EXPORTED_ALL | 5018 | 5018 | SKIP_DISTANCE_EDGES_TOO_LARGE | 12587653 | true | `[708, 737, 1284, 1306, 1900, 1945, 2206, 2216, 2301, 2337, 2446, 2491, 2602, 2610, 2641, 3235, 3242, 3282, 3462, 3513, 3740, 3744, 3751, 3770, 3829, 3869, 3882, 3887, 3891, 3938, 3976, 4091, 4171, 4283, 4393, 4402, 4498, 4502, 4532, 4544, 4566, 4610, 4611, 4687, 4739, 4818, 4844]` |
| 66 | EXPORTED_ALL | 9994 | 9994 | SKIP_DISTANCE_EDGES_TOO_LARGE | 49935021 | true | `[2301, 2326, 2996, 3038, 3127, 3133, 3499, 3548, 3586, 3637, 4067, 4400, 4552, 4710, 4733, 4872, 5701, 6688, 6937, 7479, 7516, 7550, 7585, 7888, 7952, 7966, 7970, 8474, 8479, 8808, 8945, 8946, 9142, 9250, 9284, 9354, 9379, 9385, 9638, 9642]` |
| 67 | EXPORTED_ALL | 19418 | 19418 | SKIP_DISTANCE_EDGES_TOO_LARGE | 188519653 | true | `[4107, 4176, 4297, 5391, 5527, 5538, 5651, 5763, 5765, 5772, 5793, 5795, 6528, 7513, 7517, 7529, 7531, 7597, 7664, 7671, 8013, 8126, 8178, 8210, 8248, 8517, 8532, 8541, 8587, 8598, 8858, 8895, 9028, 9206, 9250, 10729, 11528, 11641, 11705, 11775, 12277, 12341, 12465, 12528, 12860, 12940, 12960, 13004, 13025, 13034, 13065, 13093, 13401, 13746, 13749, 13764, 13828, 13844, 13996, 14053, 14111, 14141, 14169, 14173, 14199, 14204, 14213, 14238, 14286, 15049, 15058, 15122, 15478, 15566, 15592, 15669, 15952, 15974, 16011, 16022, 16061, 16268, 16353, 16362, 16373, 16397, 16487, 16769, 16775, 16780, 16969, 17271, 17406, 17419, 17589, 17659, 17691, 17696, 17712, 17717, 17825, 17842, 18002, 18055, 18066, 18116, 18279, 18290, 18345, 18486, 18529, 18634, 19108, 19137, 19250, 19273, 19340, 19376, 19411, 19416]` |
| 68 | EXPORTED_ALL | 36234 | 36234 | SKIP_DISTANCE_EDGES_TOO_LARGE | 656433261 | true | `[9698, 9834, 9844, 10160, 10187, 11554, 11567, 11692, 11802, 11865, 11916, 13053, 13291, 13310, 13465, 13494, 13496, 13525, 13536, 14579, 14645, 15032, 15292, 15298, 15324, 15357, 15592, 15648, 16008, 16058, 16112, 16209, 16251, 16477, 16680, 20757, 20991, 20996, 21043, 21059, 22182, 22347, 22406, 23521, 23653, 24227, 24299, 24362, 24365, 25081, 25566, 25696, 26088, 27256, 27475, 28410, 28535, 28547, 29211, 29242, 29682, 29777, 29805, 29815, 29840, 29849, 29955, 30211, 30309, 30550, 31019, 31870, 31933, 32284, 32439, 32470, 32763, 32872, 33154, 33185, 33543, 33587, 34073, 34163, 34247, 34473, 34692, 34693, 34810, 34861, 34867, 35038, 35041, 35042, 35363, 35376, 35461, 35572, 35590, 35702, 35741, 35822, 35877, 35975, 36032, 36042, 36067, 36105, 36187]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01_RESULTS.sha256 | Integrity checksum |
