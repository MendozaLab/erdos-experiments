# EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01 — PMF Transfer-Operator Parity Scan

## Identification

| Field | Value |
|---|---|
| Experiment ID | EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01 |
| Erdős Problem | #30 — finite Sidon set rigidity |
| Data integrity | REAL_COMPUTATION — exact transfer-state enumeration, no sampling |
| Scan window | n = 71 through n = 71 |
| Frontier k | 10 |
| Prune deficiency | 0 |
| Ground-face export cap | 225000 |
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
| 71 | 10 | 311297246 | 203840 |

## Parity Gate

Checked 1 n-values. h(n) matched in 1. Maximizer counts matched in 1. Mismatches: [].

## Spectral / Ground-State Summary

| n | h(n) | degeneracy | entropy ln | h-1 count | h-2 count | gap | terminal states | peak states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 71 | 10 | 203840 | 12.225091 | 0 | 0 | NA | 203840 | 0 | 0.001424 |

## Zero-Temperature Field Tilts

After the cardinality field selects the exact h(n) ground-state layer, the top-k frontiers are the zero-temperature response to small prefix and mass fields.

| n | prefix-field witness | mass-field witness | joint-field witness | joint score |
|---|---|---|---|---:|
| 71 | `[0, 4, 13, 21, 40, 45, 56, 68, 70, 71]` | `[0, 1, 6, 31, 44, 53, 55, 63, 67, 70]` | `[0, 4, 13, 23, 34, 51, 63, 65, 66, 71]` | 0.001424 |

## Small Ground-Face Export

| n | export status | exact maximizers | exported | edge status | edge total | field split | pareto minima |
|---|---|---:|---:|---|---:|---|---|
| 71 | EXPORTED_ALL | 203840 | 203840 | SKIP_DISTANCE_EDGES_TOO_LARGE | 20775270880 | true | `[57047, 57160, 57198, 57505, 57621, 57635, 57943, 57986, 58057, 58230, 58259, 58532, 58697, 58741, 58859, 58910, 59015, 59354, 59382, 64570, 64664, 64730, 65007, 65077, 65080, 65086, 65514, 65520, 66010, 66109, 66117, 66134, 66242, 66314, 66317, 66324, 66328, 66378, 66383, 70580, 70686, 70714, 71003, 71306, 71331, 71335, 71611, 71640, 71647, 71652, 71672, 71791, 71839, 71940, 71976, 71987, 72082, 72177, 72178, 72402, 72404, 75330, 75556, 75660, 75797, 75805, 76264, 76266, 76471, 78530, 78538, 78774, 78838, 79015, 79199, 79227, 79261, 79270, 79449, 79474, 79492, 80780, 80806, 80830, 80839, 80857, 81013, 81042, 81050, 81176, 81280, 81281, 82451, 82452, 82520, 82664, 82697, 82798, 82858, 82922, 82940, 83218, 83265, 83312, 83694, 83709, 83712, 83721, 83742, 83858, 84144, 84170, 84562, 84569, 84639, 84643, 84873, 84884, 84885, 85074, 85105, 85455, 85651, 85653, 85774, 85972, 86072, 86596, 114481, 114528, 114865, 114914, 114916, 114940, 114942, 114965, 115195, 115335, 115346, 115557, 115636, 115738, 115766, 115770, 115794, 115880, 115886, 115899, 116012, 116038, 116072, 120431, 120937, 121000, 121184, 121247, 121362, 121420, 121456, 121546, 121643, 121832, 121851, 121860, 121947, 122017, 124894, 124942, 124974, 124980, 125405, 125416, 125425, 125441, 125574, 125585, 125610, 125611, 125756, 125781, 126178, 128664, 128742, 128753, 129426, 129546, 129704, 131173, 131174, 131195, 131329, 131788, 131922, 131960, 132900, 133130, 133420, 133478, 133479, 133511, 133606, 134450, 134538, 134661, 134817, 135371, 135396, 135591, 136236, 136263, 136378, 136825, 137147, 149342, 149387, 149638, 149649, 149656, 149895, 149994, 150006, 150010, 150021, 150025, 150075, 150201, 150231, 150384, 150579, 150632, 150689, 150737, 154024, 154188, 154197, 154252, 154262, 154456, 154507, 154630, 154887, 154970, 154975, 154995, 155051, 157598, 157784, 157806, 157808, 157839, 157843, 158187, 158200, 158445, 158528, 158578, 158735, 158745, 158792, 160299, 160508, 160662, 160772, 160871, 160994, 161031, 162531, 162740, 163329, 164188, 164198, 164299, 164358, 164380, 164446, 164558, 165035, 165549, 165869, 166037, 166184, 166403, 166421, 166530, 166564, 166604, 166626, 166678, 166802, 167176, 167290, 170852, 170895, 171120, 171259, 171263, 171265, 171363, 171464, 171707, 171739, 171811, 174346, 174621, 174663, 174690, 174855, 174863, 174904, 174992, 175005, 175064, 175078, 175107, 175175, 175238, 176948, 177268, 177282, 177429, 177665, 179168, 179264, 179478, 179562, 179686, 179702, 179723, 179782, 179811, 179827, 180604, 180629, 180647, 180812, 180993, 181021, 181844, 181859, 181879, 182034, 182180, 182218, 182242, 182526, 182528, 182549, 182630, 182729, 182780, 182789, 183124, 183317, 183504, 183597, 183605, 183612, 183971, 184021, 184210, 184278, 186296, 186319, 186433, 186534, 186667, 186688, 186753, 186772, 188342, 188428, 188524, 188668, 188704, 188747, 188766, 188775, 188776, 188795, 189979, 189982, 190003, 190071, 190098, 190771, 190796, 191172, 191208, 191738, 191936, 192007, 192356, 192434, 192703, 192863, 192991, 193285, 194709, 194854, 194958, 194962, 194966, 195069, 195088, 195747, 195801, 195845, 195960, 196033, 196074, 196810, 196868, 196982, 197280, 197444, 197528, 197860, 198126, 198128, 198203, 198612, 198614, 198644, 199242, 199245, 199310, 199372, 199444, 199535, 200006, 200030, 200112, 200329, 200366, 200515, 200666, 200887, 200938, 201219, 201430, 201715, 202239, 202244, 202584, 202859, 202865, 203029, 203066, 203198, 203454, 203693, 203826]` |

## Interpretation

This is a parity engine, not yet a compressed transfer matrix. The point is to prove that the PMF state representation can reproduce the exact ground-state surface before we trust spectral language.

The theorem language remains blocked: these packets are finite parity and near-ground evidence, not a Sidon proof.

## Artifacts

| File | Type |
|---|---|
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01_RESULTS.json | Structured transfer-operator results |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01_REPORT.md | Human-readable report |
| EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01_RESULTS.sha256 | Integrity checksum |
