# EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-79-V2-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_GROUNDFACE_BRANCH_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Sources

- `EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-74-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-74-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-75-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-75-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-76-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-76-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-77-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-77-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-78-PAR8-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-78-PAR8-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-79-PAR8-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-79-PAR8-2026-05-01_RESULTS.json`

## Result

The post-72 face does not collapse back to a single translated branch. Each subsequent row contains the previous exact face and the previous face shifted by `+1`, while new difference-skeleton branches accumulate rapidly. The field-selected Pareto face remains a small subset of a much larger branching ground face.

## Checks

- all source exports are `EXPORTED_ALL`: `True`
- all rows preserve field-selection split: `True`
- every row after `72` contains every previous witness exactly: `True`
- every row after `72` contains every previous witness shifted by `+1`: `True`
- skeleton count is strictly increasing across the window: `False`
- skeleton count is nondecreasing across the window: `True`
- growth rule used for this certificate is nondecreasing with first post-anchor plateau allowed, then strict: `True`

## Branch Window Table

| n | exact face | skeletons | field split | prev exact | prev +1 | Pareto count | joint count |
|---:|---:|---:|---|---:|---:|---:|---:|
| 72 | 4 | 2 | True | 0 | 0 | 1 | 1 |
| 73 | 8 | 2 | True | 4 | 4 | 1 | 1 |
| 74 | 34 | 13 | True | 8 | 8 | 1 | 1 |
| 75 | 84 | 25 | True | 34 | 34 | 1 | 1 |
| 76 | 214 | 65 | True | 84 | 84 | 1 | 1 |
| 77 | 482 | 134 | True | 214 | 214 | 2 | 2 |
| 78 | 970 | 244 | True | 482 | 482 | 3 | 3 |
| 79 | 1974 | 502 | True | 970 | 970 | 8 | 8 |

## Claim Boundary

Finite exact-face branch certificate only. This does not prove Erdős #30, does not establish an asymptotic phase transition, and does not upgrade PMF analogy to theorem language.
