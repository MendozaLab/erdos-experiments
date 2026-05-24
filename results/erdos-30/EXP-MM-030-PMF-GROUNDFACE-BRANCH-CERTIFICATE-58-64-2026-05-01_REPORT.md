# EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-64-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_GROUNDFACE_BRANCH_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Sources

- `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`: supplies the exact `n=58` exported face.
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01`: supplies exact exported faces for `59 <= n <= 64`.

## Result

The post-58 face does not collapse back to a single translated branch. Each subsequent row contains the previous exact face and the previous face shifted by `+1`, while new difference-skeleton branches accumulate rapidly. The field-selected Pareto face remains a small subset of a much larger branching ground face.

## Checks

- all source exports are `EXPORTED_ALL`: `True`
- all rows preserve field-selection split: `True`
- every row after `58` contains every previous witness exactly: `True`
- every row after `58` contains every previous witness shifted by `+1`: `True`
- skeleton count is strictly increasing from `58` to `64`: `True`

## Branch Window Table

| n | exact face | skeletons | field split | prev exact | prev +1 | Pareto count | joint count |
|---:|---:|---:|---|---:|---:|---:|---:|
| 58 | 10 | 2 | True | 0 | 0 | 2 | 1 |
| 59 | 18 | 4 | True | 10 | 10 | 2 | 1 |
| 60 | 54 | 18 | True | 18 | 18 | 1 | 1 |
| 61 | 152 | 49 | True | 54 | 54 | 4 | 4 |
| 62 | 398 | 123 | True | 152 | 152 | 6 | 6 |
| 63 | 1022 | 312 | True | 398 | 398 | 10 | 10 |
| 64 | 2360 | 669 | True | 1022 | 1022 | 11 | 11 |

## Claim Boundary

Finite exact-face branch certificate only. This does not prove Erdős #30, does not establish an asymptotic phase transition, and does not upgrade PMF analogy to theorem language.
