# EXP-MM-030-PMF-AFFINE-EMBEDDING-N57-2026-04-30

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_AFFINE_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Derived from `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`.

## Result

At `n = 57`, all six exported exact maximizers share one positive-difference skeleton. The exposed mass winner and joint winner are not different skeletons: the joint winner is the mass winner translated by `+1`.

`[2, 4, 16, 23, 31, 34, 47, 51, 56, 57] = [1, 3, 15, 22, 30, 33, 46, 50, 55, 56] + 1`

## Local Affine Relations

| relation | from | to | parameter |
|---|---:|---:|---:|
| translation | 0 | 2 | 1 |
| translation | 0 | 4 | 2 |
| reflection_about_n | 0 | 5 | 57 |
| translation | 1 | 3 | 1 |
| translation | 1 | 5 | 2 |
| reflection_about_n | 1 | 4 | 57 |
| translation | 2 | 0 | -1 |
| translation | 2 | 4 | 1 |
| reflection_about_n | 2 | 3 | 57 |
| translation | 3 | 1 | -1 |
| translation | 3 | 5 | 1 |
| reflection_about_n | 3 | 2 | 57 |
| translation | 4 | 0 | -2 |
| translation | 4 | 2 | -1 |
| reflection_about_n | 4 | 1 | 57 |
| translation | 5 | 1 | -2 |
| translation | 5 | 3 | -1 |
| reflection_about_n | 5 | 0 | 57 |

## Translation Chains

- Chain A: `0 -> 2 -> 4` by successive `+1` translations.
- Chain B: `1 -> 3 -> 5` by successive `+1` translations.
- The exposed handoff occurs on Chain B: mass selects `3`, joint selects `5`.

## Certificate Checks

- all witnesses share one difference skeleton: `true`
- mass-to-joint is +1 translation: `true`
- Chain A present: `true`
- Chain B present: `true`
- reflection pairs present about `n = 57`: `true`

## Interpretation

The `n = 57` handoff is an embedding-selection event inside one difference skeleton. Maxwell's Litmus points to missing face geometry; the affine certificate says that geometry is a translated embedding chain, not a new difference-memory skeleton.

## Claim Boundary

Finite affine-embedding certificate only. This does not prove the Sidon asymptotic, does not prove a phase transition, and does not upgrade PMF to proof language.
