# EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-74-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_RESET_WINDOW_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Sources

- `EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-74-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-74-2026-05-01_RESULTS.json`

## Result

The `h = 11` reset window is now visible through `n = 74`.

At `n = 72`, the face collapses to `4` exact maximizers. At `n = 73`, the face is exactly the `n = 72` face plus its `+1` translate, with no new witnesses. At `n = 74`, the same inherited-plus-translated rule still holds, but new branch production restarts: `22` additional witnesses appear and the skeleton count jumps from `2` to `13`.

## Checks

- source exports complete: `True`
- source parity matched: `True`
- `h` stays `11`: `True`
- every post-72 row contains the previous face unchanged: `True`
- every post-72 row contains the previous face shifted by `+1`: `True`

## Reset Window Table

| n | h(n) | exact maximizers | skeletons | prev exact | prev +1 | new witnesses | Pareto count | distance edges |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 72 | 11 | 4 | 2 | - | - | - | 1 | 6 |
| 73 | 11 | 8 | 2 | 4 | 4 | 0 | 1 | 28 |
| 74 | 11 | 34 | 13 | 8 | 8 | 22 | 1 | 561 |

## Claim Boundary

Finite reset-window certificate only. This does not prove Erdős #30, does not establish an asymptotic theorem, and does not upgrade PMF analogy to theorem language.
