# EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-73-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_RESET_WINDOW_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Sources

- `EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01_RESULTS.json`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01`: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01_RESULTS.json`

## Result

The reset at `n = 72` opens the `h = 11` layer with only `4` exact maximizers. At `n = 73`, the same `h = 11` layer doubles to `8` exact maximizers.

The doubling is exact persistence plus translated persistence: every `n = 72` witness appears unchanged at `n = 73`, every `n = 72` witness shifted by `+1` appears at `n = 73`, and there are no other `n = 73` witnesses.

## Checks

- source exports complete: `True`
- source parity matched: `True`
- `h` stays `11`: `True`
- `n = 73` contains the `n = 72` face unchanged: `True`
- `n = 73` contains the `n = 72` face shifted by `+1`: `True`
- `n = 73` has no additional witnesses: `True`
- skeleton count stays `2`: `True`

## Reset Window Table

| n | h(n) | exact maximizers | skeletons | Pareto count | distance edges |
|---:|---:|---:|---:|---:|---:|
| 72 | 11 | 4 | 2 | 1 | 6 |
| 73 | 11 | 8 | 2 | 1 | 28 |

## Claim Boundary

Finite reset-window certificate only. This does not prove Erdős #30, does not establish an asymptotic theorem, and does not upgrade PMF analogy to theorem language.
