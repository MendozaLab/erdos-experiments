# Erdos #242 Salez Bound Sweep

Experiment: `EXP-MATH-ERDOS242-SALEZ-BOUND-SWEEP-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `FULL_HARD_STRIP_COVERAGE_WITHIN_SWEEP`

## Meaning

This sweep asks whether the bounded seven-equation search only works at a large constant window, or whether hard-strip coverage appears steadily as the window grows.

## Coverage By Bound

| Constant bound | Covered hard-strip primes | Coverage | First uncovered targets |
|---:|---:|---:|---|
| 5 | 1168 / 1181 | 0.988992 | `[24481, 26041, 29761, 33289, 48889, 53089, 66889, 70849]` |
| 10 | 1181 / 1181 | 1.000000 | `[]` |
| 20 | 1181 / 1181 | 1.000000 | `[]` |
| 40 | 1181 / 1181 | 1.000000 | `[]` |
| 80 | 1181 / 1181 | 1.000000 | `[]` |

## Boundary

The sweep is a tractability diagnostic, not a theorem. It is useful because it shows how quickly the Salez-family reconstruction covers the hard strip under finite constants.
