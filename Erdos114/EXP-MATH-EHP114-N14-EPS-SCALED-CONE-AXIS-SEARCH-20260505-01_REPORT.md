# EHP114 n=14 Epsilon-Scaled Cone Axis Search

Experiment: `EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01`

## Meaning

This run tests the weaker scalar target after the normalized change of variables
`s = eta * eps^(1/28) * u`. The exponent is treated as a coordinate
hypothesis, not a theorem.

## Verdict

- Status: `EPS_SCALED_AXIS_FAIL`
- Largest eta certified on every axis/sign/epsilon in the coarse grid: `None`
- Worst-mode largest eta certified: `0.014`
- Oracle points evaluated: `3100`

## Coarse Eta Grid

| eta | admissible/evaluated | inadmissible | all pass | min margin |
|---:|---:|---:|:---:|---:|
| 0.002 | 140/192 | 52 | yes | 6.90517600926 |
| 0.004 | 114/192 | 78 | yes | 7.06829532794 |
| 0.006 | 114/192 | 78 | yes | 7.16112481435 |
| 0.008 | 114/192 | 78 | yes | 7.05417894193 |
| 0.010 | 114/192 | 78 | yes | 6.77022625827 |
| 0.012 | 114/192 | 78 | yes | 6.78008233384 |
| 0.014 | 104/192 | 88 | yes | 6.76222134792 |
| 0.016 | 94/192 | 98 | yes | 6.9117454141 |
| 0.018 | 94/192 | 98 | yes | 6.48522366975 |
| 0.020 | 92/192 | 100 | yes | 8.29977449084 |

## Worst Mode

Mode: index `21`, label `m6_sin_tangent`, sign `-`.

This is the `m6_sin_tangent` danger lane identified by the interval Taylor packet.

## Claim Ceiling

Finite axis-direction evidence for the epsilon-scaled scalar theorem target at n=14. This is not a full cone certificate, not a Lean proof, and not a proof of Erdős #114.

## Next Blocker

Move from axis directions to low-dimensional spectral/admissible sub-cones, or prove analytic Cauchy bounds that dominate the observed radial-base shape softening.
