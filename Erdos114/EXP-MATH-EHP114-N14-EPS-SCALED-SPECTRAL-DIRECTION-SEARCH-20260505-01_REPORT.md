# EHP114 n=14 Epsilon-Scaled Spectral-Direction Search

Experiment: `EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01`

## Meaning

This run follows the negative midpoint eigenvectors from the admissible Taylor matrices.
It tests the scalar target on mixed directions rather than only signed coordinate axes.

## Verdict

- Status: `EPS_SCALED_SPECTRAL_DIRECTIONS_PASS`
- Largest eta with all evaluated admissible spectral points passing: `0.014`
- Failure count among evaluated admissible points: `0`
- Oracle points evaluated: `60`

## Eta Grid

| eta | admissible/evaluated | inadmissible | all pass | min margin |
|---:|---:|---:|:---:|---:|
| 0.001 | 12/24 | 12 | yes | 10.4903287562 |
| 0.002 | 12/24 | 12 | yes | 10.4444521944 |
| 0.004 | 6/24 | 18 | yes | 12.0978128016 |
| 0.006 | 6/24 | 18 | yes | 12.0929803431 |
| 0.008 | 6/24 | 18 | yes | 12.0877403366 |
| 0.010 | 6/24 | 18 | yes | 12.0812770511 |
| 0.012 | 6/24 | 18 | yes | 12.0719074009 |
| 0.014 | 6/24 | 18 | yes | 12.0595154164 |

## Claim Ceiling

Finite mixed-direction evidence for the epsilon-scaled scalar theorem target at n=14. This is not a full cone certificate, not a Lean proof, and not a proof of Erdős #114.

## Next Blocker

Promote the passing admissible spectral-direction evidence to a certified low-dimensional cone or derive analytic Cauchy bounds that remove reliance on sampled eigenvectors.
