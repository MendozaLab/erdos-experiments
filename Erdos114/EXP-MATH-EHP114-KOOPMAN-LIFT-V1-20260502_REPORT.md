# EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502

## Honest Scope

This is a prototype scoping artifact, not a certified proof of EHP at any n. The Koopman lift here uses EDMD with finite observable basis; the resulting spectral gap is an empirical estimate, not an interval-arithmetic certificate. To translate this into a closure-relevant artifact, the spectral gap would need to be (1) derived analytically rather than empirically estimated, (2) shown n-invariant or growing in n, (3) lifted to interval-arithmetic certificates. None of those are achieved in this scoping run.

## Setup

- Degree: n = 3
- State slice (sym-reduced): (0.0, 0.0, -1.0) = z^3 - 1
- Snapshot grid: 7^3 = 343 states in [-0.4, 0.4]^3 around x_star
- Observable basis size: K = 11
- Euler step h = 0.1

## Koopman Spectrum

- L\* (closed form):       9.1797242223
- L\* (radial marcher):    9.1713318637  (rel. err = 9.14e-04)
- Constant-function eigenvalue distance to 1: 4.15e-11
- Dominant non-trivial |lambda|: 1.000019
- Spectral gap in unit disk:     -0.000019
- Implied continuous decay rate: -0.0002
- # eigvals inside unit disk:    7
- # eigvals outside unit disk:   3

## Step-Size Sweep

| h (flow time) | |lambda_d|_max | gap_in_disk |
|---:|---:|---:|
| 0.100 | 1.000019 | -0.000019 |
| 0.300 | 1.000102 | -0.000102 |
| 1.000 | 1.000000 | -0.000000 |
| 3.000 | 1.000000 | -0.000000 |
| 10.000 | 1.000001 | -0.000001 |
| 30.000 | 1.000004 | -0.000004 |

## Margin Probe

- Sample-grid margin (this prototype): 4.53%
- v5 preprint certified margin (n=3):   6.10%

## Notes

- See KOOPMAN_LIFT_DESIGN_2026-05-02.md for the design rationale.
- This run does NOT modify the v5 preprint or any erdos-114 results.
- This run does NOT publish externally; Cooley filter applies.
