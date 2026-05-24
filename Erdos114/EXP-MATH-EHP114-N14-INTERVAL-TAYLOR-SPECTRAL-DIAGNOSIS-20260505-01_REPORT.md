# EHP114 n=14 Interval Taylor Spectral Diagnosis

Experiment: `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01`

## Meaning

The first interval Taylor packet failed the Gershgorin test. This diagnosis
checks whether that was merely a crude row-sum bound. It recomputes the
interval Taylor matrices and applies:

```text
lambda_min(midpoint matrix) - Frobenius(interval radius).
```

## Verdict

- Status: `RADIAL_BASE_SHAPE_SOFTENING_DETECTED`
- Global interval spectral lower bound: `-93160.96568063524`
- All positive by interval spectral bound: `False`
- All pass working lambda14: `False`

## Per-Epsilon Diagnosis

| eps | midpoint lambda_min | radius | spectral lower | Gershgorin lower | lambda14 pass |
|---:|---:|---:|---:|---:|---:|
| 1e-04 | -39575.5873226 | 8.11e-07 | -39575.5873234 | -74273.966431 | no |
| 1e-03 | -93160.96568 | 6.77e-07 | -93160.9656806 | -131824.721059 | no |
| 1e-02 | -32480.9669225 | 4.79e-07 | -32480.9669229 | -42486.401566 | no |
| 1e-01 | -450.924794927 | 2.56e-07 | -450.924795183 | -589.086951476 | no |

## Consequence

The matrix failure is not just a Gershgorin artifact. The shape Hessian around
radially contracted bases develops negative directions under this Taylor
stencil. The local axis endpoints still pass the M14 budget, so the route is
not dead; it means the mixed-remainder theorem must explicitly absorb
radial-base shape softening.

The next theorem target should be phrased as:

```text
radial reserve dominates negative radial/shape curvature on the local cone.
```

not as:

```text
the shape cone remains uniformly positive after radial contraction.
```

## Claim Ceiling

This is a shadow signature, not universal law. It is a diagnostic of the
fixed-n interval Taylor strategy, not a proof or disproof of EHP #114.
