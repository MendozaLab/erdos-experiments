# EHP114 n=14 Admissible-Stencil Spectral Taylor Scan

Experiment: `EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01`

## Meaning

The previous Taylor matrix used an ambient stencil with some non-admissible
off-axis points. This scan shrinks the Taylor step per epsilon until every
diagonal and off-axis stencil point is root-admissible.

## Verdict

- Status: `ADMISSIBLE_SHAPE_SOFTENING_CONFIRMED`
- All stencil points admissible: `True`
- Max root radius: `1.000000000001`
- Global interval spectral lower bound: `-94465620.44867483`
- All positive by interval spectral bound: `False`

## Per-Epsilon Scan

| eps | admissible h | spectral lower | Gershgorin lower | positive |
|---:|---:|---:|---:|---:|
| 1e-04 | 9.44955e-06 | -94465620.4487 | -114144044.361 | no |
| 1e-03 | 9.4535e-05 | -1417389.76757 | -1810787.93178 | no |
| 1e-02 | 0.000949327 | -18629.7112196 | -20536.3168723 | no |
| 1e-01 | 0.004 | -450.924795183 | -589.086951476 | no |

## Consequence

If the lower bound remains negative here, radial-base shape softening survives
the admissibility correction. The next target is the epsilon-scaled total
deficit theorem, not a transported positive shape-cone theorem.

## Claim Ceiling

This is a diagnostic of a fixed-n local route. It does not settle Erdős #114 and
does not prove local stability.
