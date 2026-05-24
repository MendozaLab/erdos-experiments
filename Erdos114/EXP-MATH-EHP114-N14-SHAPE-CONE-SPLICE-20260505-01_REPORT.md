# EHP114 n=14 Shape-Cone Splice Packet

Experiment: `EXP-MATH-EHP114-N14-SHAPE-CONE-SPLICE-20260505-01`

## Meaning

The radial/Puiseux lane is now certified for the radial family. This packet
asks the next question: does the n=14 nonradial shape cone look strong enough
to formulate a local-stability splice theorem?

The answer is yes as a target, not as a proof. This is a shadow signature, not
universal law. The shape numbers below are floating diagnostics; the next step
is interval hardening.

## Verdict

- Status: `SHAPE_SPLICE_TARGET_READY`
- Radial compact status: `COMPACT_MIDDLE_CERTIFIED`
- Radial tail status: `RADIAL_TAIL_CERTIFIED`
- Shape lambda target: `100000.0`
- Shape lambda below diagnostic bounds: `True`
- Current blocker: `Replace the floating n=14 tensor-cone matrix with an interval matrix certificate and prove mixed-remainder absorption.`

## Shape-Cone Diagnostic

| quantity | value |
|---|---:|
| quotient basis rank | 25 |
| expected rank `2n-3` | 25 |
| shape positive symmetric deficits | 24/24 |
| shape slopes below Hessian threshold 1.5 | 24/24 |
| mixed positive eigenvalues | 24/24 |
| floating min eigenvalue | 320878.725222 |
| Gershgorin lower bound | 308540.903818 |
| worst Gershgorin row | 5 |

## Mixed Remainder Obligation

For 0 < eps <= 1e-1 and ||s|| <= eta14(eps), prove R14(eps,s) <= 12*eps^(1/14) + 0.5*lambda14*||s||^2.

This is now the live mathematical bottleneck. The radial singularity is no
longer the first obstruction; the obstruction is proving that mixed terms
cannot eat more than half the certified radial deficit.

## Lean-Shaped Target

```lean
theorem ehp114_n14_radial_shape_remainder_splice
    (eps : Real) (s : ShapeQuotient 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hs : quotientNorm s <= eta14 eps) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      + (50000 : Real) * quotientNormSq s
      <= D14 (radialMode14 eps + shapeMode14 s) := by
  -- radial certificates + interval shape cone + mixed-remainder absorption
  sorry
```

## Source Boundary

No scorecard, D1, public document, Lean file, git, or email state was changed.
