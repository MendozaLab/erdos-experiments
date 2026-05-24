# EHP114 Shape Cone Jordan Articulation

## Purpose

This note records the current mathematical interpretation of the n=14 shape
lane after the radial and shape-matrix interval artifacts. It is book-facing
guidance and theorem-target guidance. It is not a public claim, scorecard
upgrade, or Lean status change.

The right phrase is:

```text
shadow signature, not universal law
```

## Meaning

The local model now splits into two different kinds of geometry.

The radial direction is singular. Its deficit is not quadratic; it follows a
Puiseux law with exponent `1/14`. That is the boundary term.

The nonradial shape directions are spectral. After quotienting the trivial
rotation/scale directions, the interval shape matrix has a certified positive
Gershgorin lower bound. That is the cone term.

Together, the picture is:

```text
Puiseux boundary term + positive spectral shape cone + mixed remainder.
```

## Jordan-Style Reading

The Jordan-algebra language is potentially the right articulation layer because
the shape side behaves like a spectral cone:

- a self-adjoint quadratic operator on the quotient shape space;
- positive spectral directions after removing symmetry modes;
- an energy functional that behaves like a trace form;
- a cone of admissible positive directions bounded by a local radius.

The honest claim is narrower than "Jordan algebra solves EHP." The current
artifact-backed claim is only that the n=14 shape matrix looks naturally
expressible as a Jordan-style spectral cone, while the radial direction sits on
a Puiseux boundary outside ordinary Hessian geometry.

## Candidate Algebraic Models

The most plausible models are:

1. **Euclidean Jordan algebra of real symmetric matrices.**
   Use this if the certified shape matrix can be represented as a positive
   self-adjoint operator with spectral idempotents.

2. **Spin-factor-like cone.**
   Use this if the quotient modes reduce to one scalar radial coordinate plus a
   Euclidean vector of shape modes with Lorentz-cone-style positivity.

3. **Fourier block spectral cone.**
   Use this if the shape modes remain block-diagonal by Fourier frequency and
   the cone is better described as a direct product of positive blocks rather
   than one simple algebra.

The current data favors the third model as the immediate technical target and
the first two as articulation candidates.

## Theorem Target

The next useful theorem is not "build a Jordan algebra." It is a local spectral
cone theorem that can later be recognized as Jordan-like:

```lean
theorem ehp114_n14_shape_operator_positive_on_quotient
    (s : ShapeQuotient 14) :
    lambda14 * quotientNormSq s <= shapeQuadratic14 s := by
  -- certified interval matrix plus quotient basis normalization
  sorry
```

Then it must splice into the mixed-remainder theorem:

```lean
theorem ehp114_n14_local_mixed_remainder_absorption
    (eps : Real) (s : ShapeQuotient 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk (radialMode14 eps + shapeMode14 s))
    (hlocal : quotientNorm s <= min (eta14Boundary eps) (0.008 : Real)) :
    mixedRemainder14 eps s
      <= (12 : Real) * Real.rpow eps ((1 : Real) / 14)
        + (50000 : Real) * quotientNormSq s := by
  sorry
```

## Claim Ceiling

Safe internal language:

```text
The n=14 shape lane now has an interval-certified positive spectral-cone
signature, and the Jordan language is a plausible articulation of that cone.
The live mathematical blocker remains the uniform mixed-remainder splice.
```

Do not claim:

- a general EHP result;
- a Jordan-algebra closure theorem;
- any upgrade for n > 14;
- any public-facing advance claim beyond the local artifact stack.

## Artifacts This Note Depends On

- `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01`
- `EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01`
- `EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01`
- `EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01`
- `EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02`
