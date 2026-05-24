# EHP114 n=14 Local Mixed-Remainder Theorem Target

## Meaning

The n=14 local-stability problem has been narrowed to one missing theorem
family. The original strong form was:

```text
the mixed radial/shape remainder must be too small to cancel the certified
radial Puiseux deficit plus the certified positive shape-cone energy.
```

The 2026-05-05 Taylor/spectral runs show that this strong form is too
optimistic as the immediate target: the radial contraction creates genuine
shape softening. The current next theorem target is therefore weaker and more
direct:

```text
on an epsilon-scaled local admissible cone, the total deficit itself stays
above the radial Puiseux reserve floor.
```

This is a shadow signature, not universal law. It is a fixed-n local theorem
target, not a solution of Erdős #114.

## The Closure Algebra

Write the local deficit schematically as

```text
D14(eps, s) = R14(eps) + Q14(s) - M14(eps, s)
```

where:

- `R14(eps)` is the radial Puiseux deficit;
- `Q14(s)` is the positive shape-cone quadratic energy on the quotient space;
- `M14(eps, s)` is the mixed radial/shape remainder.

The available artifact-backed budgets are:

```text
R14(eps) >= 24 eps^(1/14)
Q14(s)   >= 100000 ||s||^2
```

The missing theorem should prove the conservative absorption bound:

```text
M14(eps, s) <= 12 eps^(1/14) + 50000 ||s||^2
```

Then the local deficit still has reserve:

```text
D14(eps, s) >= 12 eps^(1/14) + 50000 ||s||^2.
```

That last inequality is pure bookkeeping once the three component estimates are
available. The new Lean scratch file checks exactly this algebraic splice.

## Exact Local Domain

The theorem should not quantify over all admissible perturbations. The broad
admissible scout failed because unit-disk admissibility alone allows tangential
excursions outside the local quadratic regime.

The correct local domain is:

```text
0 < eps <= 1/10
RootsInClosedUnitDisk(radialMode14(eps) + shapeMode14(s))
||s|| <= min(eta14Boundary(eps), 0.008)
```

Here `eta14Boundary(eps)` is the analytic boundary radius that must eventually
replace the finite scout's empirical cap. The constant `0.008` is the current
certified local scout cap, not a universal constant.

## Lean-Shaped Targets

### Strong Later Target: Mixed-Remainder Absorption

The original strong target remains useful as a later route:

```lean
theorem ehp114_n14_local_mixed_remainder_absorption
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hlocal : quotientNorm s <= min (eta14Boundary eps) localConeCap14) :
    mixedRemainder14 eps s
      <= radialReserve14 * Real.rpow eps ((1 : Real) / 14)
        + shapeReserve14 * quotientNormSq s := by
  -- analytic theorem still open
  sorry
```

The scratch scaffold does not assert this theorem as proved. It packages it as
`MixedRemainderAbsorption14`, an explicit assumption that would close the local
deficit once supplied.

### Current Target: Epsilon-Scaled Scalar Deficit

The immediate target after the Taylor/spectral diagnosis is:

```lean
theorem ehp114_n14_eps_scaled_cone_deficit
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0_14 * Real.rpow eps ((1 : Real) / 28)) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= totalDeficit14 eps s := by
  -- direct total-deficit interval/Cauchy certificate target
  sorry
```

The exponent `1/28` is a coordinate hypothesis, not doctrine: squaring it lands
at the radial Puiseux scale `1/14`, but the true admissible cone law still has
to be earned by interval or analytic bounds.

## What Is Already Artifact-Backed

| component | artifact | status | role |
|---|---|---|---|
| radial theorem shape | `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01` | `CERTIFICATE_TARGET_READY` | fixes conservative constant `24` |
| radial compact middle | `EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01` | `COMPACT_MIDDLE_CERTIFIED` | covers `1e-4 <= eps <= 1e-1` |
| radial singular tail | `EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01` | `RADIAL_TAIL_CERTIFIED` | covers `0 < eps <= 1e-4` |
| shape matrix | `EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01` | `SHAPE_INTERVAL_MATRIX_CERTIFIED` | supports `lambda14 = 100000` |
| broad mixed scout | `EXP-MATH-EHP114-N14-ADMISSIBLE-MIXED-REMAINDER-SCOUT-20260505-01` | `FINITE_MIXED_SCOUT_FAIL` | shows admissibility alone is too broad |
| local mixed scout | `EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02` | `LOCAL_MIXED_SCOUT_PASS` | finite evidence for the local cap |

## What Remains Axiomatized / Open

The scratch scaffold isolates three propositions:

```lean
RadialCertificate14
ShapeConeCertificate14
MixedRemainderAbsorption14
```

The first two now have artifact-backed interval routes. They are still not
native Lean proofs of the analytic/numeric certification pipeline.

The third is the real missing theorem. It must be proved analytically or by a
formal interval remainder certificate. It cannot be replaced by the finite
scout.

## Attack Routes

## Route B First Run — 2026-05-05 Update

The first interval Taylor packet has now been executed:

```text
EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01
```

It produced a useful negative result:

```text
status = INTERVAL_TAYLOR_MATRIX_FAIL
axis endpoints pass the conservative M14 budget
ambient Taylor matrix fails the uniform positive-shape test
```

A follow-up spectral diagnosis confirms this was not merely Gershgorin
crudeness:

```text
EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01
status = RADIAL_BASE_SHAPE_SOFTENING_DETECTED
global interval spectral lower bound = -93160.96568063524
```

Interpretation: the boundary shape cone remains real, but it should not be
transported unchanged to radially contracted bases. The mixed term `M14` must
absorb **radial-base shape softening**. The correct target is therefore:

```text
radial reserve dominates the negative radial/shape curvature on the local cone
```

not:

```text
the shape cone remains uniformly positive after radial contraction
```

### Route B Second Update — Admissible Stencil + Epsilon-Scaled Search

The admissibility-contamination objection has now been tested:

```text
EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01
status = ADMISSIBLE_SHAPE_SOFTENING_CONFIRMED
all stencil points admissible = true
global interval spectral lower bound = -94465620.44867483
```

This means the negative spectral directions are not merely an artifact of
using non-root-admissible stencil points. The admissible Taylor step becomes
very small at small epsilon, so the matrix estimate is numerically harsh, but
the qualitative message survives: do not attack transported positive
shape-cone curvature as the next theorem.

The first epsilon-scaled scalar run was also executed:

```text
EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01
status = EPS_SCALED_AXIS_FAIL
all evaluated admissible axis points passed
failure reason = many outward axis points leave the root-admissible domain
worst m6_sin_tangent lane certified through eta = 0.014
```

The status is deliberately conservative: it fails the stronger "every signed
axis point is root-admissible" requirement. But under the theorem's actual
`hadm` hypothesis, every evaluated admissible axis point passed the scalar
reserve target.

Finally, the mixed-direction check followed the lowest admissible-stencil
midpoint eigenvectors:

```text
EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01
status = EPS_SCALED_SPECTRAL_DIRECTIONS_PASS
failure count among evaluated admissible spectral points = 0
largest eta with all evaluated admissible spectral points passing = 0.014
```

Interpretation: the scalar epsilon-scaled target is now the cleanest local
closure route. It has axis and spectral-direction evidence, but it is still not
a full cone certificate.

### Route A — Cauchy Bound On Local Derivatives

Treat the mixed term as the integral remainder after radial plus quadratic
shape expansion. Bound all third-order and radial/shape cross derivatives on
the local tube:

```text
0 < eps <= 1/10,
||s|| <= min(eta14Boundary(eps), 0.008).
```

This is the most direct route. It should produce explicit constants feeding the
`12` and `50000` budgets.

### Route B — Interval Taylor Model

Build a validated Taylor model over `(eps, s)` boxes:

```text
radial Puiseux model + shape quadratic model + interval remainder.
```

This is closer to the current Rust/inari and Arb stack. It should produce a
machine-verifiable certificate but may be harder to port into Lean. After the
first run, Route B should be reframed as a **softening absorption** model, not
as a proof that the boundary shape Hessian stays positive at every radial base.

### Route C — Jordan / Spectral Cone Normal Form

Diagonalize the quotient shape matrix into a positive spectral basis, then
bound mixed terms mode-by-mode. This is the cleanest articulation route, but it
should not be allowed to delay the Cauchy or interval route.

## Next Artifact For Another Agent

Attack this single proposition first:

```lean
ehp114_n14_eps_scaled_cone_deficit
```

Treat `MixedRemainderAbsorption14` as a stronger later theorem, not the
immediate blocker. The next agent should either:

1. promote the passing admissible axis/spectral-direction evidence to a
   low-dimensional certified cone, or
2. derive Cauchy-style analytic derivative bounds that imply the scalar
   epsilon-scaled theorem.

Do not reopen the radial theorem or shape-matrix certificate unless a defect is
found.

The older mixed proposition is:

```lean
MixedRemainderAbsorption14
```

It remains valuable, but it is no longer the first theorem to attack.

## Verification

Scratch file:

```text
erdos-experiments/Erdos30/scratch/Ehp114LocalMixedRemainderScratch.lean
```

Build command:

```bash
cd erdos-experiments/Erdos30
lake build Ehp114LocalMixedRemainderScratch
```

Expected meaning of a build pass: the Lean names, definitions, and closure
algebra are well-typed. It does not mean the mixed-remainder theorem has been
proved.

## Claim Ceiling

Safe internal language:

```text
The local n=14 #114 closure problem has been reduced to a single
mixed-remainder absorption theorem on a local admissible cone.
```

Unsafe language, paraphrased:

- Do not claim the Erdős #114 conjecture is settled.
- Do not claim the middle range is closed.
- Do not claim the local theorem is already established.
- Do not claim a Jordan-algebra proof.
