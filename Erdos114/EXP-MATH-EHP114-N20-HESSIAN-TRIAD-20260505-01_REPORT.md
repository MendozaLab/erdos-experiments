# EXP-MATH-EHP114-N20-HESSIAN-TRIAD-20260505-01 Report

## Status

Diagnostic n=20 Hessian/Puiseux triad for EHP114. This is not a proof and not
a publication packet. The claim ceiling remains: this is a shadow signature,
not universal law.

## Verdict

- Run status: `COMPLETE`
- Triad verdict: `SUPPORTS_STRATIFIED_SHADOW_SIGNATURE`
- Component decision: radial, Fourier, and tensor-cone machinery were runnable
  now by reusing existing local probe functions through a scoped wrapper.

The result preserves the stratified reading. The radial axis follows the
Puiseux exponent target, while the nonradial Fourier and tensor-cone probes
retain broad positive finite-difference shape signals under floating
estimators. This supports the next interval-hardening target, but it does not
establish a local or global EHP114 theorem.

## Component Summary

| component | runnable now? | main signal | diagnostic pass? |
|---|---|---|---|
| radial hypergeometric | yes | tail slope 0.0499976270055 vs target 0.05; relative error 4.74599e-05 | True |
| Fourier quotient Hessian | yes | positive modes 37/37; rank 37 vs expected 37 | True |
| tensor-cone mixed shape | yes | mixed positive eigenvalues 36/36; condition proxy 2.56314 | True |

## Radial Puiseux

The n=20 radial family lands on the expected `1/20` Puiseux exponent within the
predeclared 5 percent tolerance. The smallest tested eps is
`1e-08`, and the tail proxy
`deficit / eps^(1/20)` is `41.5978372626`.
This says the smooth radial Hessian model is still the wrong object at n=20.

## Nonradial Shape

The Fourier basis has the expected reduced rank `2n - 3 = 37`. Under the
floating marching-squares estimator, `37` of
`37` tested quotient modes are positive at all eps values. The
large reference-error flag remains `None`, so this
is shape evidence only, not an interval certificate.

The tensor-cone mixed matrix excludes the singular radial direction. Its shape
basis rank is `36`, and the finite-difference proxy has
minimum eigenvalue `601487.934073` and maximum eigenvalue
`1541696.72864`. The condition proxy is compared against the
predeclared diagnostic bound `100.0`.

## Lean-Shaped Next Target

```lean
theorem ehp114_n20_stratified_shadow_local_bound
    (r : Real) (s : ShapeQuotient 20)
    (hr_pos : 0 < abs r) (hr_small : abs r <= delta20)
    (hs_small : quotientNorm s <= eta20) :
    C20 * Real.rpow (abs r) ((1 : Real) / 20)
      + lambda20 * quotientNormSq s
      <= D20 (radialMode20 r + shapeMode20 s) + R20 r s := by
  -- interval Puiseux radial bound plus quotient shape-cone remainder bound
  sorry
```

This target is deliberately local and stratified: radial Puiseux first, then
quotient shape-cone positivity plus a mixed remainder bound. The missing piece
is interval control of the remainders, not another floating sweep.

## Source Boundary

All durable outputs from this run were written under:

`/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114`

No scorecard, D1, public docs, git, email, CLAUDE.md, or AGENTS.md was changed.
