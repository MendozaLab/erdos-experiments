# Quorum Synthesis — EHP114 Interval Taylor / M14

Date: 2026-05-05  
Scope: Internal synthesis of GPT-5.5, Claude Opus 4.7, and Gemini 3.1 Pro
reviews of the n=14 interval Taylor / M14 packet.

## Consensus

All three model reviews agree on the main interpretation:

```text
axis endpoint budget passes;
uniform positive Taylor-matrix transport fails;
the failure is not just Gershgorin looseness;
the correct replacement target is radial reserve absorbing radial-base shape
softening on a local cone.
```

This is a shadow signature, not universal law. It does not settle Erdős #114
and does not establish local stability.

## Disagreements That Matter

1. **Theorem shape.**
   One review favored a weaker scalar theorem
   `D14 >= 12 eps^(1/14)` on an epsilon-scaled cone. Another kept the
   mixed-remainder absorption framing but warned that the old `+50000||s||^2`
   reserve may be too strong after softening. The synthesis is to test the
   scalar theorem first.

2. **Trust in `eps^(1/28)`.**
   The exponent is dimensionally natural because squaring it gives the radial
   Puiseux scale `eps^(1/14)`. But it is not yet proven. Treat it as a
   coordinate hypothesis, not doctrine.

3. **Admissibility contamination.**
   The ambient Taylor matrix used off-axis stencil points, and some were not
   root-admissible. Before canonizing the negative direction story, rerun the
   spectral scan with an admissible-only stencil.

4. **Critical path.**
   High-dimensional interval boxes are likely to overinflate. The next run
   should be a targeted admissible spectral scan plus an epsilon-scaled axis
   cone search before attempting a 24-dimensional box certificate.

## Concrete Next Runs

1. `EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01`
   - Find per-epsilon Taylor step sizes whose full diagonal/off-axis stencil is
     root-admissible.
   - Recompute the spectral lower bound.
   - Decide whether negative directions persist inside the admissible stencil.

2. `EXP-MATH-EHP114-N14-EPS-SCALED-CONE-M14-SEARCH-20260505-01`
   - Test `s = eps^(1/28) u` axis endpoints.
   - Search conservative `eta0` values.
   - Track the worst known direction `m6_sin_tangent`.

## Results Added After Quorum

The admissible-only spectral rerun landed:

```text
EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01
status = ADMISSIBLE_SHAPE_SOFTENING_CONFIRMED
all stencil points admissible = true
global interval spectral lower bound = -94465620.44867483
```

The epsilon-scaled axis search then showed:

```text
EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01
status = EPS_SCALED_AXIS_FAIL
all evaluated admissible axis points pass the scalar target
worst m6_sin_tangent lane certified through eta = 0.014
```

The `FAIL` is a conservative full-axis admissibility failure, not a scalar
deficit failure: many outward signed axis points leave the root-admissible
domain and are outside the theorem's `hadm` hypothesis.

The mixed spectral-direction search then showed:

```text
EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01
status = EPS_SCALED_SPECTRAL_DIRECTIONS_PASS
failure count among evaluated admissible points = 0
largest eta passing = 0.014
```

Updated synthesis: the next theorem should be the direct scalar statement
`D14(eps,s) >= 12 eps^(1/14)` on an epsilon-scaled admissible cone. The
stronger mixed-remainder absorption theorem remains a later route.

## Claim Ceiling

Safe internal language:

```text
The interval Taylor packet falsified uniform transported shape positivity and
redirected the n=14 local theorem target toward an epsilon-scaled reserve
absorption statement.
```

Unsafe language, paraphrased:

- Do not claim the conjecture is settled.
- Do not claim local stability is established.
- Do not claim `eps^(1/28)` is the true cone law before testing.
- Do not claim full-cone certification from axis or eigenvector sampling.
