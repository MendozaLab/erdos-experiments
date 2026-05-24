# EHP114 n=14 Analytic Cauchy-Bound Packet

Date: 2026-05-05

Scope: internal fixed-n analytic route for the current n=14 local theorem
target. This packet does not certify the theorem. The analytic route is not yet
a proof.

Phrase to preserve: shadow signature, not universal law.

## Meaning

The current target is no longer transported positivity of the shape Hessian.
The admissible spectral Taylor packet shows that radial contraction creates
real shape softening, even when the Taylor stencil points are root-admissible.
The clean analytic route is therefore to bound the total loss from shape
motion after the epsilon scaling

```text
s = eps^(1/28) y.
```

The reason for the exponent is bookkeeping: a quadratic shape displacement at
scale `eps^(1/28)` naturally lands at the same scale as the radial Puiseux
reserve, `eps^(1/14)`. The theorem still has to be earned by analytic bounds.
The exponent is a coordinate hypothesis here, not a law.

## Current Theorem Target

The downstream target remains the scalar deficit theorem:

```lean
theorem ehp114_n14_eps_scaled_cone_deficit
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0_14 * Real.rpow eps ((1 : Real) / 28)) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= totalDeficit14 eps s := by
  -- open analytic/Cauchy theorem target
  sorry
```

The Cauchy route should prove this by decomposing the total deficit into a
radial Puiseux reserve plus a bounded scaled shape softening:

```text
D14(eps, s)
  = D14(eps, 0)
    + [D14(eps, s) - D14(eps, 0)]

D14(eps, 0) >= 24 eps^(1/14)
D14(eps, s) - D14(eps, 0) >= -12 eps^(1/14)
```

The first line is the radial reserve. The second line is the analytic Cauchy
closure target.

## Lean-Shaped Cauchy Closure

Introduce the scaled shape function

```text
G14(eps, y)
  = eps^(-1/14) * (totalDeficit14 eps (eps^(1/28) y)
      - totalDeficit14 eps 0).
```

The route asks for a uniform analytic envelope on a fixed scaled `y`-ball,
not a sampled eigenvector certificate.

```lean
def ScaledShape14 := ShapeQuotient14

axiom scaledShapeMode14 :
    Real -> ScaledShape14 -> ShapeQuotient14

axiom scaledShapeMode14_norm :
    forall (eps : Real) (y : ScaledShape14),
      quotientNorm (scaledShapeMode14 eps y)
        = Real.rpow eps ((1 : Real) / 28) * quotientNorm y

axiom scaledSoftening14 :
    Real -> ScaledShape14 -> Real

def CauchyEnvelope14
    (rhoY B14 : Real) : Prop :=
  0 < rhoY /\
  0 <= B14 /\
  forall eps : Real,
    0 < eps -> eps <= (1 : Real) / 10 ->
      -- analytic extension of y |-> scaledSoftening14 eps y
      -- to the complexified quotient ball of radius rhoY
      ComplexAnalyticOnScaledBall14 eps rhoY /\
      SupNormOnScaledBall14 eps rhoY <= B14 /\
      scaledSoftening14 eps 0 = 0 /\
      HasZeroShapeGradient14 eps

def CauchySofteningBudget14
    (rhoY B14 eta0 : Real) : Prop :=
  0 < eta0 /\
  eta0 < rhoY /\
  B14 * (eta0 / rhoY)^2 / (1 - eta0 / rhoY) <= 12

theorem ehp114_n14_eps_scaled_cone_deficit_from_cauchy
    (rhoY B14 eta0 : Real)
    (hRadial :
      forall eps : Real, 0 < eps -> eps <= (1 : Real) / 10 ->
        (24 : Real) * Real.rpow eps ((1 : Real) / 14)
          <= totalDeficit14 eps 0)
    (hCauchy : CauchyEnvelope14 rhoY B14)
    (hBudget : CauchySofteningBudget14 rhoY B14 eta0)
    (hEta : eta0_14 <= eta0)
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0_14 * Real.rpow eps ((1 : Real) / 28)) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= totalDeficit14 eps s := by
  -- proof shape:
  -- 1. write s = eps^(1/28) y with quotientNorm y <= eta0
  -- 2. apply Cauchy to bound scaledSoftening14 eps y >= -12
  -- 3. multiply by eps^(1/14)
  -- 4. splice with radial reserve 24 -> target reserve 12
  sorry
```

The analytic lemma inside `CauchyEnvelope14` can be stated in a cleaner generic
form:

```lean
theorem cauchy_softening_from_zero_gradient
    (f : ScaledShape14 -> Real)
    (rhoY B eta : Real)
    (hanalytic : ComplexAnalyticOnScaledBall f rhoY)
    (hsup : SupNormOnScaledBall f rhoY <= B)
    (hzero : f 0 = 0)
    (hgrad : HasZeroGradient f 0)
    (heta : 0 <= eta)
    (hlt : eta < rhoY) :
    forall y : ScaledShape14, quotientNorm y <= eta ->
      f y >= -B * (eta / rhoY)^2 / (1 - eta / rhoY) := by
  -- finite-dimensional Cauchy/Taylor estimate
  sorry
```

That generic lemma is the best next Lean target because it is independent of
EHP geometry. The EHP-specific work is then reduced to proving analyticity,
the sup envelope, and the zero-gradient normalization.

## Constants That Must Be Bounded

These constants should be made explicit before any certification claim:

| constant | meaning | current status |
|---|---|---|
| `alpha = 1/14` | radial Puiseux exponent | artifact-backed target |
| `beta = 1/28` | shape scaling exponent | coordinate hypothesis |
| `C_rad = 24` | radial reserve lower bound | artifact-backed interval target |
| `C_target = 12` | retained scalar reserve | theorem target |
| `eta0_14` | epsilon-scaled cone radius | finite evidence supports testing `0.014`, not a proof |
| `rhoY` | scaled complex Cauchy radius in `y` coordinates | open |
| `B14` | normalized sup bound for scaled shape softening on `||y|| <= rhoY` | open |
| `m` | derivative order needed in Taylor/Cauchy expansion | at least order 2; order 3 if separating Hessian plus remainder |
| `r_eps` | radial epsilon contour or annulus radius | open if differentiating in `eps`; avoidable if eps is fixed |
| `r_coeff` | coefficient-space analytic tube radius around the radial base | open |
| `r_root` | root-admissibility boundary margin to keep roots in the closed unit disk | open |
| `r_len` | contour separation for the lemniscate length integral | open |
| `N_shape = 25` | quotient coordinate dimension in scratch model | scratch-backed, final quotient definition still needed |

The budget inequality is the concrete numerical gate:

```text
B14 * (eta0 / rhoY)^2 / (1 - eta0 / rhoY) <= 12.
```

If this gate is met with `eta0 >= eta0_14`, the Cauchy envelope implies the
scalar theorem target without using sampled eigenvectors.

## Artifact-Backed Facts

The following are backed by current local artifacts:

1. The current live theorem target is the epsilon-scaled scalar deficit
   theorem:
   `D14(eps,s) >= 12 eps^(1/14)` under root admissibility and
   `||s|| <= eta0_14 eps^(1/28)`.

2. The phrase to preserve is: shadow signature, not universal law.

3. The radial packet supplies an interval target with working constant
   `C14 = 24` for the radial Puiseux reserve.

4. The shape matrix packet supports a positive boundary shape-cone constant
   in the uncontracted shape slice, recorded in the scratch algebra as
   `shapeLambda14 = 100000`.

5. The admissible-stencil spectral Taylor run confirms shape softening after
   radial contraction:

   ```text
   status = ADMISSIBLE_SHAPE_SOFTENING_CONFIRMED
   all stencil points admissible = true
   global interval spectral lower bound = -94465620.44867483
   ```

6. The epsilon-scaled axis search is finite evidence only. It did not certify
   every signed axis point because many outward points leave the admissible
   domain, but every evaluated admissible axis point passed the scalar target.
   The danger lane was `m6_sin_tangent`, with worst-mode evidence through
   `eta = 0.014`.

7. The epsilon-scaled spectral-direction search is finite evidence only. It
   followed the lowest midpoint eigenvectors and found zero failures among
   evaluated admissible points through `eta = 0.014`.

8. The Lean scratch file only checks theorem shapes and algebraic splicing.
   It intentionally leaves the analytic theorem as an axiom.

## Current Assumptions

These are not yet artifact-backed as proofs or certificates:

1. `totalDeficit14` has the needed complex analytic extension in shape
   variables after radial scaling.

2. The scaled normalized function `G14(eps,y)` has a uniform sup bound `B14`
   on a fixed scaled `y`-ball `||y|| <= rhoY` for every
   `0 < eps <= 1/10`.

3. The scaled shape gradient vanishes at `y = 0`. This is plausible from
   symmetry/criticality of the radial mode, but it must be stated and proved
   as a lemma.

4. A root-admissible real point `s` inside the theorem cone can be lifted to a
   scaled coordinate `y` without leaving the analytic tube needed by the
   Cauchy contour.

5. The constants can satisfy
   `B14 * (eta0 / rhoY)^2 / (1 - eta0 / rhoY) <= 12` with
   `eta0 >= eta0_14`.

6. The quotient coordinate model in the scratch file matches the final
   analytic quotient by scale/rotation modes.

## Feasible Route

Prove a scaled Cauchy envelope for the total deficit directly:

1. Define `G14(eps,y)` as the normalized shape softening.
2. Prove `G14(eps,0) = 0`.
3. Prove `d_y G14(eps,0) = 0` from the radial mode's first variation.
4. Prove a uniform complex analytic extension for `||y|| <= rhoY`.
5. Bound the sup norm by `B14` on that ball.
6. Apply the generic Cauchy lemma to obtain
   `G14(eps,y) >= -12` for `||y|| <= eta0`.
7. Splice with the radial reserve `24 eps^(1/14)`.

De-risking experiment or lemma:

```text
EXP-MATH-EHP114-N14-SCALED-CAUCHY-ENVELOPE-SCAN-20260505-01
```

This should estimate candidate pairs `(rhoY, B14)` by sampling the complexified
coefficient tube, not real eigenvector directions. The matching Lean lemma is:

```lean
theorem ehp114_n14_scaled_shape_first_variation_zero
    (eps : Real) (hpos : 0 < eps) (hsmall : eps <= (1 : Real) / 10) :
    HasZeroShapeGradient14 eps := by
  sorry
```

If the first-variation lemma fails, the Cauchy budget must include a linear
term and the `eps^(1/28)` cone is probably too wide.

## Risky Route

Try to recover a signed Hessian-plus-remainder theorem:

```text
D14(eps,s) - D14(eps,0)
  >= -K2(eps) ||s||^2 - K3(eps) ||s||^3
```

and then show, after `s = eps^(1/28)y`, that

```text
K2(eps) eta0^2 eps^(1/14)
  + K3(eps) eta0^3 eps^(3/28)
  <= 12 eps^(1/14).
```

This route is risky because the admissible spectral Taylor packet already
found very negative Hessian lower bounds after radial contraction. It may
still work if those Hessian constants scale in a way that is harmless on the
root-admissible cone, but the existing evidence says not to make transported
positive curvature the premise.

De-risking experiment or lemma:

```text
EXP-MATH-EHP114-N14-SCALED-HESSIAN-REMAINDER-BOUND-20260505-01
```

This should bound the scaled quantities

```text
eps^(-1/14) * lambda_min(shape Hessian at eps)
eps^(-1/14) * third_order_remainder(eps, eta eps^(1/28))
```

over admissible boxes. The Lean-side lemma would be:

```lean
theorem ehp114_n14_scaled_hessian_softening_budget
    (eps : Real) (hpos : 0 < eps) (hsmall : eps <= (1 : Real) / 10)
    (s : ShapeQuotient14)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0_14 * Real.rpow eps ((1 : Real) / 28)) :
    totalDeficit14 eps s
      >= totalDeficit14 eps 0
        - (12 : Real) * Real.rpow eps ((1 : Real) / 14) := by
  sorry
```

## Next Theorem For Another Agent

Attack this without re-reading the whole packet stack:

```lean
theorem ehp114_n14_scaled_total_deficit_softening_cauchy
    (rhoY B14 eta0 : Real)
    (hCauchy : CauchyEnvelope14 rhoY B14)
    (hBudget : CauchySofteningBudget14 rhoY B14 eta0)
    (hEta : eta0_14 <= eta0)
    (eps : Real) (y : ScaledShape14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hy : quotientNorm y <= eta0) :
    totalDeficit14 eps (scaledShapeMode14 eps y)
      >= totalDeficit14 eps 0
        - (12 : Real) * Real.rpow eps ((1 : Real) / 14) := by
  sorry
```

Then splice it with the radial reserve:

```lean
theorem ehp114_n14_eps_scaled_cone_deficit_from_scaled_softening
    (hRadial :
      forall eps : Real, 0 < eps -> eps <= (1 : Real) / 10 ->
        (24 : Real) * Real.rpow eps ((1 : Real) / 14)
          <= totalDeficit14 eps 0)
    (hSoft :
      forall eps y, 0 < eps -> eps <= (1 : Real) / 10 ->
        quotientNorm y <= eta0_14 ->
          totalDeficit14 eps (scaledShapeMode14 eps y)
            >= totalDeficit14 eps 0
              - (12 : Real) * Real.rpow eps ((1 : Real) / 14)) :
    EpsScaledDeficit14 := by
  sorry
```

This is the shortest analytic bridge from the current evidence to the current
Lean-shaped scalar target.

## Guardrails

Do not claim a full cone certificate from the axis or spectral-direction
searches. Do not claim local stability from this packet. Do not treat the
Cauchy route as certified until `rhoY`, `B14`, the analytic tube, and the
budget inequality are all bound.

## Source Artifacts Read

- `erdos-experiments/Erdos114/EHP114_N14_LOCAL_MIXED_REMAINDER_THEOREM_TARGET_2026-05-05.md`
- `erdos-experiments/Erdos114/PACKET_README_EHP114_N14_EPS_SCALED_CONE_2026-05-05.md`
- `erdos-experiments/Erdos114/THIRD_PARTY_REVIEW_EHP114_INTERVAL_TAYLOR_M14_2026-05-05.md`
- `erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01_REPORT.md`
- `erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01_REPORT.md`
- `erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01_REPORT.md`
- `erdos-experiments/Erdos30/scratch/Ehp114LocalMixedRemainderScratch.lean`
- `erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean`

