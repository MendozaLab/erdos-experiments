# EHP114 n=14 Exact-Length Lift Packet

Date: 2026-05-05

Scope: Subagent B bridge after the Python `SUBDIV8` low-dimensional cell
certificate. This packet works only from the current n=14 selected
`eps = 0.1` coefficient cell and does not claim an Erdos #114 solution.

Required phrase: shadow signature, not universal law.

## Meaning

The current low-dimensional PASS is real, but it is not yet proof-grade
lemniscate length. It certifies one selected coefficient cell for the existing
root-affine interval marching-squares oracle functional. The exact-length lift
has one load-bearing bottleneck:

```text
replace the marching-squares oracle functional by a validated exact-length
enclosure, or prove a conservative error bound connecting the two.
```

The good news is that the `SUBDIV8` pass has measurable slack. A future exact
or validated length engine does not have to reproduce the grid-oracle upper
bound exactly. It only has to keep the exact length upper bound below the
already implied cap:

```text
L_exact(C) <= Lstar_lower - 12 eps^(1/14)
```

for every one of the 64 subcells.

## Source State Read

Primary source artifacts:

- `erdos-experiments/Erdos114/low_dim_cone/EHP114_N14_LOW_DIM_CONE_CERTIFICATE_PACKET_2026-05-05.md`
- `erdos-experiments/Erdos114/cauchy_bounds/EHP114_N14_ANALYTIC_CAUCHY_BOUND_PACKET_2026-05-05.md`
- `erdos-experiments/Erdos114/EHP114_N14_ROOT_AFFINE_RUST_PORT_SPEC_2026-05-05.md`
- `erdos-experiments/Erdos114/low_dim_cone/EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.json`

Current low-dimensional result:

```text
experiment_id = EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01
status = ONE_CELL_ROOT_AFFINE_SUBDIV8_PASS
subcells = 64
failure_count = 0
Lstar_lower = 30.852910841548532
target = 12 * 0.1^(1/14) = 10.180114778928864
max marching-squares length upper = 18.110795101362747
min margin lower = 2.5620009612569206
```

Claim ceiling inherited from the source result:

```text
continuous grid-oracle cell certificate only; not exact lemniscate
certification, not a Lean theorem, and not a proof of Erdos #114.
```

## Budget Prototype Added

I added and ran:

```text
erdos-experiments/Erdos114/exact_length_lift/ehp114_exact_length_lift_budget.py
```

It produced:

- `EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.json`
- `EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_REPORT.md`
- `EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.sha256`

Status:

```text
BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE
```

This is a real derived artifact, but it is only a budget extraction. It does
not compute exact length and does not certify the marching-squares functional.

The uniform exact length cap on every subcell is:

```text
Lstar_lower - target = 20.672796062619668
```

The worst current subcell is:

```text
subcell = (6, 4)
u0_interval = [-0.00043749999999999995, -0.00021874999999999998]
u1_interval = [0.0008749999999999999, 0.0010937499999999999]
marching_length_upper = 18.110795101362747
allowed_additive_exact_length_error = 2.5620009612569206
allowed_relative_exact_length_error_vs_marching = 0.14146264407044962
```

So a uniform bridge of the form

```text
L_exact(C) <= L_ms(C) + E
```

would preserve the current scalar reserve on the full selected cell if

```text
E <= 2.5620009612569206.
```

Per-cell budgets are larger for many subcells; the mean additive budget is
`3.3689587217741392`.

## What Is Needed For Proof-Grade Exact Length

There are two viable lift routes.

### Route A: Direct Validated Exact-Length Enclosure

Replace each `length_upper` row by a proof-grade upper enclosure for the exact
lemniscate length over the root-affine subcell:

```text
for every subcell C:
  exact_length_upper(C) <= 20.672796062619668
```

or more sharply:

```text
exact_length_upper(C) <= marching_length_upper(C) + margin_lower(C).
```

Concrete requirements:

1. Keep the root-affine polynomial model:

   ```text
   p(z) = prod_i (z - r_i(a,b))
   ```

   with `a,b` intervalized on each `SUBDIV8` subcell.

2. Replace the grid contour length by a validated implicit-curve length
   enclosure for `F(z,a,b) = |p(z,a,b)|^2 - 1`.

3. For each curve patch, prove either:

   ```text
   length(patch) <= interval_arc_length_upper(patch)
   ```

   directly from interval implicit-function data, or use a collar/coarea
   overcount with explicit constants.

4. Emit the same audit rows as `SUBDIV8`, but with a field named
   `exact_length_upper` or `validated_length_upper`, not just
   `length_upper`.

Minimal data the exact validator must report per subcell:

```text
sub_i, sub_j,
u0_interval, u1_interval,
root_radius_upper,
gradient_lower_on_curve_or_collar,
curvature_or_hessian_upper,
validated_length_upper,
exact_length_cap,
margin_after_exact_length,
pass
```

### Route B: Conservative Marching-to-Exact Error Bound

Keep the current marching-squares upper functional, but prove a conservative
comparison theorem:

```text
L_exact(C) <= L_ms(C) + E(C)
```

Then require:

```text
E(C) <= margin_lower(C)
```

for all 64 subcells.

This route needs a theorem with real constants. The missing inputs are:

- a lower bound on `|grad(|p|)|`, equivalently `|p'|`, on the relevant
  `|p| = 1` curve or a collar around it;
- a validated collar width `tau` around `|p| = 1`;
- an interval upper bound on curvature/Hessian variation inside that collar;
- an explicit overcount bound showing that the grid oracle cannot undercount
  exact length by more than `E(C)`.

The coarea shape is:

```text
area({ ||p| - 1| <= tau }) =
  integral_{1-tau}^{1+tau} length({ |p| = t }) / |grad |p|| dt
```

To become useful here, it must be turned into a one-sided length upper bound,
not just a qualitative identity. The finite budget says the bound may be quite
coarse: a uniform additive error below `2.5620009612569206` already preserves
the selected-cell pass.

## Worst-Subcell Bridge Diagnostic

I ran the first Rust bridge diagnostic on the worst accepted subcell `(6,4)`.
It computes interval data for

```text
F(z,a,b) = |p(z,a,b)|^2 - 1
```

using the root-affine subcell, then tries to bound `|p'|` away from zero on
candidate level-set boxes. This is the first required input for an exact-length
bridge, because on `|p| = 1` the level-set gradient is controlled by `2|p'|`.

Artifact:

```text
erdos-experiments/scripts/erdos-114/bridge-diagnostic-worst-subcell-output-z16/
EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
```

Result:

```text
status = BRIDGE_REGULARITY_INTERVAL_UNRESOLVED
z-subdivision per active marching cell = 16
candidate level boxes = 95071
regularity unresolved boxes = 13560
min gradient lower candidate = 0.0
max Hessian upper candidate = 13614.580436750759
sum normal-drift error candidate = 149.62095139919995
sum relative-length error candidate = 423.1887593268329
available exact-length bridge budget = 2.5620009612530126
```

Interpretation: the bridge did not pass. The budget is generous, but the
current interval representation is too coarse. Even after splitting each active
marching cell into `16 x 16` `z`-boxes, the interval enclosure for `p'` still
contains zero on many boxes that may contain the level set. The simple
condition-ratio error candidates are also far over budget.

This failure is useful: it says the next exact-length proof cannot be obtained
by naive box interval arithmetic on `p'`. It needs one of:

1. Bernstein or affine arithmetic for `p'` on the level-set collar.
2. An implicit-function chart that follows the curve instead of boxing the
   ambient grid.
3. A coarea/Crofton-style overcount theorem with constants that avoid requiring
   pointwise `p'` separation on every ambient subbox.
4. Sharper interval root-affine arithmetic with subdivision driven by
   `F`-collar relevance rather than uniform `z` boxes.

The proof target remains unchanged, but the implementation route is now
clearer:

```lean
theorem ehp114_n14_subdiv8_marching_to_exact_length_error
    (C : Subdiv8Cell)
    (hRootAffine : RootAffineSubcell14 C)
    (hLevelChart :
      ValidatedLevelSetCharts C gammaC kappaC)
    (hGrid : MarchingSquaresOracleUpper C LmsC)
    (hBudget : chartContourError gammaC kappaC gridStepC <= marginLower C) :
    exactLemniscateLengthUpper C
      <= LmsC + chartContourError gammaC kappaC gridStepC := by
  sorry
```

The important change is `ValidatedLevelSetCharts`, not raw ambient
`GradientLowerOnLevelCollar`.

## Relation To The Analytic Cauchy Packet

The Cauchy packet is a broader analytic route. It tries to avoid sampled
eigenvectors and selected cells by proving a scaled shape-softening envelope:

```text
G14(eps,y) =
  eps^(-1/14) * (D14(eps, eps^(1/28)y) - D14(eps,0))
```

and then bounding:

```text
G14(eps,y) >= -12
```

on a fixed scaled `y`-ball.

That route remains open because `rhoY`, `B14`, the analytic tube, the zero
shape-gradient lemma, and the budget inequality are not yet certified. It is
cleaner for a full cone theorem. The exact-length route in this packet is
narrower: it tries to promote the already-passing selected `eps = 0.1`
SUBDIV8 cell from grid-oracle length to validated exact length.

## Smallest Next Theorem Target

The smallest target another agent can attack without re-reading the corpus is
the generic budget splice. It is intentionally independent of EHP geometry:

```lean
theorem exact_length_budget_splice
    (Lstar target Lms E Lexact : Real)
    (hBudget : Lms + E <= Lstar - target)
    (hBridge : Lexact <= Lms + E) :
    target <= Lstar - Lexact := by
  linarith
```

The EHP-specialized version should be packaged as:

```lean
theorem ehp114_n14_subdiv8_exact_length_budget_splice
    (C : Subdiv8Cell)
    (hCap :
      marchingLengthUpper C + exactLengthBridgeError C
        <= lstarLower14_eps01 - scalarTarget14_eps01)
    (hBridge :
      exactLemniscateLengthUpper C
        <= marchingLengthUpper C + exactLengthBridgeError C) :
    scalarTarget14_eps01
      <= lstarLower14_eps01 - exactLemniscateLengthUpper C := by
  -- pure arithmetic once the constants are in the row record
  linarith
```

That lemma does not solve the analytic problem, but it freezes the acceptance
criterion. The next analytic theorem is then sharply stated:

```lean
theorem ehp114_n14_subdiv8_marching_to_exact_length_error
    (C : Subdiv8Cell)
    (hRootAffine : RootAffineSubcell14 C)
    (hGrad : GradientLowerOnLevelCollar C gammaC)
    (hCurv : CurvatureUpperOnLevelCollar C kappaC)
    (hGrid : MarchingSquaresOracleUpper C LmsC)
    (hBudget : validatedContourError gammaC kappaC gridStepC <= marginLower C) :
    exactLemniscateLengthUpper C
      <= LmsC + validatedContourError gammaC kappaC gridStepC := by
  sorry
```

This is the real bottleneck theorem. Proving it for the worst subcell `(6,4)`
is the smallest meaningful analytic/exact-length lift after the current
Python certificate.

## Claim Ceiling And Bans

Safe statement:

```text
For one selected eps = 0.1 two-dimensional spectral-span coefficient cell, the
root-affine SUBDIV8 certificate has enough margin that a validated exact-length
enclosure may be allowed up to 2.5620009612569206 above the current worst
marching-squares upper value.
```

Unsafe statements:

- Do not say Erdos #114 is solved.
- Do not say exact lemniscate certification is done.
- Do not say the full cone is certified.
- Do not treat the Cauchy route as proved.
- Do not present this selected-cell result as a universal law.

## Verification Commands

Run from `/Users/kenbengoetxea/container-projects/apps/H2/Math`:

```bash
python3 erdos-experiments/Erdos114/exact_length_lift/ehp114_exact_length_lift_budget.py
```

Run from `/Users/kenbengoetxea/container-projects/apps/H2/Math`:

```bash
(cd erdos-experiments/Erdos114/exact_length_lift && shasum -a 256 -c EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.sha256)
```
