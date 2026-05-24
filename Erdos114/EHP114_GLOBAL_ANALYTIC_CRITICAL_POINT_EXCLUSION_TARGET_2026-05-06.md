# EHP114 Global Analytic Critical-Point Exclusion Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` residual branch geometry only.

Claim ceiling: theorem-target packet only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Why The Local-Collar Route Is Retired

The branch-centered moving-frame pilot tested the Class A wall failures from L18. It processed 14 candidates, beginning with ownership key `2286:0`.

Result:

```text
validated center points           = 1 / 14
moving-frame tangent-zero count   = 1 / 14
quadratic center-strip passes     = 0 / 14
remaining Class A unresolved      = 14 / 14
```

The easiest family `2286:0` failed branch centering: the normal-axis endpoint intervals at `r = 0` had the same sign. So the intended branch center is not captured by the fixed local normal axis used in the residual box.

The only validated center was in source key `2484:1`, split path `root/y0/y1`, but its quadratic center-strip budget still failed:

```text
allowed C0                 = 0.0013565586804541244
quadratic center-strip C0  = 0.006771875136516755
margin                     = -0.00541531645606263
```

That retires the current local-collar sequence:

```text
axis collar -> normal collar -> third-order collar -> branch-centered moving frame
```

as a proof-facing closure route for the current residual representation.

## Active Theorem Target

The next useful theorem is not a sharper local box collar. It is a hard-cell critical-point exclusion theorem.

For the hard-cell root-affine family:

```text
F(x,y;u) = |p_u(x + i y)|^2 - 1
```

prove that on the unresolved residual domain:

```text
|grad F(x,y;u)| >= m_global > 0
```

or prove an equivalent decomposition:

```text
residual domain = excluded regions union certified regular collars
```

where the regular collars are defined by analytic geometry rather than inherited slab residual boxes.

## Required Ingredients

A proof-facing target must supply:

- an analytic partition of the residual hard-cell domain that does not depend on the failed local normal-axis centering;
- lower bounds for `|grad F|` or for a dominant directional derivative on each partition element;
- an exclusion rule for regions where `F` cannot cross zero;
- a branch-count/ownership rule that prevents double counting;
- a length bound that only applies after critical-point exclusion and wall separation are certified.

## Recommended Next Artifact

Create a theorem-target or diagnostic packet that classifies the 64 L18 residual regions by global critical-point obstruction:

```text
EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01
```

Minimum fields:

```text
processed_region_count
gradient_lower_bound_candidate
critical_candidate_region_count
regular_region_count
excluded_region_count
first_failed_condition
claim_ceiling
```

Acceptance for a useful next artifact:

```text
regular_region_count > 0
```

or an exact theorem target that identifies why `|grad F|` cannot be bounded with current constants.

## Non-Goals

Do not rerun all `4398` slab unresolved branches.

Do not run z64.

Do not promote length.

Do not expand Lean until the analytic critical-point theorem has a stable statement and real constants.
