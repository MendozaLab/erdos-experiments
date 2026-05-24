# EHP114 Regular-Slice Wall Taylor Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` regular-slice wall separation after the L22 regular residual diagnostic.

Claim ceiling: theorem-target packet only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Input Artifact

Source diagnostic:

```text
EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01
```

Result summary:

```text
processed regular regions      = 8
certified regular regions      = 0
derivative sign failures       = 0
wall separation failures       = 8
ownership duplicates           = 0
accepted length upper          = 20.316451752723314
hard-cell cap                  = 20.672796062619668
margin                         = 0.3563443098963539
```

The first failed condition is:

```text
2286:0:root/y0/y0 failed because x-wall F intervals do not have strict opposite signs
```

## Meaning

The good news is that the regular-slice derivative condition is stable. All eight regions keep a negative `Fx` interval, so the geometry is locally a regular `x=f(y)` graph candidate.

The bad news is that the inherited x-walls are not proof walls. Direct interval evaluation of `F` on the two x-walls is still too wide; both wall intervals include zero. That is a dependency problem in the wall certificate, not a length-budget failure and not a critical-point failure.

## Next Theorem Target

Replace raw x-wall interval evaluation with a Taylor or affine wall certificate.

For each of the eight regular regions, choose a center point or center strip and prove:

```text
F(x_wall, y, u) = F(center, y0, u0)
               + Fx * dx
               + Fy * dy
               + Fu * du
               + controlled second-order remainder
```

Then use the already-certified `Fx < 0` lower bound to build two artificial proof walls:

```text
x_left_proof(y,u)  < root_x(y,u) < x_right_proof(y,u)
```

with strict signs after dependency-aware cancellation.

The target is not to widen the inherited slab blindly. The target is to produce a dependency-aware bracket whose walls are tied to the actual root location.

## Required Result Fields For The Next Artifact

```text
processed_regular_region_count
wall_taylor_certified_count
wall_taylor_fail_count
dependency_width_reduction
regular_slice_length_upper
total_validated_length_upper
exact_length_cap
margin_to_cap
first_failed_condition
claim_ceiling
```

Pass condition:

```text
wall_taylor_certified_count = 8
total_validated_length_upper <= 20.672796062619668
```

## Non-Goals

Do not process the 56 L21 critical-candidate regions in this artifact.

Do not rerun z64.

Do not rerun the global slab validator.

Do not treat L22 wall failure as a failure of regularity; it is a wall-certificate failure.
