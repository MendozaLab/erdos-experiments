# EHP114 p-Prime Root-Location Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` critical-candidate regions after direct `p'(z)` exclusion.

Claim ceiling: theorem-target packet only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Current Certified Local State

The residual work has narrowed in three steps:

```text
L21 residual split                 = 8 regular regions + 56 critical candidates
L24 monotone Taylor exclusions     = 8 / 8 regular regions closed
L27 p-prime exclusions, p=8 tiles  = 8 / 56 critical candidates closed
remaining critical candidates      = 48
```

The local length budget remains unchanged:

```text
accepted length upper              = 20.316451752723314
hard-cell cap                      = 20.672796062619668
margin                             = 0.3563443098963539
```

## Meaning

The key analytic invariant is now explicit:

```text
F(z) = |p(z)|^2 - 1
grad F(z) = 0 on F=0 only if p'(z) = 0
```

because `|p(z)|=1` implies `p(z) != 0`.

So the remaining critical-candidate problem should not be carried by the wider product `p'(z) * conjugate(p(z))`. The next proof-facing object is the location of roots of `p'`.

## What The Latest Artifact Shows

The direct affine `p'` diagnostic at parameter subdivision 8 closed exactly two ownership families:

```text
2484:1 / root/y*
3707:1 / root/y*
```

Those eight regions are now regular by `p'(z) != 0` on every parameter tile. The first remaining obstruction is:

```text
2587:3 / root/y0/y0
```

where affine `p'` intervals still contain the origin.

## Next Theorem Target

Build a derivative-root location theorem for the remaining 48 regions.

For the hard-cell root-affine family, prove one of:

```text
A. every zero of p'_u lies outside each remaining residual box;
B. each remaining residual box is separated from all p'_u zeros by a certified distance;
C. the residual box decomposes into smaller regions where p' avoids zero or F avoids zero.
```

The recommended first diagnostic is:

```text
EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01
```

Minimum result fields:

```text
processed_remaining_candidate_count
pprime_root_excluded_count
pprime_root_near_count
still_unresolved_count
minimum_root_box_distance
first_failed_condition
total_validated_length_upper
exact_length_cap
margin_to_cap
claim_ceiling
```

Useful pass condition:

```text
pprime_root_excluded_count > 0
```

Strict closure condition:

```text
still_unresolved_count = 0
total_validated_length_upper <= 20.672796062619668
```

## Non-Goals

Do not rerun z64.

Do not rerun global slab length.

Do not treat parameter tiling beyond `p=8` as the next proof idea unless the root-location theorem identifies a specific tile-level obstruction.

Do not promote the hard-cell certificate until the remaining 48 critical candidates are closed.
