# EHP114 Critical-Candidate Theorem Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` residual regions after the regular-slice monotone Taylor exclusion.

Claim ceiling: theorem-target packet only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Current Certified Local State

The L21 global critical-point diagnostic split the 64 L18 residual regions into:

```text
regular x-chart regions      = 8
critical candidates          = 56
excluded by L21              = 0
```

The L24 monotone Taylor wall diagnostic then closed the eight regular regions:

```text
processed regular regions    = 8
monotone Taylor exclusions   = 8
ownership duplicates         = 0
wall Taylor failures         = 0
new branch length promoted   = 0
```

The length budget remains unchanged:

```text
accepted length upper        = 20.316451752723314
hard-cell cap                = 20.672796062619668
margin                       = 0.3563443098963539
```

## Meaning

The residual problem is now cleaner. The easy regular slice was not a hidden branch contribution; it was an over-wide residual bucket that a dependency-aware wall bound could exclude.

The remaining work is not length. It is the 56-region critical-candidate theorem:

```text
both Fx and Fy still contain zero under full root-affine uncertainty
```

That can mean one of two things:

1. the interval representation is too dependent and must be sharpened; or
2. a genuine near-critical geometry needs a different analytic invariant.

The artifact does not distinguish these yet.

## Next Proof-Facing Theorem

For the 56 critical candidates, prove one of the following on each region:

```text
A. no critical point exists: |grad F| >= m > 0 after dependency-aware root/coordinate decomposition;
B. F avoids zero after affine/Taylor root-parameter cancellation;
C. the region decomposes into smaller owned regular collars plus excluded boxes.
```

The recommended first diagnostic is an affine root-parameter gradient test:

```text
EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01
```

Minimum result fields:

```text
processed_critical_candidate_count
gradient_affine_regular_count
affine_excluded_count
still_critical_candidate_count
ownership_duplicate_count
total_validated_length_upper
exact_length_cap
margin_to_cap
first_failed_condition
claim_ceiling
```

Useful pass condition:

```text
gradient_affine_regular_count + affine_excluded_count > 0
```

Strict closure condition:

```text
still_critical_candidate_count = 0
total_validated_length_upper <= 20.672796062619668
```

## Non-Goals

Do not rerun z64.

Do not rerun global slab length.

Do not process all `4398` slab unresolved branches until the 56-region critical-candidate theorem has a working local mechanism.

Do not promote the hard-cell certificate until the 56 critical candidates are closed.
