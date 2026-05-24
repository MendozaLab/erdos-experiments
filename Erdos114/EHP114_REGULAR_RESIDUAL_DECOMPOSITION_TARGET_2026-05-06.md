# EHP114 Regular Residual Decomposition Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` residual-region decomposition after the global critical-point exclusion diagnostic.

Claim ceiling: theorem-target packet only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Input Artifact

Source diagnostic:

```text
EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01
```

Result summary:

```text
processed residual regions       = 64
regular regions                  = 8
critical-candidate regions       = 56
excluded regions                 = 0
minimum regular gradient bound   = 33.07361176255703
accepted length upper            = 20.316451752723314
hard-cell cap                    = 20.672796062619668
margin                           = 0.3563443098963539
```

The first critical-candidate obstruction remains:

```text
source key   = 2484:1
split path   = root/y0/y0
reason       = both Fx and Fy intervals contain zero under root-affine uncertainty
```

## Meaning

The local-collar route failed because it tried to make one inherited chart carry the residual geometry. The global critical-point diagnostic changes the proof question: eight residual regions already have a robust dominant derivative in the x-direction, while fifty-six regions remain critical candidates under the full root-affine uncertainty.

That split is useful. The proof-facing route should no longer ask whether every residual tube admits the same collar. It should build a residual-domain atlas:

```text
regular slice:        x-chart theorem with ownership and length
critical slice:       separate critical-candidate exclusion/decomposition theorem
```

## Interpretation Layer

The Hessian picture is still useful as search geometry. It explains why naive quadratic positivity and fixed-axis collars keep failing near the soft directions.

The MDL / Hilbert-Riemann view is useful as architecture: regular regions are locally compressible by a short chart, while critical candidates are exactly where that local description loses compression. But this remains an organizing lens, not a certificate. The certificate must still be an explicit interval or analytic inequality.

## Next Proof-Facing Target

First prove the regular slice, because it now has real constants.

For each of the eight regular regions, prove:

```text
Fx has a fixed sign and |Fx| >= m_x > 0
F crosses the owned collar exactly once
all seams are half-open and ownership-stable
length <= width_y * sqrt(1 + sup(|Fy/Fx|)^2)
```

The proof packet should report:

```text
regular_region_count
regular_regions_certified
ownership_duplicate_count
regular_slice_length_upper
total_validated_length_upper
margin_to_cap
first_failed_condition
claim_ceiling
```

Pass condition for the next diagnostic:

```text
regular_regions_certified = 8
ownership_duplicate_count = 0
total_validated_length_upper <= 20.672796062619668
```

No critical-candidate region should be promoted through this regular-slice theorem.

## Critical-Candidate Follow-Up

The fifty-six critical candidates need a different theorem. The first obstruction is not length; it is uncertainty in both partial derivatives. The next candidate theorem should either:

- shrink root-affine uncertainty analytically enough to recover a directional derivative; or
- prove a no-critical-point statement directly from the polynomial/root geometry; or
- decompose each critical candidate into excluded regions plus smaller regular collars with a new ownership rule.

Do not treat failure of the eight-region regular theorem as failure of the critical-candidate theorem. They are separate proof obligations.

## Non-Goals

Do not rerun z64.

Do not rerun the global slab validator.

Do not process all `4398` unresolved slab branches until the eight-region regular theorem is working.

Do not expand Lean until the regular-slice theorem has stable constants and a complete artifact.
