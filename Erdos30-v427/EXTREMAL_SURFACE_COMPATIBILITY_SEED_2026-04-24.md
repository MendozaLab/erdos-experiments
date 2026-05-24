# Extremal-Surface Compatibility Seed for Erdos #30

**Date:** 2026-04-24
**Primary packet:** `EXP-MM-030-COMPATIBILITY-SUMMARY-MAXIMIZER-DIAGNOSTICS-2026-04-24`
**Stress packet:** `EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24`
**Beyond-50 packet:** `EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`
**Rust extension packet:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24`
**Rust continuation packet:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24`
**Rust boundary packet:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24`
**Rust instrumented packet:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-69-2026-04-24`
**Rust second-pinch packet:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24`
**Scope:** Finite exact evidence, not an asymptotic theorem

## Idea

The current compatibility signal appears to live on the exact extremal surface.
It is not yet a generic dense-Sidon stability phenomenon.

The April 24 compatibility-summary packet shows a strong one-sided pattern on
exact maximizers: the density-adjusted-mass-best witness usually keeps prefix
residual essentially zero, while the prefix-best witness usually pays visible
density-adjusted mass cost.

The near-maximizer pilot then tested the first obvious robustness question by
adding the `|A| = h(n)-1` layer. That layer does not preserve the pattern. So the
right theorem seed should mention the exact extremal surface explicitly.

The beyond-50 probe tested the other obvious question: whether the exact-surface
signal immediately collapses when the original `10 <= n <= 50` window is
extended. It does not. For `51 <= n <= 55`, the one-sided normalized comparison
holds in `5/5` values, and the joint witness is mass-best in `4/5`. That is
still finite evidence, but it moves the story from "interesting packet artifact"
to "small exact-surface ridge worth stress-testing."

The Rust extension `56 <= n <= 60` keeps the ridge but makes it less naive. The
direct normalized comparison holds in `4/5` values, with the exception at
`n = 57`. That is exactly the kind of local pinch a good theorem seed has to
respect: the effect is not monotone-clean, but it has not collapsed off the
original window either.

The next Rust continuation `61 <= n <= 65` rebounds to `5/5` on the direct
comparison. At that point `n = 57` looked isolated, not like a systematic
post-55 failure. The caveat remained: the joint witness is a third witness at
`n = 63` and `n = 65`, so the statement still could not be reduced to
"mass-best always solves the whole frontier."

The later first-hit `n = 70` packet adds a second direct-comparison failure,
but the top-k frontier follow-up corrects the theorem seed. At `n = 70`, the
exact face contains witnesses that simultaneously have zero prefix residual and
zero density-adjusted mass deviation. The failure is therefore a tie-selection
artifact of first-hit reporting, not a second genuine face-level pinch. The
remaining local failure to explain in the checked window is the small `n = 57`
witness-handoff pinch.

## Candidate Finite Statement

Let `S_n` be the family of exact Sidon maximizers `A subset [0,n]` with
`|A| = h(n)`.

Define two observables on `S_n`:

- `P(A)`: prefix residual after the current general drift
  `max(||A| - sqrt(n)|, 1) * sqrt(n)`.
- `M(A)`: density-adjusted mass deviation from `n (|A| + 1) / 2`.

The finite candidate is:

> On the exact extremal surface `S_n`, low `M(A)` tends to force low `P(A)` much
> more strongly than low `P(A)` forces low `M(A)`.

The current data supports this only as a finite statement. For `10 <= n <= 50`,
the mass-best prefix penalty is no larger than the prefix-best mass penalty in
`40/41` scanned values of `n`, with the lone exception a numerical tie at
`n = 14`.

For `10 <= n <= 30`, the same comparison holds on the `h(n)` layer in `20/21`
values, but only `9/21` values on the `h(n)-1` layer.

For `51 <= n <= 55`, the beyond-50 exact-maximizer probe keeps the same direct
comparison in `5/5` values. The mass-best witness has numerically zero prefix
residual in `4/5`, the prefix-best witness has zero density-adjusted mass
deviation in `0/5`, and the joint witness is mass-best in `4/5`.

For `56 <= n <= 60`, the Rust exact-maximizer extension keeps the direct
comparison in `4/5` values and identifies the first beyond-55 exception at
`n = 57`. The joint witness is mass-best in `3/5`, prefix-best in `2/5`, and a
third witness at `n = 59`.

For `61 <= n <= 65`, the direct comparison returns to `5/5`. The mass-best
witness has numerically zero prefix residual in `3/5`, while the prefix-best
witness has zero density-adjusted mass deviation in `0/5`. The joint witness is
mass-best in `3/5` and a third witness in `2/5`.

Across the combined `51 <= n <= 65` exact-surface probes, the direct comparison
holds in `14/15` values across `49,708` exact maximizers, with `n = 57` as the
only exception in that subwindow. That is enough to keep the exact-surface theorem seed alive,
but not enough to erase the local-pinch caveat.

The `66 <= n <= 68` boundary probe keeps the direct comparison in `3/3`.
Combined over `51 <= n <= 68`, the direct comparison is now `17/18` across
`115,354` exact maximizers, with `n = 57` as the only exception through `68`. This
strengthens the exact-surface seed, but the runtime and maximizer counts now
argue for instrumented scanning before pushing much farther.

The instrumented `n = 69` packet keeps the direct comparison, bringing the
combined `51 <= n <= 69` roll-up to `18/19` across `181,766` exact maximizers.
The joint witness is a third witness at `n = 69`, so the seed should still be
phrased as one-sided compatibility on the exact surface, not as an optimizer
selection theorem.

At `n = 70`, the first-hit direct comparison fails again. The combined
first-hit `51 <= n <= 70` roll-up is `18/20` across `298,968` exact maximizers,
with failures at `n = 57` and `n = 70`. The top-k frontier packet then shows
that `n = 70` is not a face-level failure: there are joint witnesses with zero
prefix residual and zero density-adjusted mass deviation. That keeps the
exact-surface seed alive, but the correct target is now a face-aware
compatibility law rather than a single-witness optimizer-selection law.

The follow-up diagnostic `PINCH_ANALYSIS_57_70_2026-04-29.md` separates the two
failures. The `n = 57` failure is a near-tie witness handoff across `56, 57,
58`; the `n = 70` failure disappears at face level once top-k witnesses are
examined, so it should be treated as first-hit tie selection rather than
mass-center overfit.

## Weak Theorem Shape

Do not try to formalize optimizer selection first. The weaker theorem shape is:

> If `A` is an exact or certified near-extremal Sidon set and `M(A)` is close to
> minimal among exact maximizers, then `P(A)` is bounded by a smaller error term
> than the generic super-floor envelope.

That statement has three advantages:

- It avoids proving that the exact optimizer is unique.
- It avoids comparing all dense Sidon sets.
- It matches the finite stress test, which says the effect weakens one layer
  below `h(n)`.

## Why The Generic Version Is Too Strong

The generic dense-Sidon version would say that the one-sided compatibility law
holds across broad size layers. The near-maximizer pilot already argues against
that naive form.

On `|A| = h(n)-1` for `10 <= n <= 30`, the mass-best prefix penalty is no
larger than the prefix-best mass penalty in only `9/21` values. The mass-best
witness has numerically zero prefix residual in only `2/21` values. The best
joint witness is usually a third set.

That does not refute all stability. It says the stability, if real, is tied to
being on or extremely close to the extremal surface.

## Lean Boundary

The current compiled Lean theorem
`sidon_in_range_superfloor_prefix_mass_joint_envelope_external` is only a joint
envelope. It is the right formal shadow for "one object, two observables," but
it does not see `h(n)`, optimizer selection, Pareto frontiers, or finite
maximizer comparison.

So the next Lean work should not start by trying to prove the finite Pareto
statement. It should first ask for a formal handle on the extremal-surface
hypothesis: either a definition of exact maximizer families, or a weaker
near-extremal predicate that can carry both observables without pretending to
classify their optimizers.

That handle now exists in the scratch file
`Erdos30_IntervalOccupancyTarget.lean`:

- `IsMaximalSidonInRange A n` says `A` is cardinality-maximal among Sidon sets
  in `[0,n]`.
- `NearExtremalSidonInRange A n δ` says `A` is within `δ` elements of that
  maximal cardinality.
- `prefixResidualAfterGeneralDrift A n t` names the prefix residual used by the
  maximizer packets.
- `densityAdjustedMassDeviation A n` names the mass observable centered at
  `n (|A| + 1) / 2`.

The file also proves the basic projection lemmas back to `SidonInRange` and the
zero-slack inclusion `IsMaximalSidonInRange.nearExtremal_zero`. These are
definitions and bookkeeping only; they do not prove the Pareto statement.

The first wrapper theorem also now compiles:
`maximal_sidon_in_range_superfloor_prefix_mass_joint_envelope_external`. It
reuses the existing joint prefix/mass envelope under the
`IsMaximalSidonInRange` hypothesis. This is still not a compatibility theorem;
it only shows that the current envelope package can be stated on the exact
extremal surface.

## Current Status

This is a theorem seed, not a theorem. It is supported by exact enumeration and
stress-tested against the first near-maximizer layer. Its value is mainly
negative discipline: it rules out the easy generic-dense version and points the
next formal attempt at the extremal surface.
