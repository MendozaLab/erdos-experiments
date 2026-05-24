# Observable-Split Pareto Note for Erdos #30

**Date:** 2026-04-24
**Discovery source:** `EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23`
**Machine-readable summary source:** `EXP-MM-030-COMPATIBILITY-SUMMARY-MAXIMIZER-DIAGNOSTICS-2026-04-24`
**Beyond-50 probe:** `EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`
**Rust probe:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24`
**Rust extension:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24`
**Rust boundary probe:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24`
**Rust single-n probe:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-69-2026-04-24`
**Rust second-pinch probe:** `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24`
**Scope:** Exact maximizers `A subset [0,n]` with `10 <= n <= 50`, plus bounded exact-surface probes for `51 <= n <= 70`
**Data integrity:** Real exact enumeration, no sampling

## Claim

The April 23 exact maximizer packet, now rerun with an April 24 compatibility
summary block, supports a finite, one-sided Pareto reading of the `#30`
rigidity branch. Prefix and density-adjusted mass are not behaving like two
measurements of one universal affine center. They behave like two observables
on the same Sidon object, with a strong asymmetry in which observable controls
the other.

This is finite evidence, not an asymptotic theorem.

## What The Packet Shows

The packet scanned all exact maximizers in the window `10 <= n <= 50`. The
prefix-best maximizer and the density-adjusted-mass-best maximizer coincide in
only `8` of the `41` scanned values of `n`:

`n = 11, 14, 17, 18, 25, 34, 44, 47`.

So the first conclusion is already negative: one witness usually does not
optimize both observables at once.

The stronger conclusion is asymmetric. With numerical zero interpreted as
`mass_best_prefix_residual <= 1e-9`, the mass-best witness sits on the
zero-prefix-residual face in `31/41` cases. In the opposite direction, the
prefix-best witness has exactly zero density-adjusted mass deviation in only
`1/41` case.

That is the one-sided Pareto signature: optimizing density-adjusted mass usually
preserves the best prefix face, but optimizing prefix usually does not preserve
the best density-adjusted mass face.

## Representative Witnesses

At `n = 24`, the split is strongest in normalized units. The prefix-best witness
is

`[6, 10, 16, 21, 23, 24]`,

while the density-adjusted-mass-best witness is

`[6, 8, 12, 13, 21, 24]`.

The prefix-best witness pays density-adjusted mass ratio
`0.20245547473002568` in `n^(11/8)` units. The mass-best witness pays prefix
ratio `0.07451341280198084` in `n^(7/8)` units. A third witness,

`[5, 6, 12, 16, 21, 24]`,

is the best joint compromise in the packet's summed normalized score.

At `n = 30`, the asymmetry is cleaner. The mass-best witness

`[5, 8, 9, 17, 23, 28, 30]`

has zero prefix residual and is also the best joint compromise. The prefix-best
witness

`[5, 9, 14, 20, 27, 28, 30]`

pays density-adjusted mass deviation `13`, or about `0.12103239495416415` in
`n^(11/8)` units.

At `n = 43`, the split persists later in the scanned window. The prefix-best
witness pays density-adjusted mass ratio `0.19011906974296766`, while the
mass-best witness pays only prefix ratio `0.012194693907375178`. The best joint
witness is again a third set.

## Regime Reading

The split is persistent but not monotone. It appears in `16/21` cases for
`10 <= n <= 30` and `17/20` cases for `31 <= n <= 50`. That makes the finite
signal stable across the scanned window, but it does not justify saying the
split is cleanly strengthening with `n`.

The right current sentence is:

> The exact finite data shows a persistent mixed split with strong one-sided
> Pareto flavor.

The wrong sentence is:

> The exact finite data already proves a large-`n` two-regime law.

## KvN / Holevo Inspiration Boundary

The useful physics analogy is a readout limitation, not a theorem import. In a
KvN-like picture, the Sidon set is treated as the underlying state and prefix
and density-adjusted mass are treated as observables on that state. In a
Holevo-like picture, no single readout should be assumed to expose all the
structure carried by the state.

That is exactly the role this analogy can play here: it motivates looking for a
compatibility law between observables. It does not justify saying that Sidon
sets obey a quantum information bound, and it does not replace the finite
enumeration or the Lean envelopes.

The next useful theorem seed, if one exists, should therefore look like a weak
compatibility statement: a bound describing what prefix control can force about
density-adjusted mass, or what density-adjusted mass control can force about
prefix. The April 24 compatibility-summary packet makes that direction
machine-readable because the implication looks one-sided in finite data.

## Finite Compatibility Candidate

The first candidate compatibility law is one-sided and finite:

> Among exact maximizers in this window, optimizing density-adjusted mass almost
> always leaves the prefix observable near its best face, while optimizing
> prefix often leaves visible density-adjusted mass cost.

In normalized units, the mass-best witness has prefix cost at most
`0.08838834764831843 · n^(7/8)` across the whole `10 <= n <= 50` window, with
mean normalized prefix cost about `0.010237449639906073`. The prefix-best
witness has density-adjusted mass cost as high as `0.20245547473002568 ·
n^(11/8)`, with mean normalized mass cost about `0.08800402748216497`.

Comparing the two natural normalized penalties directly, the mass-best prefix
penalty is no larger than the prefix-best density-adjusted mass penalty in
`40/41` scanned values of `n`; the only exception is the numerical tie at
`n = 14`, where the mass-best prefix penalty is about `8.8e-17` and the
prefix-best mass penalty is exactly zero.

The joint-score diagnostic points in the same direction. The best joint witness
equals the mass-best witness in `31/41` values of `n`, and it also equals the
prefix-best witness in `8/41`; those counts overlap when one witness optimizes
both observables. A third witness is joint-best in `10/41` values. So the finite
frontier is not simply "mass beats prefix." It is more precise: the
mass-optimal face is often already compatible with joint optimization, while
the prefix-optimal face is rarely enough by itself.

The bounded beyond-50 probe did not break that reading. In
`EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`, exact
maximizers for `51 <= n <= 55` still satisfy the direct normalized comparison
in `5/5` values. The mass-best witness has numerically zero prefix residual in
`4/5` values, the prefix-best witness has zero density-adjusted mass deviation
in `0/5`, and the joint witness equals the mass-best witness in `4/5`. The
probe is too small to turn this into an asymptotic claim, but it matters because
it shows the exact-surface signal survives past the original `10 <= n <= 50`
packet.

The Rust port then extended the exact-surface stress test to `56 <= n <= 60`.
That packet preserves the broad direction but adds useful friction: the direct
normalized comparison holds in `4/5` values, with a real exception at `n = 57`.
The joint witness equals the mass-best witness in `3/5`, equals the prefix-best
witness in `2/5`, and is a third witness at `n = 59`. So the right update is
not "the ridge is clean." It is "the ridge persists, but with local pinches that
the theorem seed must be able to explain."

The port itself also matters operationally. A Rust replay of `51 <= n <= 55`
matched the prior Python packet on the checked maximizer counts, witness
choices, and compatibility diagnostics, while reducing total runtime from about
`193.3s` to about `9.93s`. The `56 <= n <= 60` Rust packet completed in about
`23.8s`, making larger exact-surface scouting practical without changing the
mathematical contract.

The next Rust packet, `61 <= n <= 65`, made the `n = 57` exception look
isolated at first. The direct normalized comparison holds in `5/5` values, the
mass-best witness has numerically zero prefix residual in `3/5`, and the
prefix-best witness again has zero density-adjusted mass deviation in `0/5`.
The joint witness is mass-best in `3/5` and a third witness in `2/5`, so the
frontier is still not collapsing to one universal optimizer. The strongest
split in this packet is at `n = 65`.

Aggregating the bounded exact-surface probes from `51 <= n <= 65`, the direct
comparison holds in `14/15` values across `49,708` exact maximizers, with
`n = 57` the only exception in that subwindow. The mass-best witness has numerically zero
prefix residual in `9/15`, while the prefix-best witness has zero
density-adjusted mass deviation in `0/15`. That is the current best finite
sentence: the one-sided compatibility ridge is real enough to keep probing, but
not clean enough to promote into a monotone law.

The bounded `66 <= n <= 68` Rust probe keeps the ridge intact: the direct
comparison holds in `3/3`, with the joint witness mass-best at `n = 67, 68` and
a third witness at `n = 66`. The roll-up from `51 <= n <= 68` is therefore
`17/18` direct compatibility across `115,354` exact maximizers, with `n = 57`
the only failure through `68`. The engineering caveat is now visible: the `66 <= n <=
68` packet took about `76.9s`, and `n = 68` alone has `36,234` maximizers, so
the next wider scan should add progress or pruning rather than treating exact
enumeration as free.

The instrumented `n = 69` single probe keeps the direct comparison again. It
has `66,412` exact maximizers, `1.52B` recursive nodes, and took about `38.8s`
with progress telemetry enabled. The mass-best prefix penalty remains below
the prefix-best density-adjusted mass penalty, but the joint witness is a third
witness. The roll-up from `51 <= n <= 69` is now `18/19` direct compatibility
across `181,766` exact maximizers, with `n = 57` still the only failure through
`69`.
That strengthens the ridge while preserving the frontier interpretation.

The next single point, `n = 70`, is a second first-hit direct-comparison
failure. It has `117,202` exact maximizers, `1.82B` recursive nodes, and took
about `48.4s`. Here the first mass-best witness pays more prefix penalty than
the first prefix-best witness pays density-adjusted mass penalty, so the
first-hit direct comparison fails. The current first-hit roll-up from `51 <= n
<= 70` is therefore `18/20` direct compatibility across `298,968` exact
maximizers, with failures at `n = 57` and `n = 70`.

The top-k frontier follow-up changes the meaning of that failure. At `n = 70`,
the exact face contains joint witnesses with zero prefix residual and zero
density-adjusted mass deviation. So `n = 70` is not a second genuine pinch; it
is a tie-selection artifact exposed by using only one first-hit witness per
observable. The remaining face-level pinch in this checked window is `n = 57`.

## Near-Maximizer Robustness Check

The Tao-style stress test is to ask whether this is an exact-optimizer artifact.
The pilot packet
`EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24` scans both `|A| =
h(n)` and `|A| = h(n)-1` for `10 <= n <= 30`.

On the `h(n)` layer, the one-sided pattern remains strong: the mass-best prefix
penalty is no larger than the prefix-best density-adjusted mass penalty in
`20/21` values, with the lone exception again at `n = 14`. The mass-best witness
has numerically zero prefix residual in `16/21` values, while the prefix-best
witness has zero density-adjusted mass deviation in only `1/21`.

On the `h(n)-1` layer, the pattern weakens sharply. The mass-best prefix penalty
is no larger in only `9/21` values, with exceptions at
`n = 10, 14, 15, 16, 20, 21, 22, 23, 24, 28, 29, 30`. The mass-best witness has
numerically zero prefix residual in only `2/21` values, and the best joint
witness is usually a third set rather than the mass-best set.

So the current finite reading should be sharpened: the one-sided compatibility
signal appears tied to the exact extremal surface. It is not yet a generic
dense-Sidon stability law one layer below the maximum.

The corresponding theorem seed is recorded separately in
`EXTREMAL_SURFACE_COMPATIBILITY_SEED_2026-04-24.md`. The main discipline is that
future formal statements should condition on the exact extremal surface, or on a
near-extremal predicate strong enough to behave like that surface. The naive
generic dense-Sidon version is already too strong for the current finite data.

## Formal Status

The current Lean theorem
`sidon_in_range_superfloor_prefix_mass_joint_envelope_external` is the honest
formal shadow of this observation. It says one Sidon set carries simultaneous
prefix and density-adjusted mass envelopes. It does not prove a Pareto law, an
optimizer-selection theorem, or an incompatibility statement.

So the next theorem-side target should be weaker than a full tradeoff theorem:
look first for a compatibility shadow between the two observables, or keep this
as a finite Pareto note until a formal statement becomes visible.

## Consequence For The Search

The old local question was:

> Which affine center fixes both prefix and mass?

The current better question is:

> What common Sidon structure is being read out differently by prefix and
> density-adjusted mass?

That is why the next `#30` work should not simply chase another universal
recentered coordinate. The data is pointing toward observables, frontiers, and
compatibility, not one scalar correction that makes the two readouts collapse.
