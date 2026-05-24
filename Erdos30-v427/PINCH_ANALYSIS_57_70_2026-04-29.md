# Pinch Analysis for Erdos #30 Exact-Surface Compatibility

**Date:** 2026-04-29
**Scope:** Diagnostic comparison of existing exact-maximizer packets, not a new theorem
**Base packets:** `EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`, `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24`, `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24`, `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24`, `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-69-2026-04-24`, `EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-70-2026-04-24`
**Top-k packets:** `EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29`, `EXP-MM-030-RUST-TOPK-FRONTIER-70-2026-04-29`
**Data integrity:** Real exact enumeration from saved JSON packets

## Status

The top-k frontier scan corrected the interpretation of the pinches.

Across `51 <= n <= 70`, the direct normalized comparison holds in `18/20`
values across `298,968` exact maximizers. The two first-hit failures are
`n = 57` and `n = 70`.

But `n = 70` is not a face-level incompatibility. The top-k frontier packet
finds joint witnesses with prefix residual `0` and density-adjusted mass
deviation `0`. So the old `n = 70` failure is a first-hit tie-selection
artifact: the scanner's single stored mass-best witness over-specialized, but
the exact maximizer face itself contains perfectly compatible witnesses.

The remaining genuine pinch is `n = 57`, and even that is a tiny handoff
failure. The theorem target should therefore be face-aware compatibility, not
single-witness optimizer selection.

## Failure Table

| n | h(n) | maximizers | density slope n/h | first-hit delta | top-k face reading | joint witness |
|---|---|---:|---:|---:|---:|---:|---|
| 57 | 10 | 6 | 5.7 | +0.000192 | true small handoff pinch | prefix-best |
| 70 | 10 | 117202 | 7.0 | +0.020902 | compatible joint face exists | third witness in first-hit packet |

Here `delta = mass_best_prefix_ratio - prefix_best_mass_ratio`. Positive delta
means the first-hit direct comparison fails.

## What The Two Pinches Mean

The `n = 57` pinch is almost a tie. It sits in a tiny six-maximizer transition
zone:

| n | h(n) | maximizers | delta | joint witness |
|---|---|---:|---:|---|
| 56 | 10 | 4 | -0.011840 | mass-best and prefix-best |
| 57 | 10 | 6 | +0.000192 | prefix-best |
| 58 | 10 | 10 | -0.016488 | mass-best |

The witness handoff is visible. The mass-best witness at `n = 57` is the same
set that was the shared best witness at `n = 56`, while the prefix-best witness
at `n = 57` becomes the mass-best witness at `n = 58`. So `n = 57` looks like a
boundary crossing between adjacent exact surfaces, not a deep collapse of the
ridge.

The old `n = 70` pinch is different. It is not a near-tie in the first-hit
diagnostic, and it occurs in a much larger maximizer family:

| n | h(n) | maximizers | delta | joint witness |
|---|---|---:|---:|---|
| 68 | 10 | 36234 | -0.172256 | mass-best |
| 69 | 10 | 66412 | -0.157558 | third witness |
| 70 | 10 | 117202 | +0.020902 | third witness |

The top-k frontier scan changes the interpretation. At `n = 70`, the exact
frontier contains witnesses with prefix residual `0`, density-adjusted mass
deviation `0`, and joint score `0`. The first-hit packet stored a mass-best
witness that paid prefix cost, but that was not representative of the whole
mass-best face.

So `n = 70` is not evidence that perfect mass centering can lose prefix
compatibility. It is evidence that storing one arbitrary first-hit optimizer is
too brittle when an observable has a large tied face.

## Candidate Structural Variable

The emerging variable is not just the value of `n`, and not just whether
`n/h(n)` is an integer. The integer-slope case `n = 60` succeeds, and the
face-level `n = 70` frontier also succeeds once ties are handled correctly.

The better candidate is a frontier-selection variable:

> The finite signal lives on observable faces, not arbitrary first-hit
> optimizers. A true pinch occurs only when the relevant faces fail to intersect
> near the joint frontier.

This separates the two apparent failures:

- At `n = 57`, the failure is a handoff pinch: adjacent witnesses exchange
  roles across `56, 57, 58`, and the delta is tiny.
- At `n = 70`, the first-hit failure disappears at face level: the top-k joint
  frontier contains exact zero/zero witnesses.

## Consequence For Theorem Shape

The old seed was:

> On the exact extremal surface, mass-best usually controls prefix better than
> prefix-best controls mass.

The sharper seed is:

> On the exact extremal surface, the mass and prefix observable faces usually
> intersect or lie close on the joint frontier, except at local handoff pinches
> such as `n = 57`.

This is still finite evidence. It is not a SOTA improvement and not an
asymptotic theorem.

## Recommended Next Move

Do not jump straight to a wide `71 <= n <= 75` sweep. The better next action is
now face-aware:

1. Re-run a small top-k window around the original failure, e.g. `56 <= n <=
   58`, with `--frontier-k 10`, to confirm that `n = 57` is the only true
   face-level handoff pinch in that local band.
2. Then run `n = 71` with `--frontier-k 10 --progress`, because `n = 70` no
   longer looks like the start of a failure band once the tied face is visible.

The top-k diagnostic should stay enabled for all future exact-surface packets.
Without it, first-hit optimizer selection can manufacture false pinches.
