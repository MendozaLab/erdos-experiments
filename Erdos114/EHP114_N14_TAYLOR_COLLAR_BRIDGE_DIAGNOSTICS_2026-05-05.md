# EHP114 n=14 Taylor/Collar Bridge Diagnostics

Date: 2026-05-05

## Meaning

After root-factored `p'` failed to improve the bridge, two successor routes
were tested on the same hard `n=14`, `eps=0.1`, worst subcell `(6,4)`:

1. Taylor/Lipschitz lower bounds for `p'`.
2. A validated level-set collar filter that rejects boxes whose center value
   cannot reach the level set under the local gradient upper bound.

The combined route is the first one that materially moves the bridge. It
resolves regularity on the retained boxes and cuts the normal error/budget ratio
to about `2.19`. It is still over budget, so this is not a proof.

This remains a shadow signature, not universal law. It is not an exact
lemniscate-length certificate, not a Lean theorem, and not a proof of Erdős
#114.

## Verified Local Artifacts

All three runs use:

```text
engine = erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_bridge_diagnostic.rs
experiment_id = EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
degree = 14
eps = 0.1
subcell = (6,4)
z_subdivision = 16
derivative_mode = recurrent
```

Hash checks:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114
for d in bridge-diagnostic-taylor-lipschitz-output-z16 bridge-diagnostic-level-set-collar-output-z16 bridge-diagnostic-taylor-collar-output-z16; do
  (cd "$d" && shasum -a 256 -c *_RESULTS.sha256)
done
```

All three returned `OK`.

## Result Matrix

```text
ambient z16 baseline:
  unresolved cells        = 13560
  normal error / budget   = 58.40007122458721
  relative error / budget = 165.17900097389062

taylor-lipschitz:
  unresolved cells        = 0
  active boxes            = 95071
  collar rejected boxes   = 0
  min gradient lower      = 1.6950326992351197
  normal error / budget   = 9.720058113201167
  relative error / budget = 27.491330529471572

level-set-collar:
  unresolved cells        = 4443
  active boxes            = 20935
  collar rejected boxes   = 74136
  min gradient lower      = 0.0
  normal error / budget   = 11.545178669554264
  relative error / budget = 32.654380979830194

taylor-collar:
  unresolved cells        = 0
  active boxes            = 20935
  collar rejected boxes   = 74136
  min gradient lower      = 2.4845471318753294
  normal error / budget   = 2.1882768069742764
  relative error / budget = 6.189092230763413
```

## Interpretation

Taylor/Lipschitz arithmetic solves the regularity obstruction but not the
budget. The level-set collar removes most ambient boxes but cannot resolve
regularity alone. Together, they turn the blocker from "regularity unresolved"
into "constants still too loose."

That is progress. The next attack should not return to root-factored `p'` or
blind subdivision. It should tighten the Taylor-collar constants:

1. Replace the crude global `p''` upper bound in each retained box with a
   sharper local Taylor/Bernstein bound.
2. Split the retained collar boxes adaptively by condition ratio, not uniformly.
3. Recompute the normal-error budget first; the relative-length budget is still
   farther away.

The immediate target is:

```text
taylor-collar normal error / budget < 1
```

on the same `(6,4)` cell before scaling to the 164-cell packet.

