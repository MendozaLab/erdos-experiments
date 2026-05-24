# Face-Aware Frontier Backfill for Erdos #30

Date: 2026-04-29
Packet: `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`
Scanner: `rust-exact-maximizer`
Mode: exact enumeration, no sampling
Window: `51 <= n <= 70`

## What Changed

The earlier exact-surface scans reported one witness per observable. That was
good enough to see the ridge, but it mixed two different things:

1. whether the exact maximizer face contains compatible witnesses;
2. which tied witness the scanner happened to report first.

The new Rust scanner keeps the old first-hit fields for packet parity, but also
stores bounded top-k frontier witnesses for:

- prefix-best witnesses;
- density-adjusted-mass-best witnesses;
- joint-best witnesses.

This gives a face-aware diagnostic instead of a single-witness diagnostic.

## Backfill Result

The full `51 <= n <= 70` backfill preserves the old first-hit statistic:

- First-hit direct compatibility: `18/20`.
- First-hit failures: `n = 57` and `n = 70`.
- Exact maximizers scanned: `298,968`.
- Total packet runtime: `284.640270375` seconds.

But the top-k frontier layer changes the interpretation:

- `n = 70` has joint witnesses with near-zero prefix ratio and zero
  density-adjusted mass ratio.
- Therefore `n = 70` is a first-hit tie-selection artifact, not a true
  face-level incompatibility.
- `n = 57` remains the important handoff pinch: the best joint witness has
  near-zero prefix ratio but nonzero density-adjusted mass ratio
  `0.02888905507479718`.

## Joint Frontier Shape

Near-zero joint witnesses appear at:

`51, 52, 53, 54, 60, 62, 64, 66, 68, 70`

Nonzero joint-best scores appear at:

`55, 56, 57, 58, 59, 61, 63, 65, 67, 69`

The largest joint score in this backfill occurs at `n = 58`:

- joint score: `0.03616352164892123`;
- prefix ratio: `0.028641813045612422`;
- mass ratio: `0.007521708603308812`;
- witness: `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]`.

That matters because it blocks the too-strong statement "every exact face has a
perfect zero-zero joint witness." The actual ridge is softer: exact maximizer
faces often contain perfect or near-perfect joint witnesses, and the top-k
frontier shows when apparent failures are merely tied-witness selection.

## Theorem Consequence

Do not aim first for:

> the mass-best optimizer is always prefix-compatible.

That is too sensitive to tied faces and scanner selection.

Aim instead for:

> on the exact Sidon maximizer face, low mass deviation and low prefix residual
> have a shared small-joint frontier, with isolated handoff pinches.

In the current checked window, `n = 57` is the real handoff pinch. `n = 70` is
not.

## Next Probe

The next probe was run as:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30/rust-exact-maximizer
./target/release/erdos30-exact-maximizer \
  --n-min 71 \
  --n-max 71 \
  --frontier-k 10 \
  --progress \
  --experiment-id EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29 \
  --output-dir /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30
```

Result:

- Exact maximizers scanned: `203,840`.
- Runtime: `56.872648209` seconds.
- Search nodes: `2,185,808,767`.
- First-hit direct compatibility: pass.
- Mass-best prefix ratio: effectively zero.
- Prefix-best mass ratio: `0.05268554766780052`.
- Best joint score: `0.0014239337207517064`.

Interpretation: `n = 71` rebounds. It is not a second `57`-style handoff
failure. The top-k joint frontier again finds witnesses with near-zero prefix
cost and very small density-adjusted mass cost. This strengthens the
face-aware ridge reading: the visible obstruction is local handoff structure,
not a broad post-70 collapse.
