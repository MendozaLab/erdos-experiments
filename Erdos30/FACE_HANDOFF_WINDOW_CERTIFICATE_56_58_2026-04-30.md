# Face-Handoff Window Certificate: n = 56..58

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_WINDOW_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Source exact packet:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Derived certificate packet:

```text
EXP-MM-030-PMF-FACE-HANDOFF-WINDOW-CERTIFICATE-56-58-2026-04-30
```

Artifact path:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/
```

SHA256 checks:

- source Maxwell V2 packet: `OK`
- derived window certificate packet: `OK`

## Question

Is the `n = 57` Maxwell Face-Handoff just a middle-row artifact, or does the
whole `56 -> 57 -> 58` exported face window support the handoff reading?

## Answer

The whole window supports the finite handoff reading.

The verified pattern is:

```text
n = 56: one exposed face
n = 57: split appears
n = 58: split persists
```

The important transition is not a raw count spike. It is a change in the
field-selected geometry of the exact maximizer face.

## Window Table

| n | h(n) | exact maximizers | exported | field split | prefix winners | mass winners | joint winners | Pareto minima |
|---:|---:|---:|---:|---|---|---|---|---|
| 56 | 10 | 4 | 4 | false | `[3]` | `[3]` | `[3]` | `[3]` |
| 57 | 10 | 6 | 6 | true | `[4,5]` | `[3]` | `[5]` | `[3,5]` |
| 58 | 10 | 10 | 10 | true | `[2,8,9]` | `[7]` | `[7]` | `[7,9]` |

All exported witnesses in all three rows pass a direct Sidon check. The
recomputed prefix, mass, joint, and Pareto roles match the source packet.

## Transition Structure

From `56` to `57`:

- all `4` maximizers from `n = 56` persist into `n = 57`;
- `2` new endpoint-shift maximizers enter at `n = 57`;
- mass selects an inherited remnant;
- joint selects an endpoint-shift state.

From `57` to `58`:

- all `6` maximizers from `n = 57` persist into `n = 58`;
- `4` new endpoint-shift maximizers enter at `n = 58`;
- the field split persists.

At `n = 57`, the exposed mass and joint candidates have symmetric-difference
distance `18`. This is a face-to-face switch, not a small local perturbation.

## Meaning

The finite object is now sharper:

> a one-face ground-state selection at `56` becomes a remnant-vs-endpoint
> exposed-face split at `57`, and that split persists at `58`.

This strengthens Maxwell's Litmus as an internal research heuristic: first-hit
accounting was missing the full face geometry. It does not upgrade the public
claim ceiling.

## Claim Boundary

Safe:

- `finite exact-face handoff`
- `remnant-to-endpoint field switch`
- `window-certified face split`
- `first-hit accounting loses face geometry`

Unsafe:

- `PMF proves Sidon`
- `physics solves Erdos #30`
- `Mendoza Limit proves the handoff`
- `phase transition proved`
- `SOTA theorem result`

## Next Gate

The structure gate has now been run:

```text
DIFFSET_SKELETON_WINDOW_56_58_2026-04-30.md
EXP-MM-030-PMF-DIFFSET-SKELETON-WINDOW-56-58-2026-04-30
```

It falsifies the over-broad skeleton-split reading at `n = 57`: all six
exported exact maximizers at `57` share the same positive-difference skeleton.
The handoff is a field/embedding switch inside one skeleton. A second local
difference skeleton first appears at `n = 58`.

The next gate is affine/reflection embedding classification of the shared
`n = 57` skeleton.
