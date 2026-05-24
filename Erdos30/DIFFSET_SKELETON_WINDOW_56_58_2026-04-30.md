# Difference-Set Skeleton Scout: n = 56..58

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_SKELETON_DIAGNOSTIC / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Source exact packet:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Derived packet:

```text
EXP-MM-030-PMF-DIFFSET-SKELETON-WINDOW-56-58-2026-04-30
```

Artifact path:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/
```

SHA256 check for the derived packet returned `OK`.

## Question

Does the `n = 57` face handoff correspond to a split in the positive-difference
skeleton, or only to a field/embedding switch inside the same skeleton?

## Answer

It is not a difference-set skeleton split at `n = 57`.

All six exported exact maximizers at `n = 57` share the same positive-difference
set. Their site sets can be far apart, and the exposed mass and joint winners
have site symmetric-difference distance `18`, but their difference-set distance
is `0`.

## Skeleton Table

| n | exact maximizers | unique difference skeletons | field split | max site distance | max difference-set distance |
|---:|---:|---:|---|---:|---:|
| 56 | 4 | 1 | false | 18 | 0 |
| 57 | 6 | 1 | true | 18 | 0 |
| 58 | 10 | 2 | true | 20 | 16 |

The first difference-skeleton branch in this local window appears at `n = 58`,
after the `n = 57` field handoff has already occurred.

## Meaning

This narrows the mechanism.

The prior tempting story was:

```text
n = 57 handoff = new difference-memory skeleton takes over
```

The packet-backed reading is sharper:

```text
n = 57 handoff = field/embedding selection inside one difference skeleton
n = 58        = first local appearance of a second difference skeleton
```

So Maxwell's Litmus still holds, but the missing information is not a new
difference skeleton at `57`. The missing information is the full exposed-face
geometry of different embeddings of the same skeleton.

## Claim Boundary

Safe:

- `n = 57 is an embedding/field-selection handoff inside one skeleton`
- `difference skeleton split appears at n = 58 in this window`
- `the skeleton scout falsifies an over-broad structural reading`

Unsafe:

- `n = 57 proves a phase transition`
- `the handoff is a difference-set bifurcation`
- `PMF proves Sidon`
- `SOTA theorem result`

## Next Gate

The next #30 gate should test whether this one-skeleton handoff is a known
Golomb-ruler/Sidon-family symmetry phenomenon:

```text
classify the shared n=57 difference skeleton up to affine/reflection
embeddings and compare its endpoint placements across indices 3 and 5.
```

That gate has now been run:

```text
AFFINE_EMBEDDING_HANDOFF_N57_2026-04-30.md
EXP-MM-030-PMF-AFFINE-EMBEDDING-N57-2026-04-30
```

The two exposed candidates are different translated embeddings of the same
skeleton:

```text
index 5 = index 3 + 1
```

The theorem-facing finite statement should therefore be about embedding
selection under fields, not skeleton replacement.
