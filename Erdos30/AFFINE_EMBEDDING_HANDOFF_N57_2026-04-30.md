# Affine Embedding Handoff at n = 57

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_AFFINE_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Source exact packet:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Derived packet:

```text
EXP-MM-030-PMF-AFFINE-EMBEDDING-N57-2026-04-30
```

SHA256 check for the derived packet returned `OK`.

## Question

If all six `n = 57` maximizers share one difference skeleton, what exactly is
the field handoff?

## Answer

The exposed handoff is a `+1` translation inside one shared difference
skeleton.

Mass winner:

```text
index 3
[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
```

Joint winner:

```text
index 5
[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
```

Certificate relation:

```text
index 5 = index 3 + 1
```

So the handoff is not:

```text
new difference skeleton replaces old difference skeleton
```

It is:

```text
field selection shifts the embedding of the same skeleton by one site
```

## Translation Chains

The six exact maximizers at `n = 57` form two local translation chains:

```text
Chain A: 0 -> 2 -> 4 by successive +1 translations
Chain B: 1 -> 3 -> 5 by successive +1 translations
```

The exposed handoff happens on Chain B:

```text
mass field  selects index 3
joint field selects index 5 = index 3 + 1
```

The certificate also verifies reflection pairing about `n = 57`:

```text
0 <-> 5
1 <-> 4
2 <-> 3
```

## Meaning

This is the cleanest current #30 finite mechanism.

The full story is now:

```text
first-hit scan:
  sees n=57 as a failure

face export:
  shows n=57 has two exposed Pareto candidates

difference-skeleton scout:
  shows both candidates share one difference skeleton

affine certificate:
  shows the handoff is a +1 translated embedding of that skeleton
```

In Maxwell's Litmus language, the missing information was not merely "more
states." It was the normal-fan response of translated embeddings of the same
finite Sidon skeleton.

## Claim Boundary

Safe:

- `finite affine-embedding handoff`
- `translated embedding selected by field response`
- `same difference skeleton, different exposed embedding`
- `first-hit accounting misses translated face geometry`

Unsafe:

- `PMF proves Sidon`
- `physics solves Erdos #30`
- `phase transition proved`
- `SOTA theorem result`

## Next Gate

Extend the affine-embedding check across the nearest face-split rows:

```text
test whether n = 58 Pareto candidates index 7 and 9 are also linked by
translation/reflection or whether n = 58 is the first true skeleton-branch row.
```

This gate has now been run:

```text
N58_BRANCH_CERTIFICATE_2026-04-30.md
EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30
```

Result: `n = 58` mixes translated embeddings and a second difference skeleton,
but the second skeleton is prefix-side only. The Pareto face remains on the
original translated chain:

```text
n=58 index 7 = n=57 index 5
n=58 index 9 = n=58 index 7 + 1
```

The `56 -> 57 -> 58` window now has a three-stage structure:

```text
56: single exposed embedding
57: translated-embedding handoff inside one skeleton
58: first local skeleton branch appears, but Pareto remains on translated chain
```
