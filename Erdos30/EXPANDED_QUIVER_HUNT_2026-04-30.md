# Expanded Quiver Hunt

Date: 2026-04-30
Problem: Erdos #30, with cross-problem flags for #755 and #166
Status: working hunt note, evidence-bound

## Why This Exists

The current quiver is larger than "does the Sidon scan pass?"

The useful arrows now are:

- field response on an exact maximizer face;
- entropy release from `h` to `h-1` and `h-2`;
- geometric observables as projectors onto a relevant Hilbert sector;
- tensor / low-rank compression tests on that sector;
- cross-problem controls using `B_2[g]` and sum-free models.

The public claim ceiling is unchanged. This is finite evidence and proof-target
generation, not a Sidon proof.

## Arrow 1: Information-Release Ratios

The verified `n = 70`, `d = 2` packet is:

```text
EXP-MM-030-PMF-TRANSFER-PRUNED-D2-70-2026-04-30
```

Layer counts:

```text
h:     117,202
h-1:   22,801,688
h-2:   150,565,250
```

Entropy-release ratios:

| step | ratio | ln ratio | bits |
|---|---:|---:|---:|
| `(h-1) / h` | 194.550332 | 5.270691 | 7.604000 |
| `(h-2) / (h-1)` | 6.603250 | 1.887562 | 2.723176 |
| `(h-2) / h` | 1284.664511 | 7.158253 | 10.327176 |

Interpretation:

> The first defect releases about 7.6 bits of configurational freedom, while
> the second releases about 2.7 bits. That is an information-theoretic jamming
> signature: the exact face is rigid, the first near-ground bath is huge, and
> the second step adds freedom much more slowly.

This is a combinatorial / information-theoretic shadow, not a Ramanujan
mock-theta shadow.

## Arrow 2: Is `70` a Spectral Defect?

The D1 first-defect entropy releases are:

| n | ground | h-1 | `(h-1)/ground` | bits |
|---:|---:|---:|---:|---:|
| 69 | 66,412 | 16,623,922 | 250.315033 | 7.967601 |
| 70 | 117,202 | 22,801,688 | 194.550332 | 7.604000 |
| 71 | 203,840 | 31,058,406 | 152.366591 | 7.251403 |

The first-defect entropy release decreases smoothly across `69, 70, 71`. That
supports the current reading:

> `70` is not a raw near-ground spectral spike. The visible event is still a
> zero-temperature field-selection artifact on the exact maximizer face.

## Arrow 3: Geometry as Hilbert-Space Projector

The naive Hilbert space is too large:

```text
H_n = span{|A> : A subset [0,n]}
```

The relevant sector is smaller:

```text
H_relevant = span{|A> :
  A is Sidon-admissible,
  |A| in {h, h-1, h-2},
  prefix residual and density-adjusted mass are small
}
```

The geometric observables are therefore not decoration. They act like
projectors:

```text
all subsets -> Sidon-admissible states -> near-ground shells -> low-joint sector
```

Safe paper language:

> The Hilbert-space view becomes useful when geometric observables act as
> projectors onto a much smaller relevant sector. Prefix balance and
> density-adjusted mass do not merely score witnesses; they identify the
> low-description subspace where field response is visible.

## Arrow 4: Tensor / Low-Rank Test

The next serious compression test is tensorial.

Candidate construction:

```text
T_layer[source_state, removed_or_added_site, target_state]
```

or a bipartite incidence operator between adjacent layers:

```text
I_{h,h-1}(A,B) = 1 if B is obtained from A by deleting one occupied site
I_{h-1,h-2}(B,C) = 1 if C is obtained from B by deleting one occupied site
```

Question:

> Do prefix/mass/joint-projected slices have lower effective rank than random
> or unprojected slices?

If yes, Hilbert/tensor language is compression evidence, not vocabulary.

## Arrow 5: #755 Reset / Plateau Control

The verified `B_2[3]` high-capacity lane through `n = 47`:

| transition | ground counts | ratio | ln ratio | bits |
|---|---:|---:|---:|---:|
| `45 -> 46` | `8 -> 142` | 17.750000 | 2.876386 | 4.149747 |
| `46 -> 47` | `142 -> 2160` | 15.211268 | 2.722036 | 3.927068 |

Interpretation:

> The #755 lane shows plateau expansion after a sparse reset. The ground-face
> expansion ratio is far smaller than the #30 first-defect entropy release, but
> it has the same finite laboratory flavor: a rigid reset followed by controlled
> degeneracy growth.

Next gate:

```text
Run or reuse a #755 field-response table over n = 45,46,47:
- split or non-split?
- joint score trend?
- does field response stabilize as the plateau expands?
```

Existing verified packets answer the first pass:

| n | h | ground | field split? | joint follows | joint score |
|---:|---:|---:|---|---|---:|
| 45 | 18 | 8 | yes | prefix | 0.145432 |
| 46 | 18 | 142 | no | prefix = mass = joint | 0.116887 |
| 47 | 18 | 2160 | yes | prefix | 0.101378 |

Interpretation:

> The #755 reset/plateau lane is not a universal split story. It is a
> field-response flag with exceptions: the sparse reset row splits, the next
> plateau row is non-split, and the larger plateau row splits again. That is
> more useful than a clean slogan because it gives the next tensor/Hilbert test
> a concrete target: explain when field selections collapse to one witness and
> when they separate.

## Current Best Hunt Targets

1. **#30 tensor incidence test at n = 70.**
   Use the `h`, `h-1`, and `h-2` layers from the D2 packet. Start with sampled
   or projected incidence if full materialization is too large.

2. **#755 field-response mechanism.**
   Explain the split / non-split / split pattern across `n = 45,46,47`, with
   joint following prefix in the split rows.

3. **#166 negative control.**
   Keep it as the rigid-face comparator. If sum-free stays unsplit while Sidon
   and `B_2[g]` split, that strengthens the field-response story.

4. **Entropy-release smoothness check.**
   If affordable, get `d = 2` for `n = 69` or `71` to see whether the second
   defect release is also smooth around `70`.

## Claim Boundary

Safe:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Safe:

> geometry appears to project the huge Hilbert space onto a smaller relevant
> sector where field response is visible.

Unsafe:

> these packets prove a Sidon theorem.

Unsafe:

> tensor language is already a proof simplifier.
