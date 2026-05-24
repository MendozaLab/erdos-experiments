# B_2[g] Capped Near-Ground Scout

Date: 2026-04-29
Problem: Erdős #755 candidate lane
Status: CAPPED_LOWER_BOUND_SCOUT

## Question

Can a capped near-ground scan keep both capacity lanes moving without paying
for full D2 enumeration?

## Source Packets

- `EXP-MM-755-PMF-B2G2-TRANSFER-CAPPED-D2-41-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-CAPPED-D2-39-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-41-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-LAYER-CAPPED-D2-39-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-42-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-40-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-IMPORTED-GROUND-LAYER-CAPPED-D2-40-FIXED-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-43-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-41-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-44-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-IMPORTED-GROUND-LAYER-CAPPED-D2-41-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-45-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-42-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-46-2026-04-30`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-43-2026-04-30`
- `EXP-MM-755-PMF-B2G2-TRANSFER-GROUNDONLY-47-2026-04-30`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-44-2026-04-30`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-45-2026-04-30`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-46-2026-04-30`

Both packets have verified SHA256 sidecars.

## Result

The global cap works as an engineering guard, but it is not a good ratio
estimator.

Both rows hit the `10,000,000` retained-terminal cap. That confirms the D2 bath
is already large in the next rows while preserving exact ground-face
enumeration.

| g | n | h | exact ground | h-1 lower | h-2 lower | retained cap hit? | near-ground counts exact? |
|---:|---:|---:|---:|---:|---:|---|---|
| 2 | 41 | 13 | 31,290 | 175,439 | 9,824,546 | yes | no |
| 3 | 39 | 16 | 156,522 | 38,420 | 9,961,580 | yes | no |

The important caution: because the cap is global, the traversal order can fill
the cap mostly from the `h-2` bath before `h-1` is adequately sampled. Therefore
`h-1 / ground` from capped packets is a lower bound, not a compression-ratio
measurement.

## Interpretation

The capped scout is useful for runtime control and lower-bound evidence:

```text
ground face remains exact
D2 bath is already at least cap-sized
near-ground ratios are censored lower bounds
```

It is not yet the right replacement for exact D2 enumeration. The next sampler
needs a per-layer cap so that `h-1` and `h-2` each get their own lower-bound
budget.

The MDL reading is deliberately narrow. Once the near-ground bath crosses the
enumeration boundary, the object stops looking like a list of sets and starts
looking like a compressed state generator. That is why transfer operators and
tensors become plausible next tools. But the claim remains computational:
physics intuition points to the compression test; exact packets decide whether
the compression is real.

## Layer-Capped Follow-Up

The scanner now supports:

```text
--near-ground-layer-cap N
```

This stops the second pass only after each retained non-ground deficiency layer
reaches its own cap. With a layer cap of `1,000,000`, both lanes hit the cap:

| g | n | h | exact ground | h-1 lower | h-2 lower | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---|---|---:|
| 2 | 41 | 13 | 31,290 | 1,000,000 | 1,000,000 | yes | no | 684,077,732 |
| 3 | 39 | 16 | 156,522 | 1,000,000 | 1,000,000 | yes | no | 2,654,116,969 |
| 2 | 42 | 13 | 108,218 | 1,000,000 | 1,000,000 | yes | no | 684,077,733 |

This is the better capped diagnostic. It still does not give exact ratios, but
it avoids the global-cap failure mode where `h-2` can consume almost the whole
budget before `h-1` is sampled.

The `B_2[3]` ground-only row at `n=40` found the next sparse reset:

| g | n | h | exact ground | field split? | second-pass nodes |
|---:|---:|---:|---:|---|---:|
| 3 | 40 | 17 | 8 | no | 0 |

This row did not run D2; it exists to avoid repeating an expensive ground pass
just to learn whether a new handoff occurred.

The scanner now also supports importing an already verified ground result:

```text
--known-h H --known-ground-count C
```

That allowed the `B_2[3] n=40` D2 shell to run without repeating the
8-billion-node ground pass. The first imported-ground packet exposed a
diagnostic bug: if one retained deficiency layer hit its per-layer cap while
another layer remained below cap, the old code could label the packet exact.
The scanner was patched so a skipped terminal state in any capped layer marks
the packet censored. The corrected packet is the citation target:

| g | n | h | exact ground | h-1 | h-2 lower | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---|---|---:|
| 3 | 40 | 17 | 8 | 810,894 | 1,000,000 | yes | no | 14,733,757,365 |

This is the cleanest post-reset lower-bound row so far. The ground face is
sparse again, the first excited shell is already huge, and the second shell is
censored at the layer cap.

The `B_2[2]` lane then produced its next reset at `n=43`:

| g | n | h | exact ground | h-1 | h-2 lower | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---|---|---:|
| 2 | 43 | 14 | 18 | 352,990 | 1,000,000 | yes | no | 2,809,606,044 |

So the bounded-sum lanes now show paired reset rows:

```text
B_2[3]: h=17 reset at n=40, ground=8
B_2[2]: h=14 reset at n=43, ground=18
```

The `B_2[3]` row at `n=41` shows the reset beginning to expand into a plateau:

| g | n | h | exact ground | field split? | joint-best score | second-pass nodes |
|---:|---:|---:|---:|---|---:|---:|
| 3 | 41 | 17 | 246 | yes | 0.11052594744024911 | 0 |

This is not a near-zero joint-face row. It is still field-sensitive, but the
joint optimum is comparatively costly. That keeps the story honest: the
handoff/reset structure persists, while the observable frontier quality varies
substantially across the plateau.

The matching `B_2[2]` plateau row is:

| g | n | h | exact ground | h-1 | h-2 lower | joint-best score | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| 2 | 44 | 14 | 122 | 997,280 | 1,000,000 | 0.14477567886658796 | yes | no | 4,115,933,005 |

That is the lower-capacity analogue of `B_2[3] n=41`: after the reset, the
ground face expands, the first-excited shell nearly reaches the million cap,
and the field frontier remains costly rather than near-perfect.

The imported-ground D2 row for `B_2[3] n=41` confirms that the first two
near-ground shells are already beyond the million-per-layer scout budget:

| g | n | h | exact ground | h-1 lower | h-2 lower | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---|---|---:|
| 3 | 41 | 17 | 246 | 1,000,000 | 1,000,000 | yes | no | 17,374,431,898 |

So the `g=3` plateau expands in both layers immediately after the sparse
reset. We no longer have an exact `h-1` ratio at `n=41`; we have a censored
lower bound showing the first-excited shell crossed the scout ceiling.

The `B_2[2]` plateau reaches the same censored regime by `n=45`:

| g | n | h | exact ground | h-1 lower | h-2 lower | joint-best score | layer cap hit? | counts exact? | second-pass nodes |
|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| 2 | 45 | 14 | 724 | 1,000,000 | 1,000,000 | 0.10793650793650794 | yes | no | 4,247,219,068 |

This is now the operational boundary for exact D2 ratios in the `g=2` lane:
after `n=44`, the million-per-layer scout becomes censored in both layers.

The next `B_2[3]` ground-only row continues the same plateau growth:

| g | n | h | exact ground | field split? | joint-best score | second-pass nodes |
|---:|---:|---:|---:|---|---:|---:|
| 3 | 42 | 17 | 3,134 | yes | 0.1029078297985861 | 0 |

The `g=3` ground face therefore expands:

```text
n=40: 8
n=41: 246
n=42: 3,134
n=43: 32,002
n=44: 212,586
```

That is a clean reset-to-plateau expansion sequence at the exact maximizer
face. The joint-best score falls across the plateau:

```text
n=41: 0.11052594744024911
n=42: 0.1029078297985861
n=43: 0.08724906944930486
n=44: 0.07477515799708315
```

The rows are still not near-zero joint witnesses, but the field tension is
compressing as the plateau thickens.

The next row, `B_2[3] n=45`, is the next handoff:

| g | n | h | exact ground | field split? | joint-best score | second-pass nodes |
|---:|---:|---:|---:|---|---:|---:|
| 3 | 45 | 18 | 8 | yes | 0.1454320987654321 | 0 |

So the high-capacity lane now shows two full reset cycles:

```text
h=17 reset: n=40, ground=8
h=17 plateau: n=41..44, ground=246 -> 3,134 -> 32,002 -> 212,586
h=18 reset: n=45, ground=8
h=18 plateau begins: n=46, ground=142
```

This is the cleanest finite handoff evidence in the #755 lane so far.

The `n=46` row is also an important exception: the top-k frontier did not split
for the tracked fields, even though the ground face expanded. That keeps the
claim scoped to reset/plateau mechanics rather than universal field splitting.

The `B_2[2]` lane now has the matching post-reset exact-face sequence:

| g | n | h | exact ground | h-1 lower | h-2 lower | joint-best score | counts exact? |
|---:|---:|---:|---:|---:|---:|---:|---|
| 2 | 43 | 14 | 18 | 352,990 | 1,000,000 | n/a | no |
| 2 | 44 | 14 | 122 | 997,280 | 1,000,000 | 0.14477567886658796 | no |
| 2 | 45 | 14 | 724 | 1,000,000 | 1,000,000 | 0.10793650793650794 | no |
| 2 | 46 | 14 | 3,504 | 1,000,000 | 1,000,000 | 0.10889819065622469 | no |
| 2 | 47 | 14 | 16,036 | n/a | n/a | 0.09422492401215805 | no |

Ground-face expansion:

```text
n=43: 18
n=44: 122
n=45: 724
n=46: 3,504
n=47: 16,036
```

## Claim Boundary

Safe internal claim:

> The capped D2 scouts show that the next `B_2[2]` and `B_2[3]` rows already
> have large near-ground baths while preserving exact ground-face counts.

Unsafe:

> Capped lower bounds preserve the exact compression ratios.

Unsafe:

> The capped packets extend the exact D2 ratio table.

## Next Gate

The next bottleneck is now the exact ground pass, especially for `g=3`.

```text
B_2[2] layered-capped D2 at n = 43
B_2[3] ground-only at n = 47
B_2[2] ground-only at n = 48
```
