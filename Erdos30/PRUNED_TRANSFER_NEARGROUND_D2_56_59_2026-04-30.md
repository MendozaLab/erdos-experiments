# Pruned Transfer Near-Ground D2 Probe, 56-59

Date: 2026-04-30
Problem: Erdos #30
Lane: PMF transfer operator / h=10 birth-window handoff
Status: D2_NEARGROUND_PASS for `56 <= n <= 59`

## Question

Does the `h = 10` birth window create a genuine field-handoff onset at
`n = 57`?

## Answer

Yes for the tested packet, at the finite-evidence level.

The D2 transfer operator retained:

```text
|A| = h(n)
|A| = h(n)-1
|A| = h(n)-2
```

for `n = 56,57,58,59`. The parity gate passed for all four rows: `h(n)` and
maximizer counts match the exact reference. The near-ground bath is already
large at `n = 56`, before the face-level split appears. The actual qualitative
change is on the exact maximizer face:

```text
56: prefix = mass = joint
57: prefix, mass, and joint split apart
```

That makes `57` a stronger field-handoff onset candidate than `70`.

## Packet

- `EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30`

SHA256 sidecar verified `OK`.

## Layer Summary

| n | h(n) | ground states | h-1 states | h-2 states | terminal retained | joint score |
|---|---:|---:|---:|---:|---:|---:|
| 56 | 10 | 4 | 69,564 | 4,813,066 | 4,882,634 | 0.011840 |
| 57 | 10 | 6 | 120,704 | 6,511,012 | 6,631,722 | 0.028889 |
| 58 | 10 | 10 | 200,946 | 8,684,372 | 8,885,328 | 0.036164 |
| 59 | 10 | 18 | 330,056 | 11,500,722 | 11,830,796 | 0.016531 |

## Entropy-Release Ratios

| n | `(h-1)/h` | `(h-2)/(h-1)` |
|---|---:|---:|
| 56 | 17,391.000000 | 69.189035 |
| 57 | 20,117.333333 | 53.941974 |
| 58 | 20,094.600000 | 43.217442 |
| 59 | 18,336.444444 | 34.844760 |

The near-ground bath is enormous throughout the birth window. That means the
onset signal is not simply "near-ground states become numerous." They are
already numerous at `56`.

## Frontier Witnesses

| n | prefix-field witness | mass-field witness | joint-field witness | reading |
|---|---|---|---|---|
| 56 | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | same | same | rigid top-1 field response |
| 57 | `[2, 3, 8, 12, 25, 28, 36, 43, 55, 57]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | field-handoff onset |
| 58 | `[0, 2, 15, 21, 22, 32, 46, 50, 55, 58]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | same as mass | split persists |
| 59 | `[0, 3, 7, 19, 36, 37, 46, 51, 57, 59]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | `[0, 5, 14, 16, 29, 37, 49, 55, 56, 59]` | three-way split returns |

## Interpretation

The D2 window gives a cleaner story than the earlier "57 and 70 are both
punctures" framing.

```text
56: huge near-ground bath, but rigid field winner
57: same h=10 layer, tiny ground face, field winners split
58: mass and joint align, prefix differs
59: three-way split returns
```

So the likely finite phenomenon is a field-selection handoff at the birth of
the `h = 10` face. It is not a raw spectral singularity, and it is not the same
phenomenon as `70`, where a zero-joint witness exists on a much larger exact
face.

The honest theorem-language candidate remains:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe: PMF proves Sidon, physics solves Erdos #30, or this constitutes
SOTA theorem progress.

## Next Gate

The next #30 gate should be symbolic/combinatorial, not just another count:

```text
compare the six exact h=10 maximizers at n=57 against the four exact h=10
maximizers at n=56 and the ten exact h=10 maximizers at n=58
```

The question is whether the `57` handoff can be described as a small finite
face flip: one field keeps the `56`-type structure, while another field selects
a newly available endpoint/prefix structure.
