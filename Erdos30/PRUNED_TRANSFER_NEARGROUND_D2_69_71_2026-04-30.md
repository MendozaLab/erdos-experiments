# Pruned Transfer Near-Ground D2 Probe, 69-71

Date: 2026-04-30
Problem: Erdos #30
Lane: PMF transfer operator / near-ground lattice-gas spectrum
Status: D2_NEARGROUND_PASS for `69 <= n <= 71`

## Question

Does the `n = 70` field-selection reading survive when the same `d = 2`
near-ground transfer probe is run uniformly across `69 <= n <= 71`?

## Answer

Yes for the tested packet.

The D2 transfer operator retained:

```text
|A| = h(n)
|A| = h(n)-1
|A| = h(n)-2
```

for `n = 69,70,71`. The parity gate passed for all three rows: `h(n)` matched
the exact reference and maximizer counts matched the exact reference. The
second near-ground layer is large and smooth across the window. It does not
turn `n = 70` into a visible near-ground degeneracy defect.

The sharper reading is:

> `n = 70` is special in the zero-temperature field-selected joint score, not
> in the raw near-ground state bath.

This strengthens the current view that `70` was a first-hit / field-selection
artifact on a degenerate exact maximizer face.

## Packet

- `EXP-MM-030-PMF-TRANSFER-PRUNED-D2-69-71-2026-04-30`

SHA256 sidecar verified `OK`.

## Layer Summary

| n | h(n) | ground states | h-1 states | h-2 states | terminal retained | pruned states | joint score |
|---|---:|---:|---:|---:|---:|---:|---:|
| 69 | 10 | 66,412 | 16,623,922 | 122,996,090 | 139,686,424 | 165,499,468 | 0.001481 |
| 70 | 10 | 117,202 | 22,801,688 | 150,565,250 | 173,484,140 | 189,156,721 | ~0 |
| 71 | 10 | 203,840 | 31,058,406 | 183,623,966 | 214,886,212 | 215,603,601 | 0.001424 |

## Entropy-Release Ratios

| n | `(h-1)/h` | `(h-2)/(h-1)` |
|---|---:|---:|
| 69 | 250.315033 | 7.398741 |
| 70 | 194.550332 | 6.603250 |
| 71 | 152.366591 | 5.912215 |

The first-defect and second-defect ratios both move smoothly downward across
`69,70,71`. That is not the shape of a raw near-ground spectral singularity at
`70`. It is the shape of a smooth near-ground bath with a zero-temperature
field-selection event on the exact face.

## Frontier Witnesses

| n | prefix-field witness | mass-field witness | joint-field witness |
|---|---|---|---|
| 69 | `[0, 3, 11, 20, 33, 43, 62, 64, 68, 69]` | `[0, 1, 5, 28, 43, 49, 57, 60, 67, 69]` | `[0, 3, 11, 21, 33, 49, 62, 64, 68, 69]` |
| 70 | `[0, 4, 12, 21, 32, 39, 55, 65, 68, 70]` | `[0, 1, 5, 28, 44, 50, 58, 61, 68, 70]` | `[0, 4, 12, 31, 44, 46, 49, 60, 69, 70]` |
| 71 | `[0, 4, 13, 21, 40, 45, 56, 68, 70, 71]` | `[0, 1, 6, 31, 44, 53, 55, 63, 67, 70]` | `[0, 4, 13, 23, 34, 51, 63, 65, 66, 71]` |

## Interpretation

The D2 result strengthens the finite lattice-gas reading without upgrading it
to a theorem. The exact maximizer face is highly degenerate, the near-ground
bath is enormous, and the field-selected frontier is stable under the D2
retention rule.

The honest theorem-language candidate remains:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe: any claim that PMF supplies a Sidon proof, that physics solves
Erdos #30, or that this is a SOTA theorem result.

## Next Gate

The next useful #30 gate is no longer another raw D2 count in this same window.
The count signal is smooth. The next gate should target the field geometry:

```text
build a compact witness-distance / field-handoff table for 51 <= n <= 71
```

The key measurement should compare prefix, mass, and joint witnesses by
overlap, symmetric difference, endpoint/prefix location, mass deviation, and
joint score. The question is whether `57` and `70` are visibly different
handoff types or two instances of the same face-selection mechanism.
