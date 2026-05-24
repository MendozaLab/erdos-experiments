# Pruned Transfer Near-Ground D2 Probe, n = 70

Date: 2026-04-30
Problem: Erdos #30
Lane: PMF transfer operator / near-ground lattice-gas spectrum
Status: D2_NEARGROUND_PASS for `n = 70`

## Question

Does the `n = 70` first-hit artifact reading stay localized to the
zero-temperature exact maximizer face when the transfer operator retains the
second near-ground layer?

## Answer

Yes for the tested packet.

With `--prune-deficiency 2`, the layer-pruned DFS retains:

```text
|A| = h(n)
|A| = h(n)-1
|A| = h(n)-2
```

The exact top-1 prefix, mass, and joint ground-face winners remain identical to
the exact Rust top-k reference at `n = 70`. The second near-ground layer is
large, but it does not create a new visible spectral/degeneracy defect in this
packet. The sharper reading remains:

> the `n = 70` event is a first-hit / field-selection artifact on the exact
> maximizer face, not a standalone near-ground spectral singularity.

## Packet

- `EXP-MM-030-PMF-TRANSFER-PRUNED-D2-70-2026-04-30`

SHA256 sidecar verified `OK`.

## Layer Summary

| n | h(n) | ground states | h-1 states | h-2 states | terminal retained | pruned states | runtime sec | joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 70 | 10 | 117,202 | 22,801,688 | 150,565,250 | 173,484,140 | 189,156,721 | 95.237988 | ~0 |

The joint score is `1.2947641983328499e-16`, numerically zero at the current
precision.

## Frontier Parity

The D2 transfer packet matches the exact Rust top-k reference for all three
top-1 field selections:

| field | witness |
|---|---|
| prefix | `[0, 4, 12, 21, 32, 39, 55, 65, 68, 70]` |
| mass | `[0, 1, 5, 28, 44, 50, 58, 61, 68, 70]` |
| joint | `[0, 4, 12, 31, 44, 46, 49, 60, 69, 70]` |

## Interpretation

The D2 result strengthens the current field-response interpretation.

The near-ground bath is enormous:

```text
h:     117,202
h-1:   22,801,688
h-2:   150,565,250
```

But the top-1 field-selected ground-face witnesses do not move, and the exact
zero-joint witness remains present. That means the `n = 70` correction survives
one deeper near-ground layer: the first-hit failure is still better understood
as optimizer selection on a degenerate exact face than as a raw spectral defect.

The honest theorem-language candidate remains:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe: any claim that PMF supplies a Sidon proof.

## Next Gate

Do not run a broad `d = 2` sweep yet. The next useful gate is cross-problem:

```text
port the same field-response framing to #755 / B_2[3]
```

The preferred test is the verified reset/plateau window through `n = 47`. The
secondary control remains #166 sum-free, where the finite face was rigid in the
first cross-problem pass.
