# Pruned Transfer Near-Ground D1 Probe, 69-71

Date: 2026-04-29
Problem: Erdős #30
Lane: PMF transfer operator / near-ground lattice-gas spectrum
Status: D1_NEARGROUND_PASS for `69 <= n <= 71`

## Question

Does the 69-71 field-frontier result survive when the transfer operator retains
the first near-ground layer, not just the exact maximizer face?

## Answer

Yes for frontier parity, with a narrower spectral interpretation.

With `--prune-deficiency 1`, the layer-pruned DFS retains both:

```text
|A| = h(n)
|A| = h(n)-1
```

The exact top-1 prefix, mass, and joint ground-face winners remain unchanged
from the exact Rust maximizer packets. The exact-vs-transfer top-1 diff over
`69 <= n <= 71` was empty.

The new near-ground information is that the first-excited layer is huge and
smoothly increasing. It does not show a visible singularity at `70` or `71`.
That supports the current reading: the interesting PMF signal is a
zero-temperature field response on the exact maximizer face, not a raw
cardinality-layer spike.

## Packets

- `EXP-MM-030-PMF-TRANSFER-PRUNED-D1-69-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-D1-70-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-D1-71-2026-04-29`

All three SHA256 sidecars verified `OK`.

## Layer Summary

| n | h(n) | ground states | first-excited states | terminal retained | pruned states | runtime sec | joint score |
|---|---:|---:|---:|---:|---:|---:|---:|
| 69 | 10 | 66,412 | 16,623,922 | 16,690,334 | 244,310,914 | 35.101241 | 0.001481 |
| 70 | 10 | 117,202 | 22,801,688 | 22,918,890 | 288,495,558 | 38.715517 | ~0 |
| 71 | 10 | 203,840 | 31,058,406 | 31,262,246 | 339,721,971 | 46.397141 | 0.001424 |

## Interpretation

The D1 run does not weaken the handoff story. It makes the story more precise:

> the cardinality spectrum is large and smooth near the frontier, while the
> isolated structure appears in the zero-temperature field-selected ground
> face.

The honest theorem-language candidate remains:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe:

> PMF proves Sidon.

## Next Gate

Do not jump directly to public theorem language. The next useful gate is either:

```text
d = 2 near-ground probe on a single point, probably n = 70
```

or the stronger morphism test:

```text
port the same state-machine API to #755 or #166 and ask whether the same
field-sensitive face behavior survives a changed local exclusion rule
```

