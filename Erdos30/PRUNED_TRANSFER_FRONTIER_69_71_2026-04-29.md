# Pruned Transfer Frontier Scale Gate, 69-71

Date: 2026-04-29
Problem: Erdős #30
Lane: PMF transfer operator / lattice-gas ground-state face
Status: SCALE_GATE_PASS for top-1 frontier parity on `69 <= n <= 71`

## Question

Can the transfer operator reproduce the `69-71` top-k frontier behavior without
enumerating every low-cardinality Sidon state?

## Answer

Yes, for the tested top-1 prefix, mass, and joint frontier witnesses.

The first full HashMap transfer operator was too blunt at this scale: even with
reachability pruning, it still carried too much prefix-state mass. The working
version is stricter. It uses the same transfer state,

```text
State = {
  occupied_mask: u128,
  used_differences_mask: u128,
  cardinality: u8
}
```

but traverses depth-first with a layer gate:

```text
keep only branches that can still reach h(n)-d
```

For this scale gate, `d = 0`, so terminal states are retained only if they reach
the exact ground layer `h(n) = 10`.

## Packets

Transfer packets:

- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-69-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-70-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29`

Exact reference packets:

- `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`
- `EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29`

All three transfer packet SHA256 sidecars verified `OK`.

## Result Summary

| n | h(n) | exact maximizers retained | pruned states | transitions attempted | runtime sec | joint score |
|---|---:|---:|---:|---:|---:|---:|
| 69 | 10 | 66,412 | 218,433,851 | 1,280,651,251 | 28.298015 | 0.001481 |
| 70 | 10 | 117,202 | 260,934,836 | 1,541,652,499 | 43.711637 | ~0 |
| 71 | 10 | 203,840 | 311,297,246 | 1,853,066,947 | 53.424933 | 0.001424 |

The retained terminal count equals the exact maximizer count in every row. That
is the key scale result: the transfer traversal did not enumerate every
low-cardinality terminal state.

## Frontier Parity

Machine comparison of transfer top-1 prefix, mass, and joint witnesses against
the exact Rust maximizer packets produced an empty diff for `69 <= n <= 71`.

| n | prefix-field winner | mass-field winner | joint-field winner |
|---|---|---|---|
| 69 | `[0, 3, 11, 20, 33, 43, 62, 64, 68, 69]` | `[0, 1, 5, 28, 43, 49, 57, 60, 67, 69]` | `[0, 3, 11, 21, 33, 49, 62, 64, 68, 69]` |
| 70 | `[0, 4, 12, 21, 32, 39, 55, 65, 68, 70]` | `[0, 1, 5, 28, 44, 50, 58, 61, 68, 70]` | `[0, 4, 12, 31, 44, 46, 49, 60, 69, 70]` |
| 71 | `[0, 4, 13, 21, 40, 45, 56, 68, 70, 71]` | `[0, 1, 6, 31, 44, 53, 55, 63, 67, 70]` | `[0, 4, 13, 23, 34, 51, 63, 65, 66, 71]` |

## Interpretation

This passes the scale gate we needed. The PMF transfer operator can reproduce
the exact frontier behavior around `69-71` without materializing the full
low-cardinality lattice gas.

The honest theorem-language candidate is now:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe:

> PMF proves Sidon.

