# Field-Handoff Witness Table, 51-71

Date: 2026-04-30
Problem: Erdos #30
Lane: exact maximizer face / field-selected top-k witnesses
Status: EXACT + INTERPRETIVE

## Question

Are the two first-hit failures, `n = 57` and `n = 70`, the same kind of
field-handoff event on the exact maximizer face?

## Answer

No. They now look like different handoff types.

- `n = 57` is a wide three-way split: prefix, mass, and joint winners are all
  far apart, and the best joint score is nonzero.
- `n = 70` is a wide split with a zero-joint witness: the first-hit scan failed,
  but the exact face still contains a witness satisfying both fields to
  numerical zero.

So `70` is now best understood as a first-hit selection artifact. `57` remains
the stronger field/optimizer handoff candidate.

## Evidence Packets

- `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`
- `EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29`

SHA256 sidecars verified `OK` during the 2026-04-30 continuation pass.

## Distance Definitions

For top-1 field-selected witnesses:

```text
p = prefix-field winner
m = density-adjusted-mass-field winner
j = joint-field winner
```

The table reports symmetric-difference distances:

```text
p-m d = |p triangle m|
p-j d = |p triangle j|
m-j d = |m triangle j|
```

All rows are exact maximizer-face data. The labels `wide split` and
`first-hit failure` are interpretive flags.

## Witness-Distance Table

| n | h | exact maximizers | split | p-m d | p-j d | m-j d | joint score | joint prefix loc | joint mass dev | flags |
|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|
| 51 | 9 | 2,690 | yes | 14 | 12 | 8 | 2.278e-16 | 51 | 0 | zero-joint; wide split |
| 52 | 9 | 5,550 | yes | 14 | 8 | 14 | 0.000e+0 | 52 | 0 | zero-joint; wide split |
| 53 | 9 | 11,260 | yes | 16 | 8 | 14 | 0.000e+0 | 53 | 0 | zero-joint; wide split |
| 54 | 9 | 21,164 | yes | 10 | 8 | 10 | 0.000e+0 | 54 | 0 | zero-joint |
| 55 | 10 | 2 | no | 0 | 0 | 0 | 6.069e-3 | 55 | 1.5 |  |
| 56 | 10 | 4 | no | 0 | 0 | 0 | 1.184e-2 | 56 | 3 |  |
| 57 | 10 | 6 | yes | 16 | 16 | 18 | 2.889e-2 | 57 | 7.5 | first-hit failure; wide split |
| 58 | 10 | 10 | yes | 18 | 18 | 0 | 3.616e-2 | 57 | 2 | wide split |
| 59 | 10 | 18 | yes | 16 | 14 | 16 | 1.653e-2 | 59 | 4.5 | wide split |
| 60 | 10 | 54 | yes | 14 | 14 | 0 | 0.000e+0 | 60 | 0 | zero-joint; wide split |
| 61 | 10 | 152 | yes | 18 | 16 | 16 | 1.754e-3 | 61 | 0.5 | wide split |
| 62 | 10 | 398 | yes | 12 | 12 | 0 | 0.000e+0 | 62 | 0 | zero-joint; wide split |
| 63 | 10 | 1,022 | yes | 12 | 12 | 0 | 1.678e-3 | 63 | 0.5 | wide split |
| 64 | 10 | 2,360 | yes | 10 | 14 | 14 | 0.000e+0 | 64 | 0 | zero-joint |
| 65 | 10 | 5,018 | yes | 12 | 10 | 10 | 1.608e-3 | 65 | 0.5 | wide split |
| 66 | 10 | 9,994 | yes | 12 | 10 | 14 | 1.818e-16 | 66 | 0 | zero-joint; wide split |
| 67 | 10 | 19,418 | yes | 16 | 14 | 14 | 1.542e-3 | 67 | 0.5 | wide split |
| 68 | 10 | 36,234 | yes | 16 | 12 | 12 | 1.328e-16 | 68 | 0 | zero-joint; wide split |
| 69 | 10 | 66,412 | yes | 14 | 4 | 14 | 1.481e-3 | 69 | 0.5 | wide split |
| 70 | 10 | 117,202 | yes | 14 | 12 | 14 | 1.295e-16 | 70 | 0 | first-hit failure; zero-joint; wide split |
| 71 | 10 | 203,840 | yes | 16 | 12 | 16 | 1.424e-3 | 71 | 0.5 | wide split |

## Interpretation

The `h = 10` transition at `n = 55` starts rigid: `55` and `56` have no
top-1 split. At `57`, the face opens into a true three-way field handoff. From
`58` onward, wide field splits are common, but zero-joint witnesses return at
many rows.

That makes `57` the more mathematically interesting handoff:

```text
57: wide split + no zero-joint rescue
70: wide split + zero-joint rescue
```

The D2 near-ground packet over `69 <= n <= 71` then explains why `70` should
not be treated as a raw spectral defect: the near-ground bath is smooth, and
the exact face contains a zero-joint witness.

## Next Gate

The next #30 test should focus on `56 <= n <= 59`, not `69 <= n <= 71`.

The right question is:

```text
Does the h=10 birth window create a genuine field-handoff onset at n=57?
```

Recommended run:

```text
EXP-MM-030-PMF-TRANSFER-PRUNED-D2-56-59-2026-04-30
```

with `--prune-deficiency 2`, `--frontier-k 5`, and the same packet contract.

Claim ceiling remains finite and interpretive. This would test the onset of
the handoff, not prove Sidon.
