# EHP114 n=14 Local Mixed-Remainder Scout

Experiment: `EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02`

## Meaning

The unrestricted admissible scout failed because the unit disk allows large
tangential excursions that are not local shape-cone perturbations. This repaired
scout adds the missing local condition: `||s|| <= 0.008`.

The claim ceiling remains narrow: this is a shadow signature, not universal
law. It is a finite local scout, not a uniform theorem.

## Verdict

- Status: `LOCAL_MIXED_SCOUT_PASS`
- Signed local points: `192`
- Local shape cap: `0.008`
- All points pass: `True`
- Minimum margin: `3.5420727715548637`
- Cap active count: `114`
- Next blocker: `Prove a uniform local cone radius, then bound mixed radial/shape remainder on that cone.`

## Per-Epsilon Summary

| eps | pass | min margin | cap active |
|---:|---:|---:|---:|
| 1e-04 | 48/48 | 3.54207277155 | 22 |
| 1e-03 | 48/48 | 5.02956427232 | 22 |
| 1e-02 | 48/48 | 6.93902629942 | 22 |
| 1e-01 | 48/48 | 8.87923940353 | 48 |

Worst point: `eps:1e-04:dir:21:-` with margin `3.5420727715548637`.

## What This Teaches

The local theorem needs two gates:

```text
roots stay in the unit disk
||s|| <= eta14
```

The first gate is admissibility. The second gate is locality. The predecessor
only used the first gate and therefore failed. This run uses both and gives the
next theorem a defensible domain.

## What Remains

Prove a uniform analytic version:

```text
R14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2
for 0 < eps <= 1e-1 and ||s|| <= min(eta14_boundary(eps), 0.008).
```

No scorecard, D1, public document, Lean file, git, or email state was changed.
