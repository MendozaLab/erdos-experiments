# EHP114 n=14 Admissible Mixed-Remainder Scout

Experiment: `EXP-MATH-EHP114-N14-ADMISSIBLE-MIXED-REMAINDER-SCOUT-20260505-01`

## Meaning

The radial and shape pieces are now separately certified. This scout checks
whether the first boundary-aware mixed target survives finite admissible tests:
roots are moved radially inward, then perturbed only as far as the unit disk
allows.

The claim ceiling remains narrow: this is a shadow signature, not universal
law. It is not a uniform theorem.

## Verdict

- Status: `FINITE_MIXED_SCOUT_FAIL`
- Signed admissible points: `192`
- All points pass: `False`
- Minimum margin: `-8460.077650590973`
- Empirical min `tmax/eps`: `0.1889910112460592`
- Next blocker: `Turn the empirical admissible cone radius and finite mixed scout into a uniform analytic remainder bound.`

## Per-Epsilon Summary

| eps | pass | min margin | min tmax | max tmax |
|---:|---:|---:|---:|---:|
| 1e-04 | 46/48 | -1.52040374428 | 1.88991011246e-05 | 0.014142438685 |
| 1e-03 | 26/48 | -73.2794527346 | 0.000189070034704 | 0.0447309476073 |
| 1e-02 | 26/48 | -804.430756136 | 0.00189865336196 | 0.141725963456 |
| 1e-01 | 0/48 | -8460.07765059 | 0.0198365298716 | 0.457321685277 |

Worst point: `eps:1e-01:dir:23:+` with margin `-8460.077650590973`.

## Interpretation

The shape-cone theorem should not allow arbitrary shape amplitude independent
of the radial slack. The unit-disk boundary imposes the real cone:

```text
||s|| <= eta14(eps).
```

This scout suggests that an `eta14(eps)` linear in `eps` is the safe first
target. That is weaker than the abstract tangent cone, but it is the honest
boundary-compatible lane toward local stability.

## What Remains

Turn the finite scout into a uniform analytic bound:

```text
R14(eps,s) <= 12 eps^(1/14) + 50000 ||s||^2
for 0 < eps <= 1e-1 and ||s|| <= eta14(eps).
```

No scorecard, D1, public document, Lean file, git, or email state was changed.
