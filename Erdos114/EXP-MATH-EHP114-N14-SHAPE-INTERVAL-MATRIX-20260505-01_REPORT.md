# EHP114 n=14 Shape Interval Matrix Certificate

Experiment: `EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01`

## Meaning

This hardens the n=14 shape cone from a floating matrix into an interval matrix
certificate, using the Rust/inari length oracle. It is still local and finite:
the certified object is the finite-difference shape matrix at eps
`0.005`, not the full EHP114 theorem.

The claim ceiling remains: this is a shadow signature, not universal law.

## Verdict

- Status: `SHAPE_INTERVAL_MATRIX_CERTIFIED`
- Shape basis rank: `24`
- Oracle points: `1152`
- Gershgorin interval lower bound: `308540.9038178724`
- Working lambda14: `100000.0`
- Lambda below interval bound: `True`
- Max matrix interval width: `8.693490372024826e-08`
- Next blocker: `Prove mixed radial/shape remainder absorption for the n=14 local cone.`

## Row Bounds

| row | diagonal lower | offdiag radius upper | lower minus radius |
|---:|---:|---:|---:|
| 0 | 702839.628467 | 44712.2577564 | 658127.37071 |
| 1 | 702865.833598 | 44780.296767 | 658085.536831 |
| 2 | 610132.782779 | 300826.172055 | 309306.610724 |
| 3 | 621309.352646 | 312196.098738 | 309113.253908 |
| 4 | 621305.791265 | 312193.757553 | 309112.033711 |
| 5 | 609396.377743 | 300855.473925 | 308540.903818 |
| 6 | 565643.157147 | 253373.702281 | 312269.454866 |
| 7 | 565678.994156 | 253368.270785 | 312310.723371 |
| 8 | 566384.044908 | 253092.072053 | 313291.972855 |
| 9 | 566319.778809 | 255009.855285 | 311309.923524 |
| 10 | 526863.555521 | 161255.134361 | 365608.42116 |
| 11 | 521434.355241 | 142104.473511 | 379329.881731 |
| 12 | 521488.679484 | 141945.273581 | 379543.405903 |
| 13 | 526896.06104 | 160199.078383 | 366696.982657 |
| 14 | 498252.945094 | 96409.2129184 | 401843.732176 |
| 15 | 498254.836353 | 95282.3151429 | 402972.521211 |
| 16 | 494713.613587 | 87550.0482116 | 407163.565375 |
| 17 | 494701.23733 | 88653.5228107 | 406047.71452 |
| 18 | 478822.204853 | 27636.4818757 | 451185.722978 |
| 19 | 471816.770445 | 19520.7753316 | 452295.995113 |
| 20 | 471800.149086 | 19592.0448985 | 452208.104187 |
| 21 | 478782.146915 | 28662.2464298 | 450119.900485 |
| 22 | 474536.097863 | 19911.1352573 | 454624.962605 |
| 23 | 474536.359779 | 18031.7974796 | 456504.562299 |

## What This Changes

The shape cone is no longer only a floating diagnostic. Under the same
interval marching-squares method used in the Rust/inari finite certificate, the
n=14 shape matrix has a positive interval Gershgorin lower bound.

This makes the Jordan-style articulation sharper: the shape directions behave
like a positive spectral cone, while the radial direction is the Puiseux
boundary term.

## What Remains

Bound mixed radial/shape remainder by <= 12*eps^(1/14) + 0.5*lambda14*||s||^2 on the n=14 local cone.

That mixed-remainder theorem is now the live blocker.

## Files

- Script: `erdos-experiments/Erdos114/ehp114_n14_shape_interval_matrix_certificate.py`
- Rust oracle: `erdos-experiments/scripts/erdos-114/src/bin/ehp114_batch_interval_lengths.rs`
- Oracle input: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01_ORACLE_INPUT.json`
- Oracle output: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01_ORACLE_OUTPUT.json`
- Result SHA: `6b686b0346bceaa9db7ec9225232c4f1748366367fa061d75c005db6cd281005` for oracle output

No scorecard, D1, public document, Lean file, git, or email state was changed.
