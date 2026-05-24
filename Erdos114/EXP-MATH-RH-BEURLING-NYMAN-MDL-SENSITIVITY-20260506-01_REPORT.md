# EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01 Report

## Claim Tested

The prior visual showed a bend near `-log2(relative residual) ~= 3.2`.
This sensitivity run asks whether that bend is stable across dictionary
families, seeded log-uniform replicates, and a denser `N` sweep.

## Claim Ceiling

SUGGESTIVE_NUMERICAL_SENSITIVITY: finite quadrature sensitivity test of the ~3.2-bit bend; not RH evidence, not theorem progress, and not an asymptotic statement.

No D1, scorecard, Lean status, public page, publication surface, staging, or
commit is updated by this run.

## Method

- Source experiment: `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-03`
- Quadrature: 2048 Gauss-Legendre nodes on `[0,1]`
- Basis: `rho_theta(x) = fractional_part(theta / x)`
- Families: harmonic, geometric, and 6 seeded log-uniform schedules
- Canonical curve for sensitivity: 16-bit coefficient quantization
- Honest y-axis: `-log2(certified_upper_relative_residual)`, not the unpenalized residual

## Result

The `3.2`-bit feature is not yet a universal constant. In this finite sweep, it
is better described as the first conditioning crossover of the tested finite
families. The geometric dictionary is often best, while seeded schedules can
briefly win around intermediate `N`. The largest condition jumps cluster around
the same information range where the visual bend appeared.

## Best Certified Curve By N

| N | best series | certified info bits | certified upper residual | condition | description bits |
|---:|---|---:|---:|---:|---:|
| 4 | geometric | 2.3492 | 0.196252 | 4.985 | 129 |
| 8 | seeded_log_uniform_r4 | 2.9449 | 0.129869 | 13.19 | 241 |
| 12 | seeded_log_uniform_r4 | 3.1566 | 0.112142 | 19.23 | 305 |
| 16 | seeded_log_uniform_r0 | 3.4274 | 0.0929477 | 41.55 | 369 |
| 20 | geometric | 3.6496 | 0.0796797 | 24.5 | 385 |
| 24 | geometric | 3.4423 | 0.0919961 | 36.39 | 449 |
| 28 | geometric | 3.9231 | 0.0659223 | 37.5 | 513 |
| 32 | geometric | 4.0227 | 0.0615262 | 44.67 | 577 |
| 40 | seeded_log_uniform_r4 | 4.3494 | 0.0490562 | 210.7 | 753 |
| 48 | geometric | 4.2612 | 0.0521498 | 73.31 | 833 |

## Rows Crossing The 3.2-Bit Band

| series | N step | certified info bits | marginal info / description bit | condition |
|---|---:|---:|---:|---:|
| geometric | 12->16 | 3.0858->3.2046 | 0.00185648 | 14.41->22.61 |
| harmonic | 40->48 | 3.1807->3.2296 | 0.000382302 | 66.46->85.1 |
| seeded_log_uniform_r0 | 12->16 | 2.4028->3.4274 | 0.0160093 | 26.65->41.55 |
| seeded_log_uniform_r0 | 20->24 | 3.2044->3.0391 | -0.00258376 | 39.74->43.13 |
| seeded_log_uniform_r0 | 32->40 | 2.5168->3.3221 | 0.00629137 | 111.8->1248 |
| seeded_log_uniform_r1 | 24->28 | 3.1024->3.3053 | 0.00316999 | 76.42->315.9 |
| seeded_log_uniform_r2 | 24->28 | 2.9035->3.4356 | 0.00831374 | 32.58->64.02 |
| seeded_log_uniform_r2 | 28->32 | 3.4356->2.8119 | -0.00974448 | 64.02->74.02 |
| seeded_log_uniform_r2 | 32->40 | 2.8119->3.6467 | 0.00652165 | 74.02->80.86 |
| seeded_log_uniform_r3 | 20->24 | 2.9703->3.3569 | 0.00603944 | 83.91->83.86 |
| seeded_log_uniform_r3 | 24->28 | 3.3569->2.8696 | -0.00761344 | 83.86->42.04 |
| seeded_log_uniform_r3 | 28->32 | 2.8696->3.3059 | 0.00681732 | 42.04->943.7 |
| seeded_log_uniform_r4 | 16->20 | 2.5032->3.2022 | 0.0109219 | 63.82->94.64 |
| seeded_log_uniform_r5 | 20->24 | 3.0293->3.2260 | 0.00307289 | 38.07->59.05 |
| seeded_log_uniform_r5 | 40->48 | 3.2179->3.1105 | -0.000839309 | 279.8->305.8 |

## Largest Condition Jumps

| series | N step | condition ratio | certified info bits | marginal info / description bit |
|---|---:|---:|---:|---:|
| seeded_log_uniform_r3 | 28->32 | 22.45 | 2.8696->3.3059 | 0.00681732 |
| seeded_log_uniform_r1 | 40->48 | 16.54 | 3.8603->3.5649 | -0.00230807 |
| seeded_log_uniform_r0 | 32->40 | 11.16 | 2.5168->3.3221 | 0.00629137 |
| seeded_log_uniform_r5 | 32->40 | 5.687 | 3.4636->3.2179 | -0.00191927 |
| seeded_log_uniform_r3 | 16->20 | 4.606 | 3.0784->2.9703 | -0.00168819 |
| seeded_log_uniform_r1 | 24->28 | 4.134 | 3.1024->3.3053 | 0.00316999 |

## Interpretation

The bend has a plausible operational meaning: after roughly three bits of
certified residual information, the finite dictionaries begin paying a larger
conditioning tax. That does not make `3.2` a mathematical constant. It makes it
a candidate finite-N diagnostic to stress-test.

Safe phrasing:

```text
first observed finite-N MDL conditioning crossover
```

Unsafe phrasing:

```text
new information-theoretic constant
```

## Next Work

If this remains interesting, the next run should vary quadrature resolution and
dictionary parametrization while preserving the same bit-cost model. A true
phenomenon should survive those changes; a numerical artifact will move.
