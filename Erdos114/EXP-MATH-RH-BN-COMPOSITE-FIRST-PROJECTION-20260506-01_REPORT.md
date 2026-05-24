# EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01 Report

## Verdict

- Status: `STABLE_FINITE_PROJECTION_ASYMMETRY`
- Rows: `24`
- Dominant order: `prime_first`
- Order counts: `{'composite_first': 5, 'prime_first': 19, 'tie': 0}`
- Standalone bulk counts: `{'composite': 6, 'prime': 18, 'tie': 0}`
- Median absolute order gap in relative residual: `0.0387626`

## Meaning

This experiment tests a finite projection model for the RH-MDL lane. The
question is whether composite-indexed Beurling-Nyman columns absorb bulk
residual first, while prime-indexed columns supply harder marginal directions.

The result is a finite linear-algebra diagnostic. It is not a zeta
formalization, not a theorem about primes, and not a public RH claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite Beurling-Nyman projection diagnostic; no RH claim, no zeta formalization, and no asymptotic claim.

## Model

- Target: `chi_(0,1] = 1` on `(0,1]`.
- Basis: `rho_theta(x) = fractional_part(theta / x)`.
- Column index split: `1` is unit, prime indices are prime columns, composite
  indices greater than `1` are composite columns.
- Grids: `legendre_2048, midpoint_2048`.
- Dictionary sizes: `[8, 16, 24, 32]`.
- Dictionaries: `harmonic, geometric, seeded_log_uniform`.

## Selected Geometric Rows

| grid | N | order winner | composite-then-prime residual | prime-then-composite residual | all-column residual | composite-minus-prime gap |
|---|---:|---|---:|---:|---:|---:|
| legendre_2048 | 8 | prime_first | 0.228306 | 0.168456 | 0.141771 | +5.985e-02 |
| legendre_2048 | 16 | prime_first | 0.162212 | 0.152112 | 0.108417 | +1.010e-02 |
| legendre_2048 | 32 | composite_first | 0.10448 | 0.122577 | 0.0614301 | -1.810e-02 |
| midpoint_2048 | 8 | prime_first | 0.227404 | 0.16766 | 0.139968 | +5.974e-02 |
| midpoint_2048 | 16 | prime_first | 0.162062 | 0.15291 | 0.108931 | +9.152e-03 |
| midpoint_2048 | 32 | composite_first | 0.103949 | 0.121876 | 0.0610429 | -1.793e-02 |

## Stability Checks

The protocol required the order asymmetry to be stable across at least three
dictionary sizes, two quadrature grids, and both orderings. The observed status
is:

```text
STABLE_FINITE_PROJECTION_ASYMMETRY
```

Grid counts:

```json
{
  "legendre_2048": {
    "composite_first": 2,
    "prime_first": 10,
    "tie": 0
  },
  "midpoint_2048": {
    "composite_first": 3,
    "prime_first": 9,
    "tie": 0
  }
}
```

N counts:

```json
{
  "8": {
    "composite_first": 0,
    "prime_first": 6,
    "tie": 0
  },
  "16": {
    "composite_first": 2,
    "prime_first": 4,
    "tie": 0
  },
  "24": {
    "composite_first": 0,
    "prime_first": 6,
    "tie": 0
  },
  "32": {
    "composite_first": 3,
    "prime_first": 3,
    "tie": 0
  }
}
```

## Spectral Diagnostic Boundary

Each row stores finite projection-spectrum data for the second-stage residual:
capturable energy, the fraction of residual energy captured, and top-mode
concentration. This is only a finite matrix diagnostic. It should not be
identified with the Riemann zeta function or with a Hilbert-Polya operator.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01_RESULTS.json`
- `EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01_REPORT.md`
- `EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher
surface was updated.
