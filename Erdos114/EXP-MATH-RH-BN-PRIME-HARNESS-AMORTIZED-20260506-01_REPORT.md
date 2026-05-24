# EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01 Report

## Verdict

- Status: `NO_AMORTIZED_PRIME_HARNESS_SIGNAL`
- Rows: `2040`
- Bundles: `18`
- Shared-harness wins: `8`
- Reusable-without-setup positive bundles: `8`
- Median shared-harness savings vs flat: `-80.5` bits
- Median reusable savings before setup: `-21.0` bits

Best shared-harness bundle:

```json
{
  "N": 48,
  "common_task_count": 16,
  "dictionary": "seeded_log_uniform",
  "factor_shared_savings_vs_flat": 1714,
  "grid": "legendre_2048",
  "harness_setup_bits": 90
}
```

## Meaning

This experiment tests whether a prime-address harness becomes useful when
its setup cost is shared across a bundle of related finite approximation
tasks. The targets are simple weighted functions over the same
fractional-part dictionary, including the standard constant target.

The result is still finite encoding accounting. It is not a theorem about
primes, not a zeta formalization, and not a quantum-mechanical claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite amortized encoding-cost diagnostic; no RH claim, no zeta formalization, and no quantum claim.

## Bundle Results

| grid | dictionary | N | tasks | harness bits | flat bits | shared factor bits | savings | win |
|---|---|---:|---:|---:|---:|---:|---:|---|
| legendre_2048 | geometric | 16 | 8 | 30 | 2508 | 2390 | 118 | True |
| legendre_2048 | geometric | 32 | 15 | 66 | 7592 | 6791 | 801 | True |
| legendre_2048 | geometric | 48 | 18 | 90 | 11328 | 9722 | 1606 | True |
| legendre_2048 | harmonic | 16 | 8 | 30 | 2296 | 2379 | -83 | False |
| legendre_2048 | harmonic | 32 | 12 | 66 | 5072 | 5351 | -279 | False |
| legendre_2048 | harmonic | 48 | 12 | 90 | 5072 | 5643 | -571 | False |
| legendre_2048 | seeded_log_uniform | 16 | 7 | 30 | 2576 | 2665 | -89 | False |
| legendre_2048 | seeded_log_uniform | 32 | 2 | 66 | 560 | 638 | -78 | False |
| legendre_2048 | seeded_log_uniform | 48 | 16 | 90 | 11784 | 10070 | 1714 | True |
| midpoint_2048 | geometric | 16 | 8 | 30 | 2508 | 2281 | 227 | True |
| midpoint_2048 | geometric | 32 | 15 | 66 | 7592 | 7253 | 339 | True |
| midpoint_2048 | geometric | 48 | 19 | 90 | 12632 | 11075 | 1557 | True |
| midpoint_2048 | harmonic | 16 | 8 | 30 | 2296 | 2379 | -83 | False |
| midpoint_2048 | harmonic | 32 | 12 | 66 | 5072 | 5594 | -522 | False |
| midpoint_2048 | harmonic | 48 | 12 | 90 | 5072 | 5643 | -571 | False |
| midpoint_2048 | seeded_log_uniform | 16 | 8 | 30 | 3016 | 3147 | -131 | False |
| midpoint_2048 | seeded_log_uniform | 32 | 2 | 66 | 648 | 744 | -96 | False |
| midpoint_2048 | seeded_log_uniform | 48 | 16 | 90 | 11432 | 9978 | 1454 | True |

## Interpretation Boundary

A win here would mean only that this finite cost model rewards a shared
factorization-address layer across several approximation tasks. A loss
means the naive prime-harness code is still too expensive under this
finite accounting model.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01_RESULTS.json`
- `EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01_REPORT.md`
- `EXP-MATH-RH-BN-PRIME-HARNESS-AMORTIZED-20260506-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.
