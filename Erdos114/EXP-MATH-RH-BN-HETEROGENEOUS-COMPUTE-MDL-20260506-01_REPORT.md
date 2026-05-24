# EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01 Report

## Verdict

- Status: `STABLE_HETEROGENEOUS_COMPUTE_SIGNAL`
- Rows: `1440`
- Tolerance rows: `118`
- Heterogeneous reusable same-support wins: `288` / `408`
- Heterogeneous with-harness same-support wins: `48` / `408`
- Median reusable objective savings vs best baseline: `11.0`
- Median with-harness objective savings vs best baseline: `-40.0`

Tolerance winner counts:

```json
{
  "composite_index": 18,
  "heterogeneous_compute_reusable": 100
}
```

## Meaning

This experiment tests the TDP-inspired correction to the finite RH-MDL
cost model: when flat bits and factorization strings fail, charge
arithmetic operations as heterogeneous compute channels.

The result is a finite operation-weighted diagnostic. It is not a theorem
about primes, not a zeta formalization, and not a quantum-mechanical
claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite operation-weighted Beurling-Nyman MDL diagnostic; no RH claim, no zeta formalization, and no quantum claim.

## Objective

```text
total_objective_bits = compute_cost_bits - certified_residual_information_bits
```

## Tolerance Winners

| grid | dictionary | N | tolerance | winner | support | objective | residual |
|---|---|---:|---:|---|---|---:|---:|
| legendre_2048 | geometric | 8 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 226.182 | 0.1418 |
| legendre_2048 | geometric | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 143.596 | 0.189005 |
| legendre_2048 | geometric | 8 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 143.596 | 0.189005 |
| legendre_2048 | geometric | 16 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 304.797 | 0.108592 |
| legendre_2048 | geometric | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 221.081 | 0.132237 |
| legendre_2048 | geometric | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 221.081 | 0.132237 |
| legendre_2048 | geometric | 16 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 143.808 | 0.218866 |
| legendre_2048 | geometric | 24 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 309.629 | 0.0966558 |
| legendre_2048 | geometric | 24 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 309.629 | 0.0966558 |
| legendre_2048 | geometric | 24 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 226.162 | 0.139832 |
| legendre_2048 | geometric | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 226.162 | 0.139832 |
| legendre_2048 | geometric | 24 | 0.25 | composite_index | composite_prefix | 146.691 | 0.201757 |
| legendre_2048 | geometric | 32 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 394.061 | 0.065202 |
| legendre_2048 | geometric | 32 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 309.422 | 0.0837122 |
| legendre_2048 | geometric | 32 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 309.422 | 0.0837122 |
| legendre_2048 | geometric | 32 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 226.160 | 0.13967 |
| legendre_2048 | geometric | 32 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 226.160 | 0.13967 |
| legendre_2048 | geometric | 32 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 146.989 | 0.248031 |
| legendre_2048 | geometric | 48 | 0.06 | heterogeneous_compute_reusable | integer_prefix | 559.767 | 0.053166 |
| legendre_2048 | geometric | 48 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 394.300 | 0.0769641 |
| legendre_2048 | geometric | 48 | 0.10 | composite_index | composite_prefix | 327.589 | 0.094029 |
| legendre_2048 | geometric | 48 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 0.25 | composite_index | composite_prefix | 152.844 | 0.224411 |
| legendre_2048 | harmonic | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 127.614 | 0.191318 |
| legendre_2048 | harmonic | 8 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 127.614 | 0.191318 |
| legendre_2048 | harmonic | 16 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 289.155 | 0.139212 |
| legendre_2048 | harmonic | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 127.614 | 0.191318 |
| legendre_2048 | harmonic | 16 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 127.614 | 0.191318 |
| legendre_2048 | harmonic | 24 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 294.155 | 0.139212 |
| legendre_2048 | harmonic | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | harmonic | 24 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | harmonic | 32 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 551.933 | 0.119335 |
| legendre_2048 | harmonic | 32 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 294.155 | 0.139212 |
| legendre_2048 | harmonic | 32 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | harmonic | 32 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | harmonic | 48 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 703.845 | 0.11226 |
| legendre_2048 | harmonic | 48 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 294.155 | 0.139212 |
| legendre_2048 | harmonic | 48 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | harmonic | 48 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 130.614 | 0.191318 |
| legendre_2048 | seeded_log_uniform | 8 | 0.20 | heterogeneous_compute_reusable | integer_prefix | 258.268 | 0.150482 |
| legendre_2048 | seeded_log_uniform | 8 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 175.871 | 0.228666 |
| legendre_2048 | seeded_log_uniform | 16 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 420.573 | 0.0929477 |
| legendre_2048 | seeded_log_uniform | 16 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 336.750 | 0.105089 |
| legendre_2048 | seeded_log_uniform | 16 | 0.15 | composite_index | composite_prefix | 258.249 | 0.14851 |
| legendre_2048 | seeded_log_uniform | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 253.344 | 0.158624 |
| legendre_2048 | seeded_log_uniform | 16 | 0.25 | composite_index | composite_prefix | 178.845 | 0.224489 |
| legendre_2048 | seeded_log_uniform | 24 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 341.987 | 0.123851 |
| legendre_2048 | seeded_log_uniform | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 258.280 | 0.151826 |
| legendre_2048 | seeded_log_uniform | 24 | 0.25 | composite_index | composite_prefix | 178.939 | 0.239578 |
| legendre_2048 | seeded_log_uniform | 32 | 0.20 | heterogeneous_compute_reusable | integer_prefix | 262.672 | 0.199104 |
| legendre_2048 | seeded_log_uniform | 32 | 0.25 | composite_index | composite_prefix | 181.961 | 0.243265 |
| legendre_2048 | seeded_log_uniform | 48 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 592.177 | 0.0706367 |
| legendre_2048 | seeded_log_uniform | 48 | 0.10 | composite_index | composite_prefix | 359.626 | 0.0964507 |
| legendre_2048 | seeded_log_uniform | 48 | 0.12 | composite_index | composite_prefix | 359.626 | 0.0964507 |
| legendre_2048 | seeded_log_uniform | 48 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 337.952 | 0.120873 |
| legendre_2048 | seeded_log_uniform | 48 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 258.267 | 0.150379 |
| legendre_2048 | seeded_log_uniform | 48 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 258.267 | 0.150379 |
| midpoint_2048 | geometric | 8 | 0.15 | heterogeneous_compute_reusable | integer_prefix | 226.163 | 0.139997 |
| ... | ... | ... | ... | 58 more rows | ... | ... | ... |

## Boundary

A positive result here means only that this finite operation-weighted
cost model beats the listed finite baselines under the stated objective.
It does not establish an asymptotic law, dictionary invariance, or a
result about zeta zeros.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_RESULTS.json`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_REPORT.md`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.
