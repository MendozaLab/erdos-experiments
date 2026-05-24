# EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01 Report

## Verdict

- Status: `NO_PRIME_HARNESS_MDL_SIGNAL`
- Rows: `864`
- Same-support comparable rows: `144`
- Factorized reusable wins on same support: `12`
- Factorized with-harness wins on same support: `0`
- Median reusable savings vs flat direct index: `-28.0` bits
- Median charged savings vs flat direct index: `-91.5` bits

Tolerance winner counts:

```json
{
  "composite_index": 47,
  "factorized_reusable_harness": 20,
  "flat_index": 51
}
```

## Meaning

This experiment tests the revised prime-harness interpretation after the
composite-first projection diagnostic reversed. The question is whether
prime factorization works as a cheaper finite addressing layer for
Beurling-Nyman dictionary indices than flat integer addressing.

The result is an encoding-cost diagnostic. It is not a theorem about
primes, not a zeta formalization, and not a quantum-mechanical claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite Beurling-Nyman encoding-cost diagnostic; no RH claim, no zeta formalization, and no quantum claim.

## Cost Models

- `flat_index`: direct index address for selected dictionary columns.
- `composite_index`: unit plus composite-rank addressing only.
- `factorized_reusable_harness`: factorization addressing after the prime harness is already shared.
- `factorized_with_harness`: factorization addressing plus setup cost for the prime harness.

## Tolerance Winners

| grid | dictionary | N | tolerance | winner | support | bits |
|---|---|---:|---:|---|---|---:|
| legendre_2048 | geometric | 8 | 0.15 | flat_index | integer_prefix | 232 |
| legendre_2048 | geometric | 8 | 0.20 | flat_index | integer_prefix | 152 |
| legendre_2048 | geometric | 8 | 0.25 | flat_index | integer_prefix | 152 |
| legendre_2048 | geometric | 16 | 0.12 | flat_index | integer_prefix | 324 |
| legendre_2048 | geometric | 16 | 0.15 | flat_index | integer_prefix | 240 |
| legendre_2048 | geometric | 16 | 0.20 | composite_index | composite_prefix | 229 |
| legendre_2048 | geometric | 16 | 0.25 | composite_index | composite_prefix | 149 |
| legendre_2048 | geometric | 24 | 0.10 | flat_index | integer_prefix | 324 |
| legendre_2048 | geometric | 24 | 0.12 | flat_index | integer_prefix | 324 |
| legendre_2048 | geometric | 24 | 0.15 | flat_index | integer_prefix | 240 |
| legendre_2048 | geometric | 24 | 0.20 | composite_index | composite_prefix | 229 |
| legendre_2048 | geometric | 24 | 0.25 | composite_index | composite_prefix | 149 |
| legendre_2048 | geometric | 32 | 0.08 | flat_index | integer_prefix | 424 |
| legendre_2048 | geometric | 32 | 0.10 | flat_index | integer_prefix | 336 |
| legendre_2048 | geometric | 32 | 0.12 | composite_index | composite_prefix | 320 |
| legendre_2048 | geometric | 32 | 0.15 | composite_index | composite_prefix | 236 |
| legendre_2048 | geometric | 32 | 0.20 | composite_index | composite_prefix | 236 |
| legendre_2048 | geometric | 32 | 0.25 | composite_index | composite_prefix | 152 |
| legendre_2048 | geometric | 48 | 0.06 | flat_index | integer_prefix | 600 |
| legendre_2048 | geometric | 48 | 0.08 | flat_index | integer_prefix | 424 |
| legendre_2048 | geometric | 48 | 0.10 | composite_index | composite_prefix | 331 |
| legendre_2048 | geometric | 48 | 0.12 | factorized_reusable_harness | factor_cost_ordered | 250 |
| legendre_2048 | geometric | 48 | 0.15 | composite_index | composite_prefix | 243 |
| legendre_2048 | geometric | 48 | 0.20 | composite_index | composite_prefix | 243 |
| legendre_2048 | geometric | 48 | 0.25 | composite_index | composite_prefix | 155 |
| legendre_2048 | harmonic | 8 | 0.20 | factorized_reusable_harness | factor_cost_ordered | 139 |
| legendre_2048 | harmonic | 8 | 0.25 | flat_index | integer_prefix | 136 |
| legendre_2048 | harmonic | 16 | 0.15 | flat_index | integer_prefix | 308 |
| legendre_2048 | harmonic | 16 | 0.20 | factorized_reusable_harness | factor_cost_ordered | 139 |
| legendre_2048 | harmonic | 16 | 0.25 | factorized_reusable_harness | factor_cost_ordered | 139 |
| legendre_2048 | harmonic | 24 | 0.15 | flat_index | integer_prefix | 308 |
| legendre_2048 | harmonic | 24 | 0.20 | factorized_reusable_harness | factor_cost_ordered | 142 |
| legendre_2048 | harmonic | 24 | 0.25 | flat_index | integer_prefix | 140 |
| legendre_2048 | harmonic | 32 | 0.12 | factorized_reusable_harness | factor_cost_ordered | 641 |
| legendre_2048 | harmonic | 32 | 0.15 | flat_index | integer_prefix | 320 |
| legendre_2048 | harmonic | 32 | 0.20 | factorized_reusable_harness | factor_cost_ordered | 142 |
| legendre_2048 | harmonic | 32 | 0.25 | factorized_reusable_harness | factor_cost_ordered | 142 |
| legendre_2048 | harmonic | 48 | 0.12 | flat_index | integer_prefix | 760 |
| legendre_2048 | harmonic | 48 | 0.15 | flat_index | integer_prefix | 320 |
| legendre_2048 | harmonic | 48 | 0.20 | factorized_reusable_harness | factor_cost_ordered | 142 |
| ... | ... | ... | ... | 78 more rows | ... | ... |

## Interpretation Boundary

If the reusable factorized code wins but the charged code loses, the
finite message is amortization: the prime-address layer helps only
when the harness is already part of the shared dictionary protocol.
That is useful MDL structure, but not a free-description claim.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01_RESULTS.json`
- `EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01_REPORT.md`
- `EXP-MATH-RH-BN-PRIME-HARNESS-MDL-20260506-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.
