# EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01 Report

## Verdict

- Status: `HARDENED_HETEROGENEOUS_COMPUTE_SIGNAL`
- Rows: `77112`
- Tolerance rows: `13853`
- Base reusable tolerance wins: `1546` / `1979`
- Base same-support wins: `2307` / `2754`
- Base median reusable savings vs best baseline: `25.0`
- Implementation audit: `PASS`

Classification reasons:

```json
[
  "base_passes=True; profile_pass_count=4/7",
  "audit_status=PASS; audit_error_count=0"
]
```

## Meaning

This hardening gate tests whether the heterogeneous-compute signal
survives regularized active-set supports, operation-weight sensitivity,
and an independent implementation audit.

The result is still a finite operation-weighted diagnostic. It is not
a theorem about primes, not a zeta formalization, and not a quantum
claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite heterogeneous-compute Beurling-Nyman MDL hardening gate; no RH claim, no zeta formalization, and no quantum claim.

## Weight Profiles

| profile | tolerance wins | tolerance rows | same-support wins | comparable | median savings |
|---|---:|---:|---:|---:|---:|
| `all_unit` | `1895` | `1979` | `2589` | `2754` | `42.0` |
| `base` | `1546` | `1979` | `2307` | `2754` | `25.0` |
| `free_arithmetic` | `1979` | `1979` | `2754` | `2754` | `96.0` |
| `high_ops` | `557` | `1979` | `636` | `2754` | `-27.0` |
| `low_ops` | `1925` | `1979` | `2643` | `2754` | `45.0` |
| `operator_heavy` | `377` | `1979` | `432` | `2754` | `-72.0` |
| `prime_heavy` | `0` | `1979` | `0` | `2754` | `-47.0` |

## Implementation Audit

```json
{
  "audit_status": "PASS",
  "checked_rows": 77112,
  "error_count": 0,
  "errors": []
}
```

## Base Tolerance Winners

| support | grid | dictionary | N | coeff bits | tolerance | winner | objective | residual |
|---|---|---|---:|---:|---:|---|---:|---:|
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 8 | 0.12 | heterogeneous_compute_reusable | 208.937 | 0.119689 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 8 | 0.15 | heterogeneous_compute_reusable | 162.001 | 0.125109 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 8 | 0.20 | flat_index | 121.548 | 0.182818 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 8 | 0.25 | flat_index | 121.548 | 0.182818 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 16 | 0.12 | heterogeneous_compute_reusable | 225.907 | 0.117188 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 16 | 0.15 | heterogeneous_compute_reusable | 225.907 | 0.117188 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 16 | 0.20 | flat_index | 153.515 | 0.178657 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 16 | 0.25 | flat_index | 153.515 | 0.178657 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 24 | 0.12 | heterogeneous_compute_reusable | 289.907 | 0.117157 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 24 | 0.15 | heterogeneous_compute_reusable | 289.907 | 0.117157 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 24 | 0.20 | flat_index | 185.515 | 0.17864 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 16 | 24 | 0.25 | flat_index | 185.515 | 0.17864 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.08 | heterogeneous_compute_reusable | 257.338 | 0.0789759 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.10 | heterogeneous_compute_reusable | 211.499 | 0.0883485 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.12 | heterogeneous_compute_reusable | 165.887 | 0.115551 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.15 | heterogeneous_compute_reusable | 165.887 | 0.115551 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.20 | heterogeneous_compute_reusable | 119.662 | 0.19776 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 8 | 0.25 | heterogeneous_compute_reusable | 119.662 | 0.19776 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.08 | heterogeneous_compute_reusable | 307.285 | 0.076155 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.10 | heterogeneous_compute_reusable | 307.285 | 0.076155 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.12 | heterogeneous_compute_reusable | 229.782 | 0.107437 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.15 | heterogeneous_compute_reusable | 229.782 | 0.107437 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.20 | heterogeneous_compute_reusable | 151.631 | 0.193608 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 16 | 0.25 | heterogeneous_compute_reusable | 151.631 | 0.193608 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.08 | heterogeneous_compute_reusable | 403.284 | 0.0761073 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.10 | heterogeneous_compute_reusable | 403.284 | 0.0761073 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.12 | heterogeneous_compute_reusable | 293.781 | 0.107405 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.15 | heterogeneous_compute_reusable | 293.781 | 0.107405 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.20 | heterogeneous_compute_reusable | 183.631 | 0.193592 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 32 | 24 | 0.25 | heterogeneous_compute_reusable | 183.631 | 0.193592 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.08 | heterogeneous_compute_reusable | 259.349 | 0.0796126 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.10 | heterogeneous_compute_reusable | 215.541 | 0.0909499 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.12 | heterogeneous_compute_reusable | 162.861 | 0.113507 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.15 | heterogeneous_compute_reusable | 162.861 | 0.113507 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.20 | heterogeneous_compute_reusable | 119.554 | 0.1835 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 8 | 0.25 | heterogeneous_compute_reusable | 119.554 | 0.1835 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.06 | heterogeneous_compute_reusable | 568.754 | 0.0526982 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.08 | heterogeneous_compute_reusable | 311.331 | 0.0786411 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.10 | heterogeneous_compute_reusable | 311.331 | 0.0786411 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.12 | heterogeneous_compute_reusable | 226.754 | 0.105386 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.15 | heterogeneous_compute_reusable | 226.754 | 0.105386 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.20 | heterogeneous_compute_reusable | 151.520 | 0.179231 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 16 | 0.25 | heterogeneous_compute_reusable | 151.520 | 0.179231 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.06 | heterogeneous_compute_reusable | 760.751 | 0.0526099 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.08 | heterogeneous_compute_reusable | 407.331 | 0.078593 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.10 | heterogeneous_compute_reusable | 407.331 | 0.078593 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.12 | heterogeneous_compute_reusable | 290.753 | 0.105354 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.15 | heterogeneous_compute_reusable | 290.753 | 0.105354 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.20 | heterogeneous_compute_reusable | 183.520 | 0.179214 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 48 | 24 | 0.25 | heterogeneous_compute_reusable | 183.520 | 0.179214 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.08 | heterogeneous_compute_reusable | 269.272 | 0.0754536 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.10 | heterogeneous_compute_reusable | 220.514 | 0.0892516 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.12 | heterogeneous_compute_reusable | 181.897 | 0.116362 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.15 | heterogeneous_compute_reusable | 181.897 | 0.116362 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.20 | flat_index | 129.535 | 0.181112 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 8 | 0.25 | flat_index | 129.535 | 0.181112 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.06 | heterogeneous_compute_reusable | 396.925 | 0.0593495 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.08 | heterogeneous_compute_reusable | 316.307 | 0.0773152 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.10 | heterogeneous_compute_reusable | 316.307 | 0.0773152 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.12 | heterogeneous_compute_reusable | 245.791 | 0.108137 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.15 | heterogeneous_compute_reusable | 245.791 | 0.108137 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.20 | flat_index | 161.501 | 0.17688 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 16 | 0.25 | flat_index | 161.501 | 0.17688 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.06 | heterogeneous_compute_reusable | 524.924 | 0.0592865 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.08 | heterogeneous_compute_reusable | 412.306 | 0.0772685 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.10 | heterogeneous_compute_reusable | 412.306 | 0.0772685 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.12 | heterogeneous_compute_reusable | 309.791 | 0.108105 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.15 | heterogeneous_compute_reusable | 309.791 | 0.108105 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.20 | flat_index | 193.501 | 0.176863 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 64 | 24 | 0.25 | flat_index | 193.501 | 0.176863 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.08 | heterogeneous_compute_reusable | 279.320 | 0.0780384 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.10 | heterogeneous_compute_reusable | 238.455 | 0.0857029 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.12 | flat_index | 188.901 | 0.116707 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.15 | flat_index | 188.901 | 0.116707 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.20 | flat_index | 129.531 | 0.180564 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 8 | 0.25 | flat_index | 129.531 | 0.180564 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 16 | 0.06 | heterogeneous_compute_reusable | 566.451 | 0.0427082 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 16 | 0.08 | heterogeneous_compute_reusable | 334.240 | 0.0737868 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 16 | 0.10 | heterogeneous_compute_reusable | 334.240 | 0.0737868 |
| omp_ridge_alpha_1e-06 | legendre_2048 | geometric | 96 | 16 | 0.12 | flat_index | 252.796 | 0.108536 |
| ... | ... | ... | ... | ... | ... | 1899 more rows | ... | ... |

## Boundary

A hardened result here means only that the finite reusable
operation-weighted cost model survived this stricter finite gate.
It does not establish an infinite-dimensional theorem, dictionary
invariance, or a result about zeta zeros.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01_RESULTS.json`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01_REPORT.md`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.
