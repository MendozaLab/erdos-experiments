# EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01 Report

## Verdict

- Status: `REPLICATED_HETEROGENEOUS_COMPUTE_SIGNAL`
- Rows: `40014`
- Tolerance rows: `3745`
- Base heterogeneous tolerance wins: `634` / `749`
- Base same-support wins: `2268` / `3051`
- Base median objective savings vs best baseline: `20.0`

Classification reasons:

```json
[
  "base_replicates=True; ablation_pass_count=4/4"
]
```

## Meaning

This replication tries to break the positive heterogeneous-compute
signal by adding larger dictionary sizes, a shifted grid, coefficient
precision sensitivity, and one-channel-at-a-time ablations.

The result is still a finite operation-weighted diagnostic. It is not
a theorem about primes, not a zeta formalization, and not a
quantum-mechanical claim.

## Claim Ceiling

INTERNAL / METHOD-SHAPING ONLY: finite heterogeneous-compute Beurling-Nyman MDL replication gate; no RH claim, no zeta formalization, and no quantum claim.

## Base Scenario

```json
{
  "heterogeneous_encoding": "heterogeneous_compute_reusable",
  "heterogeneous_same_support_win_count": 2268,
  "heterogeneous_same_support_win_fraction": 0.7433628318584071,
  "heterogeneous_tolerance_win_count": 634,
  "heterogeneous_tolerance_win_fraction": 0.8464619492656876,
  "heterogeneous_wins_by_N": {
    "16": 89,
    "24": 78,
    "32": 90,
    "48": 114,
    "64": 108,
    "96": 155
  },
  "heterogeneous_wins_by_coefficient_bits": {
    "16": 230,
    "24": 229,
    "8": 175
  },
  "heterogeneous_wins_by_grid": {
    "legendre_2048": 207,
    "midpoint_2048": 212,
    "shifted_midpoint_2048": 215
  },
  "median_objective_savings_vs_best_baseline": 20.0,
  "same_support_comparable_count": 3051,
  "scenario": "base",
  "tolerance_row_count": 749,
  "winner_counts": {
    "composite_index": 115,
    "heterogeneous_compute_reusable": 634
  }
}
```

## Ablation Scenarios

| scenario | tolerance wins | tolerance rows | same-support wins | comparable | median savings |
|---|---:|---:|---:|---:|---:|
| `drop_exponent` | `692` | `749` | `2754` | `3051` | `33.0` |
| `drop_factor_tree_depth` | `652` | `749` | `2322` | `3051` | `22.0` |
| `drop_multiply` | `707` | `749` | `2916` | `3051` | `41.0` |
| `drop_prime_lookup` | `722` | `749` | `2754` | `3051` | `43.0` |

## Base Tolerance Winners

| grid | dictionary | N | coeff bits | tolerance | winner | support | objective | residual |
|---|---|---:|---:|---:|---|---|---:|---:|
| legendre_2048 | geometric | 16 | 8 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 208.937 | 0.119689 |
| legendre_2048 | geometric | 16 | 8 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 157.162 | 0.139862 |
| legendre_2048 | geometric | 16 | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 157.162 | 0.139862 |
| legendre_2048 | geometric | 16 | 8 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 111.836 | 0.223159 |
| legendre_2048 | geometric | 16 | 16 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 304.797 | 0.108592 |
| legendre_2048 | geometric | 16 | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 221.081 | 0.132237 |
| legendre_2048 | geometric | 16 | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 221.081 | 0.132237 |
| legendre_2048 | geometric | 16 | 16 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 143.808 | 0.218866 |
| legendre_2048 | geometric | 16 | 24 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 400.796 | 0.108549 |
| legendre_2048 | geometric | 16 | 24 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 285.081 | 0.132208 |
| legendre_2048 | geometric | 16 | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 285.081 | 0.132208 |
| legendre_2048 | geometric | 16 | 24 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 175.808 | 0.218849 |
| legendre_2048 | geometric | 24 | 8 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 213.797 | 0.108579 |
| legendre_2048 | geometric | 24 | 8 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 162.237 | 0.147357 |
| legendre_2048 | geometric | 24 | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 162.237 | 0.147357 |
| legendre_2048 | geometric | 24 | 8 | 0.25 | composite_index | composite_prefix | 114.720 | 0.205964 |
| legendre_2048 | geometric | 24 | 16 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 309.629 | 0.0966558 |
| legendre_2048 | geometric | 24 | 16 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 309.629 | 0.0966558 |
| legendre_2048 | geometric | 24 | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 226.162 | 0.139832 |
| legendre_2048 | geometric | 24 | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 226.162 | 0.139832 |
| legendre_2048 | geometric | 24 | 16 | 0.25 | composite_index | composite_prefix | 146.691 | 0.201757 |
| legendre_2048 | geometric | 24 | 24 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 405.628 | 0.0966092 |
| legendre_2048 | geometric | 24 | 24 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 405.628 | 0.0966092 |
| legendre_2048 | geometric | 24 | 24 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 290.161 | 0.139803 |
| legendre_2048 | geometric | 24 | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 290.161 | 0.139803 |
| legendre_2048 | geometric | 24 | 24 | 0.25 | composite_index | composite_prefix | 178.691 | 0.201741 |
| legendre_2048 | geometric | 32 | 8 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 213.619 | 0.0960202 |
| legendre_2048 | geometric | 32 | 8 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 213.619 | 0.0960202 |
| legendre_2048 | geometric | 32 | 8 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 162.239 | 0.147477 |
| legendre_2048 | geometric | 32 | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 162.239 | 0.147477 |
| legendre_2048 | geometric | 32 | 8 | 0.25 | composite_index | composite_prefix | 117.718 | 0.205546 |
| legendre_2048 | geometric | 32 | 16 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 394.061 | 0.065202 |
| legendre_2048 | geometric | 32 | 16 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 309.422 | 0.0837122 |
| legendre_2048 | geometric | 32 | 16 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 309.422 | 0.0837122 |
| legendre_2048 | geometric | 32 | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 226.160 | 0.13967 |
| legendre_2048 | geometric | 32 | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 226.160 | 0.13967 |
| legendre_2048 | geometric | 32 | 16 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 146.989 | 0.248031 |
| legendre_2048 | geometric | 32 | 24 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 522.060 | 0.0651411 |
| legendre_2048 | geometric | 32 | 24 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 405.421 | 0.0836641 |
| legendre_2048 | geometric | 32 | 24 | 0.12 | heterogeneous_compute_reusable | integer_prefix | 405.421 | 0.0836641 |
| legendre_2048 | geometric | 32 | 24 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 290.160 | 0.139639 |
| legendre_2048 | geometric | 32 | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 290.160 | 0.139639 |
| legendre_2048 | geometric | 32 | 24 | 0.25 | heterogeneous_compute_reusable | factor_cost_ordered | 178.988 | 0.248014 |
| legendre_2048 | geometric | 48 | 8 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 368.281 | 0.0759563 |
| legendre_2048 | geometric | 48 | 8 | 0.10 | heterogeneous_compute_reusable | integer_prefix | 266.579 | 0.0933849 |
| legendre_2048 | geometric | 48 | 8 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 209.891 | 0.11591 |
| legendre_2048 | geometric | 48 | 8 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 162.002 | 0.125192 |
| legendre_2048 | geometric | 48 | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 162.002 | 0.125192 |
| legendre_2048 | geometric | 48 | 8 | 0.25 | composite_index | composite_prefix | 120.871 | 0.228665 |
| legendre_2048 | geometric | 48 | 16 | 0.06 | heterogeneous_compute_reusable | integer_prefix | 559.767 | 0.053166 |
| legendre_2048 | geometric | 48 | 16 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 394.300 | 0.0769641 |
| legendre_2048 | geometric | 48 | 16 | 0.10 | composite_index | composite_prefix | 327.589 | 0.094029 |
| legendre_2048 | geometric | 48 | 16 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 225.906 | 0.117104 |
| legendre_2048 | geometric | 48 | 16 | 0.25 | composite_index | composite_prefix | 152.844 | 0.224411 |
| legendre_2048 | geometric | 48 | 24 | 0.06 | heterogeneous_compute_reusable | integer_prefix | 751.764 | 0.053077 |
| legendre_2048 | geometric | 48 | 24 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 522.299 | 0.0769 |
| legendre_2048 | geometric | 48 | 24 | 0.10 | composite_index | composite_prefix | 423.589 | 0.0939825 |
| legendre_2048 | geometric | 48 | 24 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 289.905 | 0.117072 |
| legendre_2048 | geometric | 48 | 24 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 289.905 | 0.117072 |
| legendre_2048 | geometric | 48 | 24 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 289.905 | 0.117072 |
| legendre_2048 | geometric | 48 | 24 | 0.25 | composite_index | composite_prefix | 184.844 | 0.224395 |
| legendre_2048 | geometric | 64 | 8 | 0.08 | heterogeneous_compute_reusable | integer_prefix | 377.206 | 0.0721142 |
| legendre_2048 | geometric | 64 | 8 | 0.10 | composite_index | composite_prefix | 287.586 | 0.0938494 |
| legendre_2048 | geometric | 64 | 8 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 220.913 | 0.1177 |
| legendre_2048 | geometric | 64 | 8 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 220.913 | 0.1177 |
| legendre_2048 | geometric | 64 | 8 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 169.267 | 0.150378 |
| legendre_2048 | geometric | 64 | 8 | 0.25 | composite_index | composite_prefix | 120.926 | 0.237511 |
| legendre_2048 | geometric | 64 | 16 | 0.06 | heterogeneous_compute_reusable | integer_prefix | 568.625 | 0.0481861 |
| legendre_2048 | geometric | 64 | 16 | 0.08 | composite_index | composite_prefix | 415.321 | 0.0780948 |
| legendre_2048 | geometric | 64 | 16 | 0.10 | composite_index | composite_prefix | 327.608 | 0.0952694 |
| legendre_2048 | geometric | 64 | 16 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 316.766 | 0.106305 |
| legendre_2048 | geometric | 64 | 16 | 0.15 | heterogeneous_compute_reusable | factor_cost_ordered | 233.186 | 0.142167 |
| legendre_2048 | geometric | 64 | 16 | 0.20 | heterogeneous_compute_reusable | factor_cost_ordered | 233.186 | 0.142167 |
| legendre_2048 | geometric | 64 | 16 | 0.25 | composite_index | composite_prefix | 152.900 | 0.233301 |
| legendre_2048 | geometric | 64 | 24 | 0.06 | heterogeneous_compute_reusable | integer_prefix | 760.622 | 0.0480926 |
| legendre_2048 | geometric | 64 | 24 | 0.08 | composite_index | composite_prefix | 543.320 | 0.0780333 |
| legendre_2048 | geometric | 64 | 24 | 0.10 | composite_index | composite_prefix | 423.607 | 0.0952211 |
| legendre_2048 | geometric | 64 | 24 | 0.12 | heterogeneous_compute_reusable | factor_cost_ordered | 412.766 | 0.106261 |
| ... | ... | ... | ... | ... | 669 more rows | ... | ... | ... |

## Boundary

A replicated result here means only that the finite reusable
operation-weighted cost model survived this robustness gate against
the listed finite baselines and ablations. It does not establish
an infinite-dimensional theorem, dictionary invariance, or a
result about zeta zeros.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01_RESULTS.json`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01_REPORT.md`
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01_RESULTS.sha256`

No D1, scorecard, public page, git staging, commit, Zenodo, arXiv, or publisher surface was updated.
