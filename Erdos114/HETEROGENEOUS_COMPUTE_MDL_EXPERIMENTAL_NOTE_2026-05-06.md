# Heterogeneous Compute MDL Experimental Note

Date: 2026-05-06
Status: internal experimental companion candidate
Primary experiment: `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01`
Replication experiment:
`EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01`

## Claim Ceiling

This note records finite operation-weighted Beurling-Nyman MDL diagnostics. It
is not a Riemann Hypothesis result, not a zeta formalization, and not a
quantum-mechanical result.

The publication split remains:

- `MS-103`: finite Lean stability theorem and certificate infrastructure;
- this note: experimental heterogeneous-compute diagnostic, internal unless
  replicated and separately approved.

## Why This Experiment Was Run

The earlier RH-MDL experiments gave useful failures:

```text
flat index bits
  -> too crude;

prime/composite projection ordering
  -> stable finite order effect, but not the expected composite-first story;

factorization harness bits
  -> no broad positive signal;

amortized factorization harness bits
  -> local wins, but no median positive signal.
```

The Transdimensional Painter precedent suggests the corrective move. In TDP,
raw-coordinate MDL failed, but mode-aware accounting recovered part of the
structure. The analogous RH-MDL question is whether finite Beurling-Nyman
supports are better described by arithmetic operations than by flat addresses
or symbol strings.

## Model

The objective is:

```text
total_objective_bits = compute_cost_bits - certified_residual_information_bits
```

The residual term uses the existing finite certificate logic:

```text
certified_residual_upper
  = fitted_residual + operator_norm * sqrt(active_count) * quantization_step / 2
```

The heterogeneous compute rows charge these finite channels:

- `prime_lookup`;
- `multiply`;
- `exponent`;
- `factor_tree_depth`;
- `coefficient_move`;
- optional `shared_harness`.

The decisive tested row family is `heterogeneous_compute_reusable`, which
charges selected arithmetic operations but treats the grammar as reusable
rather than paying a full setup penalty every row.

## Result

Run:

```text
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01
```

Verdict:

```text
STABLE_HETEROGENEOUS_COMPUTE_SIGNAL
```

Key finite facts:

- total rows: `1440`;
- tolerance rows: `118`;
- `heterogeneous_compute_reusable` won `100 / 118` tolerance rows;
- `composite_index` won the remaining `18 / 118` tolerance rows;
- the reusable heterogeneous row won `50` tolerance rows on `legendre_2048`
  and `50` on `midpoint_2048`;
- it won at every tested dictionary size: `N = 8, 16, 24, 32, 48`;
- same-support objective wins: `288 / 408`;
- median objective savings against the best listed baseline: `11.0` bits;
- adding the explicit shared-harness setup cost reversed the signal, with
  median savings `-40.0` bits.

## Interpretation

The useful reading is not "primes are residual columns" and not "factorization
strings are always shorter." Both earlier framings were too flat.

The useful reading is:

```text
integer supports look cheaper when arithmetic construction is charged as
reusable heterogeneous compute rather than as direct addressing.
```

That is a finite MDL statement about the tested design matrices, quadrature
rules, support schedules, and cost channels. It does not identify an
infinite-dimensional operator.

## Artifact Contract

Generated artifacts:

- `rh_bn_heterogeneous_compute_mdl.py`;
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_RESULTS.json`;
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_REPORT.md`;
- `EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_RESULTS.sha256`.

Verification:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114
python3 rh_bn_heterogeneous_compute_mdl.py
shasum -a 256 -c EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01_RESULTS.sha256
```

## Next Replication Gate

The first replication gate has now been run:

```text
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-REPLICATION-20260507-01
```

Verdict:

```text
REPLICATED_HETEROGENEOUS_COMPUTE_SIGNAL
```

Replication conditions:

- dictionary sizes: `N = 16, 24, 32, 48, 64, 96`;
- grids: `legendre_2048`, `midpoint_2048`, and `shifted_midpoint_2048`;
- coefficient precision settings: `8`, `16`, and `24` bits;
- ablations: remove one channel at a time from `prime_lookup`, `multiply`,
  `exponent`, and `factor_tree_depth`.

Key finite facts:

- total rows: `40014`;
- tolerance rows across all scenarios: `3745`;
- base reusable heterogeneous wins: `634 / 749` tolerance rows;
- base same-support wins: `2268 / 3051`;
- base median objective savings against the best listed baseline: `20.0`
  bits;
- base wins occurred on all three grids;
- base wins occurred at every tested dictionary size, including `64` and `96`;
- all four single-channel ablations passed the internal robustness rule.

Interpretation:

The first run was not a one-grid or small-`N` accident under this finite setup.
The reusable operation-weighted accounting signal survived the planned
replication gate. The ablation result does not mean that any one dropped
channel is irrelevant in an absolute sense; it means the finite advantage is
not carried entirely by one of the four tested arithmetic-operation charges.

## Next Robustness Gate

The second robustness gate has now been run:

```text
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01
```

Verdict:

```text
HARDENED_HETEROGENEOUS_COMPUTE_SIGNAL
```

Hardening conditions:

- support schedules were changed from prefix/factor-cost supports to
  regularized active-set supports;
- regularized supports included two ridge top-coefficient selectors and one
  ridge-refit OMP selector;
- dictionary sizes: `N = 16, 32, 48, 64, 96`;
- grids: `legendre_2048`, `midpoint_2048`, and `shifted_midpoint_2048`;
- coefficient precision settings: `8`, `16`, and `24` bits;
- operation-weight profiles: `base`, `all_unit`, `low_ops`, `high_ops`,
  `prime_heavy`, `operator_heavy`, and `free_arithmetic`;
- implementation audit independently recomputed operation-channel costs for
  every emitted row.

Key finite facts:

- total rows: `77112`;
- tolerance rows: `13853`;
- implementation audit: `PASS`, with `0` errors across `77112` checked rows;
- base reusable heterogeneous wins: `1546 / 1979` tolerance rows;
- base same-support wins: `2307 / 2754`;
- base median reusable objective savings against the best listed baseline:
  `25.0` bits;
- base wins occurred on all three grids;
- base wins occurred at every tested dictionary size, including `64` and `96`;
- base wins occurred under all three regularized support strategies;
- `4 / 7` weight profiles passed the internal robustness rule.

Important caveat:

The result is hardened, but not weight-invariant. The signal survives the base,
`all_unit`, `low_ops`, and `free_arithmetic` profiles. It fails under
`prime_heavy`, `high_ops`, and `operator_heavy` profiles. That means the finite
advantage is real under the tested reusable arithmetic accounting, but it is
not a free theorem about every possible operation-weight surface.

Interpretation:

The earlier positive result was not only a prefix-support artifact. It survived
regularized support search and an implementation audit. The next scientific
question is now calibration: what operation weights are principled, externally
defensible, and not tuned to the Beurling-Nyman residual table?

## Next Calibration Gate

Before this becomes an external companion paper, add:

- an externally motivated operation-weight schedule, such as arithmetic circuit
  cost, prefix-code decoding cost, or measured symbolic-computation cost;
- a negative-control dictionary where heterogeneous arithmetic accounting
  should not win;
- an independent reimplementation of the hardening script by another agent or
  in another language;
- a parameter-sweep report that makes the weight-sensitive failure modes as
  visible as the base-profile wins.

Until then, the result is a replicated internal finite diagnostic, not a public
mathematical claim.
