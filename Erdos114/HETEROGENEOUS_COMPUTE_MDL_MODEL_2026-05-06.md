# Heterogeneous Compute MDL Model

Date: 2026-05-06
Status: internal research model

## Claim Ceiling

This note defines a finite operation-weighted MDL diagnostic for
Beurling-Nyman approximants. It is not a Riemann Hypothesis claim, not a zeta
formalization, and not a quantum-mechanical claim.

## Why This Model Exists

The previous finite diagnostics narrowed the problem:

```text
composite-first projection
  -> stable finite asymmetry, but mostly prime_first;

prime-harness factorization bits
  -> no broad positive signal;

amortized prime-harness bits
  -> local wins, but no median positive signal.
```

Those failures used mostly homogeneous bit accounting: each index or factor
description was charged as a symbol string. The Transdimensional Painter
precedent suggests a better test. In the N-body/Koopman work, raw-coordinate MDL
failed, but mode-aware MDL recovered some structure. For RH-MDL, the analogous
move is to test arithmetic compute channels instead of flat index strings.

## Revised Hypothesis

```text
flat address bits failed;
factorization bits failed;
test arithmetic operations as heterogeneous compute channels.
```

The finite question is:

```text
Do Beurling-Nyman supports look cheaper when dictionary indices are charged as
arithmetic compute graphs rather than direct addresses?
```

## Cost Channels

The v1 heterogeneous cost surface reports:

- `direct_index`: direct finite lookup of selected dictionary columns;
- `prime_lookup`: lookup of distinct prime tokens used by selected indices;
- `multiply`: factor-composition operations needed to build composite indices;
- `exponent`: repeated-factor / prime-power operations;
- `factor_tree_depth`: depth penalty for composing factors;
- `coefficient_move`: finite coefficient movement at the same 16-bit precision
  used by the residual certificate;
- `shared_harness`: optional setup cost for making the prime harness available.

These are bit-equivalent compute units, not hardware measurements. The point is
to test whether a structured arithmetic operation surface is more explanatory
than flat address length.

## Objective

Each row stores:

```text
certified_residual_information_bits = -log2(certified_upper_relative_residual)
compute_cost_bits
total_objective_bits = compute_cost_bits - certified_residual_information_bits
```

Lower `total_objective_bits` is better. This is equivalent to adding a residual
log-cost to the compute cost: better residuals reduce the objective, while more
compute increases it.

The experiment also keeps the stricter tolerance tables from the earlier
prime-harness diagnostics, so a low objective cannot hide a useless residual.

## Baselines

The experiment compares the new heterogeneous rows against the existing finite
encodings:

- `flat_index`;
- `composite_index`;
- `factorized_reusable_harness`;
- `heterogeneous_compute_reusable`;
- `heterogeneous_compute_with_harness`.

## Success Criterion

A positive v1 signal requires stable finite behavior across:

- both quadrature grids;
- at least three dictionary sizes;
- the existing baseline encodings.

The result is still only a finite diagnostic. It does not imply an asymptotic
law or dictionary-invariant MDL curve.

## Safe / Unsafe Language

Safe:

- heterogeneous compute MDL diagnostic;
- arithmetic operation cost surface;
- finite Beurling-Nyman support accounting;
- TDP-inspired mode-aware cost model.

Unsafe:

- claiming RH support;
- claiming zeta is an operator produced by the cost model;
- claiming a quantum result;
- treating a finite positive result as an infinite-dimensional theorem.

## Run 01 Result

Experiment:

```text
EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-20260506-01
```

Status:

```text
STABLE_HETEROGENEOUS_COMPUTE_SIGNAL
```

Summary:

- total rows: `1440`;
- tolerance rows: `118`;
- tolerance winners: `heterogeneous_compute_reusable` won `100`, while
  `composite_index` won `18`;
- grid stability: `heterogeneous_compute_reusable` won `50` tolerance rows on
  `legendre_2048` and `50` on `midpoint_2048`;
- dictionary-size stability: `heterogeneous_compute_reusable` won at every
  tested size `N = 8, 16, 24, 32, 48`;
- same-support comparable rows: `heterogeneous_compute_reusable` improved the
  objective in `288 / 408` cases;
- median reusable objective savings against the best listed baseline: `11.0`
  bits;
- the explicit shared-harness variant was negative, with median savings
  `-40.0` bits.

Interpretation:

The finite signal is not that every arithmetic setup cost should be charged
up front. The signal is narrower and more useful: once the arithmetic grammar
is treated as reusable, operation-weighted compute accounting beats the flat
and factorized baselines under the stated finite objective. That is exactly the
TDP lesson translated into this setting: raw coordinate strings can miss a
mode-aware representation.

Boundary:

This result is an internal finite diagnostic. It does not establish an
infinite-dimensional theorem, a dictionary-invariant MDL curve, or a result
about zeta zeros.
