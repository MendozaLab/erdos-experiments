# RH Beurling-Nyman MDL Quantization Stability Audit

Date: 2026-05-06

Status: local theorem/audit packet only. This is not RH evidence, not theorem progress, and not an asymptotic statement.

## The Modest Theorem

For a frozen finite Nyman-Beurling dictionary, weighted design matrix `A`,
target vector `y`, coefficient vector `a`, and componentwise rounded coefficient
vector `q_b(a)`,

```text
||A q_b(a)-y|| <= ||A a-y|| + sigma_max(A) sqrt(k) Delta_b / 2.
```

Here `k` is the active coefficient count and

```text
Delta_b = max(1, ||a||_infty) 2^(-b).
```

This is a stability theorem for a finite approximation artifact. It says when
reported MDL residuals survive a stated coefficient-bit grid and conditioning
penalty.

## Audit Requirement

Every finite row used in the MDL curve must expose:

- quantization step;
- quantization penalty;
- certified upper residual;
- condition number;
- claim ceiling.

Rows that do not expose those fields are not eligible for theorem-facing
interpretation.

## Interpretation Guard

The `~3.2` bit bend is named only as:

```text
first observed finite-N MDL conditioning crossover
```

It is not a universal constant. It is not evidence for the Riemann hypothesis.
It is an experimental place to test whether residual information starts paying
a visible Gram-conditioning tax after the first coarse digit of approximation.

## Operator Grammar

The operator overlay is useful as grammar, not proof:

```text
ratio       = invariant / normalization
square      = error, energy, dimensional ledger
square root = observable scale, amplitude, fluctuation scale
```

This makes Beurling-Nyman the legitimate Hilbert-space arena for the MDL idea,
while EHP114 remains a separate local proof lane.
