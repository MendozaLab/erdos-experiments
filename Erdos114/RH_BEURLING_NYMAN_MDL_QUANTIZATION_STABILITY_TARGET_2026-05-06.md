# RH Beurling-Nyman MDL Quantization Stability Target

Date: 2026-05-06
Depends on: `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-02`
Status: THEOREM TARGET, not RH evidence

## Meaning

The overlook is now past prose. We have a frozen finite approximation object:
Nyman-Beurling basis functions, fixed dictionaries, a bit-cost model, numerical
least-squares residuals, coefficient quantization, and stored singular-value
bounds. The next useful theorem is not about RH. It is the stability theorem
that says when a finite MDL curve is mathematically auditable rather than a
floating-point artifact.

## Frozen Finite Object

For a fixed dictionary `Theta_N = {theta_1, ..., theta_N}`, define

```text
rho_j(x) = {theta_j / x}
```

on `(0,1]`, where `{.}` is fractional part. Let the numerical quadrature nodes
and weights be `(x_i, w_i)`. The weighted design matrix and target are

```text
A_ij = sqrt(w_i) rho_j(x_i)
y_i  = sqrt(w_i).
```

The finite least-squares problem is

```text
min_a ||A a - y||_2.
```

The stored experiment reports `sigma_max(A)`, `sigma_min(A)`, Gram eigenvalue
bounds, unquantized residual, and quantized residuals for multiple coefficient
bit-depths.

## Theorem Target

Let `a` be any coefficient vector and let `q_b(a)` be componentwise rounding to
a grid of width

```text
Delta_b = max(1, ||a||_infty) 2^(-b).
```

Then

```text
||A q_b(a) - y||_2
  <= ||A a - y||_2 + sigma_max(A) ||q_b(a) - a||_2
  <= ||A a - y||_2 + sigma_max(A) sqrt(k) Delta_b / 2,
```

where `k` is the number of active coefficients.

If `||y||_2 = 1` under normalized quadrature weights, the same inequality is a
relative residual bound.

## Proof Sketch

Use the triangle inequality:

```text
A q_b(a) - y = (A a - y) + A(q_b(a) - a).
```

Then

```text
||A(q_b(a) - a)||_2 <= ||A||_2 ||q_b(a) - a||_2.
```

By definition, `||A||_2 = sigma_max(A)`. Componentwise nearest-grid rounding
gives `|q_b(a_j) - a_j| <= Delta_b / 2` on each active coordinate, hence

```text
||q_b(a) - a||_2 <= sqrt(k) Delta_b / 2.
```

That proves the target inequality.

## Concrete Check From The `-02` Run

For the best `N=32` geometric dictionary row:

```text
relative residual        = 0.06143014218981282
condition number         = 44.66946219129291
sigma_max(A)             = 2.2266407891972744
sigma_min(A)             = 0.04984704717647792
gram eigen max           = 4.957929204117061
gram eigen min           = 0.0024847281122140153
active coefficient count = 32
```

The 8-bit, 16-bit, and 24-bit quantized residuals are essentially unchanged in
the smoke run. The theorem above explains what must be checked before treating
that as signal: the quantization error term must be small relative to the
unquantized residual after multiplying by `sigma_max(A)`.

## What This Theorem Would Support

Safe:

- "The finite MDL curve is stable under stated coefficient quantization for a
  frozen dictionary and quadrature rule."
- "Residual improvements are reported together with a deterministic
  quantization penalty."

Unsafe:

- "This supports RH."
- "This defines a dictionary-invariant MDL curve."
- "This proves a Nyman-Beurling asymptotic rate."

## Next Implementation Target

Add the quantization penalty explicitly to the results schema:

```text
quantization_penalty_bound =
  sigma_max * sqrt(active_count) * quantization_step / 2

certified_upper_relative_residual =
  relative_residual_unquantized + quantization_penalty_bound / ||y||_2
```

Then rerun as a new immutable experiment ID. Do not overwrite `-02`.
