# Composite-First Projection Model

Date: 2026-05-06
Status: internal research model

## Claim Ceiling

This note defines a finite projection diagnostic for Beurling-Nyman MDL
experiments. It is not a Riemann Hypothesis claim, not a zeta formalization, and
not a public spectral-zeta story.

## Meaning

The safe version of "composite numbers fall by the wayside" is a finite linear
algebra question:

```text
Do composite-indexed columns absorb a large low-cost bulk component first,
while prime-indexed columns supply harder marginal residual directions?
```

The finite Hilbert-space arena is the weighted Beurling-Nyman design matrix

```text
A_ij = sqrt(w_i) * rho_j(x_i),
rho_j(x) = {theta_j / x},
y_i = sqrt(w_i).
```

The diagnostic splits dictionary columns by the natural index `j`:

- `j = 1` is the unit column;
- prime `j` are prime-indexed columns;
- composite `j > 1` are composite-indexed columns.

The unit column is included in the first-stage standalone fits so that prime
and composite comparisons are not dominated by the special role of `1`.

## Projection Tests

For each frozen dictionary, grid, and dictionary size `N`, compute:

1. unit-plus-composite fit;
2. unit-plus-prime fit;
3. all-column fit;
4. composite-first sequential residual fit;
5. prime-first sequential residual fit.

The sequential fits are intentionally order-sensitive:

```text
composite first:
  y -> project onto span(unit, composites)
  residual -> project onto span(primes)

prime first:
  y -> project onto span(unit, primes)
  residual -> project onto span(composites)
```

This is not a theorem about primes. It is a finite diagnostic for whether the
dictionary geometry contains a stable projection asymmetry.

## Metrics

Each row reports:

- residual norm and relative residual;
- residual drop per added column;
- residual drop per 16-bit coefficient bit;
- singular values and condition numbers;
- Gram eigenvalue bounds;
- quantization penalty bound from the compiled Lean stability theorem;
- finite projection spectrum diagnostics for the second-stage residual.

The spectral diagnostic is finite-dimensional only. It asks whether the energy
captured by the second-stage projection is concentrated in a small number of
singular-vector coordinates.

## Success Criterion

The composite-first story survives only if the same order asymmetry is stable
across:

- at least three dictionary sizes;
- at least two quadrature grids;
- both composite-first and prime-first orderings.

If the sign flips across grids or dictionary families, the result remains a
useful conditioning diagnostic but not a stable composite-first phenomenon.

## Run 01 Result

Immutable experiment:

```text
EXP-MATH-RH-BN-COMPOSITE-FIRST-PROJECTION-20260506-01
```

The first run found a stable finite projection asymmetry, but the dominant
order was `prime_first`, not `composite_first`.

```text
prime_first rows:      19 / 24
composite_first rows:   5 / 24
tie rows:               0 / 24
```

The reversal matters. The naive statement "composite columns absorb the easy
bulk before primes matter" is not the result of this first finite test. The
more defensible research statement is:

```text
the frozen Beurling-Nyman dictionaries show an order-sensitive projection
asymmetry, and the hard split is measurable by finite residual accounting.
```

At `N = 32` some rows flip toward `composite_first`, so the next experiment
should test larger `N`, more grids, and active-set regularization before any
interpretive language is strengthened.

## Safe Language

Safe:

- finite projection spectrum;
- composite-first residual diagnostic;
- prime-indexed marginal gain;
- finite Beurling-Nyman MDL projection experiment.

Unsafe:

- identifying zeta with the finite projector;
- describing the experiment as a proof mechanism;
- claiming any RH support;
- claiming an asymptotic law.
