# EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-02 Report

## Claim Tested

Can the EHP114 -> RH overlook be moved past prose into a frozen numerical
object: a Nyman-Beurling description-length curve that reports residual decay,
coefficient quantization, and Gram conditioning?

## Claim Ceiling

SUGGESTIVE / METHOD-SHAPING ONLY: finite quadrature probe of a frozen Nyman-Beurling MDL cost model; not RH evidence and not theorem progress.

This version supersedes `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-01`
because the first smoke run did not store the singular-value bounds needed for
the quantization-stability theorem target.

This report does not support any public RH claim. It does not update D1,
scorecards, Lean status, or publication status.

## Frozen Model

- Target: `chi_(0,1] = 1` on `(0,1]`.
- Basis: `rho_theta(x) = fractional_part(theta / x)`.
- Norm: numerical `L2(0,1)` using 2048 Gauss-Legendre nodes.
- Dictionaries: harmonic, geometric, and seeded log-uniform.
- Bit cost: dictionary header + support-pattern bits + coefficient bits.
- Quantization: coefficient grid step `max(1, ||a||_inf) * 2^-b`.

## Result

The probe runs and produces a real MDL-style curve, but the first result is a
warning rather than a breakthrough. At `N=32`, the best unquantized
dictionary is `geometric` with relative residual
`0.0614301` and condition number
`44.66946219129291`. With 16-bit coefficient quantization, the
best `N=32` curve has relative residual
`0.0614301` at
`577` description bits.

Interpretation: the curve is computable, but any later theorem must control
conditioning. Residual improvements are not meaningful unless they survive
quantization and a stated Gram-condition penalty.

## Best Unquantized Fits

| N | dictionary | relative residual | -log2 residual | condition number | active |
|---:|---|---:|---:|---:|---:|
| 4 | geometric | 0.196236 | 2.3493 | 4.985478869606194 | 4 |
| 8 | geometric | 0.141771 | 2.8184 | 9.356113269560419 | 8 |
| 12 | geometric | 0.11774 | 3.0863 | 14.406601292592795 | 12 |
| 16 | seeded_log_uniform | 0.0928897 | 3.4283 | 41.5481223218339 | 16 |
| 24 | geometric | 0.0919214 | 3.4435 | 36.393461195238956 | 24 |
| 32 | geometric | 0.0614301 | 4.0249 | 44.66946219129291 | 32 |

## Best Quantized Fits

| N | coeff bits | dictionary | relative residual | -log2 residual | description bits | condition number |
|---:|---:|---|---:|---:|---:|---:|
| 4 | 8 | geometric | 0.196243 | 2.3493 | 97 | 4.985478869606194 |
| 4 | 16 | geometric | 0.196236 | 2.3493 | 129 | 4.985478869606194 |
| 4 | 24 | geometric | 0.196236 | 2.3493 | 161 | 4.985478869606194 |
| 8 | 8 | geometric | 0.141773 | 2.8183 | 129 | 9.356113269560419 |
| 8 | 16 | geometric | 0.141771 | 2.8184 | 193 | 9.356113269560419 |
| 8 | 24 | geometric | 0.141771 | 2.8184 | 257 | 9.356113269560419 |
| 12 | 8 | geometric | 0.117743 | 3.0863 | 161 | 14.406601292592795 |
| 12 | 16 | geometric | 0.11774 | 3.0863 | 257 | 14.406601292592795 |
| 12 | 24 | geometric | 0.11774 | 3.0863 | 353 | 14.406601292592795 |
| 16 | 8 | seeded_log_uniform | 0.0929266 | 3.4278 | 225 | 41.5481223218339 |
| 16 | 16 | seeded_log_uniform | 0.0928897 | 3.4283 | 353 | 41.5481223218339 |
| 16 | 24 | seeded_log_uniform | 0.0928897 | 3.4283 | 481 | 41.5481223218339 |
| 24 | 8 | geometric | 0.0919257 | 3.4434 | 257 | 36.393461195238956 |
| 24 | 16 | geometric | 0.0919214 | 3.4435 | 449 | 36.393461195238956 |
| 24 | 24 | geometric | 0.0919214 | 3.4435 | 641 | 36.393461195238956 |
| 32 | 8 | geometric | 0.0614505 | 4.0244 | 321 | 44.66946219129291 |
| 32 | 16 | geometric | 0.0614301 | 4.0249 | 577 | 44.66946219129291 |
| 32 | 24 | geometric | 0.0614301 | 4.0249 | 833 | 44.66946219129291 |

## What This Supports

- A concrete `K_N(epsilon)` object now exists for the frozen dictionaries and
  bit-cost model.
- The experiment surfaces the right artifact risks: support size, coefficient
  precision, and condition number.
- This is enough to design a first theorem about quantization stability.

## What This Does Not Support

- It does not support RH.
- It does not show a new RH-equivalent criterion.
- It does not show the MDL curve is dictionary-invariant.
- It does not show the residual decay is asymptotic rather than finite-N
  numerical behavior.

## Next Theorem Target

For a fixed dictionary matrix `A`, target `y`, least-squares coefficient vector
`a`, and quantized vector `q_b(a)`, prove a deterministic bound of the form:

```text
||A q_b(a) - y||_2
  <= ||A a - y||_2 + ||A||_2 * ||q_b(a) - a||_2
```

Then express `||A||_2` and the allowed quantization step in terms of the stored
Gram singular values and coefficient norm. This theorem would not mention RH;
it would make the numerical MDL curve auditable.

## Artifact Boundary

Generated artifacts:

- `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-02_RESULTS.json`
- `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-02_REPORT.md`
- `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-02_RESULTS.sha256`

No D1, scorecard, public page, Lean status, git staging, commit, or publication
surface was updated.
