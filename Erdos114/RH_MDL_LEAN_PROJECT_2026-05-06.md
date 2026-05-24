# RH-MDL Lean Project

Date: 2026-05-06
Status: Lean-first internal project, finite stability theorem compiled
Home lane: `erdos-experiments/Erdos114/`
Lean module: `Lean4/erdos-lean4-v427/ErdosLean4V427/RH/BeurlingNymanMDL.lean`
Technical note: `RH_MDL_FINITE_STABILITY_TECHNICAL_NOTE_2026-05-06.md`

## Purpose

This project turns the "RH in Hilbert space" synthesis into a small compiled
Lean target. The load-bearing mathematical object is not RH itself. It is the
finite Hilbert-space stability theorem behind the Beurling-Nyman MDL
quantization packet:

```text
||A q - y|| <= ||A a - y|| + ||A|| ||q - a||.
```

In the finite MDL probe, `A` is the weighted design operator, `a` is the
unquantized coefficient vector, `q` is a perturbed or quantized coefficient
vector, and `y` is the weighted target. The lemma says that a coefficient
quantization error becomes a residual penalty only through the operator norm of
the design map.

## Source Artifacts

- Downloads synthesis note:
  `/Users/kenbengoetxea/Downloads/of rh is in hilbert space how to start___The clean.md`
- Theorem target:
  `RH_BEURLING_NYMAN_MDL_QUANTIZATION_STABILITY_TARGET_2026-05-06.md`
- Main finite probe:
  `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-03_RESULTS.json`
- Sensitivity run:
  `EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01_RESULTS.json`
- Finite audit:
  `EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_RESULTS.json`
- Prior-art memo:
  `EHP114_RH_MDL_PRIOR_ART_PERPLEXITY_2026-05-06.md`

## Claim Ceiling

Safe claims:

- The v1 Lean file proves a finite Hilbert-space residual perturbation bound.
- The theorem supplies the formal algebraic core of the current MDL
  quantization-stability packet.
- The result is compatible with the Beurling-Nyman finite least-squares
  experiment because finite weighted design matrices are bounded linear maps.

Unsafe claims:

- This bears on RH itself, either as proof, support, or numerical indication.
- This proves an asymptotic Beurling-Nyman rate.
- This defines an invariant MDL curve across dictionaries.
- This models zeta, zeros, prime distribution, fractional-part functions, or
  componentwise rounding in Lean.

## Compiled Theorem Interfaces

Namespace: `Erdos.RH_MDL`

```lean
residual_bound
```

For normed real spaces and a bounded linear map `A : E ->L[Real] F`, prove:

```text
||A q - y|| <= ||A a - y|| + ||A|| * ||q - a||
```

```lean
residual_bound_of_coeff_error
```

Assuming an external coefficient error estimate `||q - a|| <= delta`, prove:

```text
||A q - y|| <= ||A a - y|| + ||A|| * delta
```

```lean
finite_coordinate_error_bound
```

For `x : EuclideanSpace Real (Fin k)`, assuming `Delta >= 0` and each
coordinate is bounded by `Delta / 2`, prove:

```text
||x|| <= sqrt(k) * Delta / 2
```

```lean
CoordinatewiseQuantized
```

A paper-facing predicate recording the certified output of a quantizer:

```text
for every i, ||(q - a)_i|| <= Delta / 2
```

```lean
residual_bound_of_coordinate_error
```

Combining the two pieces, assuming each coordinate of `q - a` is bounded by
`Delta / 2`, prove:

```text
||A q - y|| <= ||A a - y|| + ||A|| * sqrt(k) * Delta / 2
```

```lean
residual_bound_of_coordinatewise_quantized
```

The same residual certificate, using `CoordinatewiseQuantized` as the input
interface instead of a raw coordinate hypothesis.

## Build Instructions

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/Lean4/erdos-lean4-v427
lake env lean ErdosLean4V427/RH/BeurlingNymanMDL.lean
lake build ErdosLean4V427.RH.BeurlingNymanMDL
rg -n "sorry|admit|axiom" ErdosLean4V427/RH/BeurlingNymanMDL.lean
```

Build PASS plus an empty integrity scan means v1 has achieved its local goal.
It does not update D1, the master scorecard, Zenodo, arXiv, or any public
surface.

## Next Targets

1. Nearest-grid rounding theorem:

```text
nearest_grid_quantize a Delta satisfies CoordinatewiseQuantized
```

2. Componentwise rounding/grid model:

```text
Delta_b = max(1, ||a||_infty) * 2^(-b)
```

3. Beurling-Nyman specialization:

```text
A_ij = sqrt(w_i) * rho_j(x_i)
y_i  = sqrt(w_i)
```

4. Composite-first prime-residual experiment:

Restrict the finite dictionary to composite-indexed basis functions, measure
the residual, then add prime-indexed basis functions back by marginal
description-length gain.

## Operator-Grammar Reading

The project keeps the interpretive grammar as commentary, not proof data:

- ratios encode invariant comparison;
- squares encode Hilbert-space error and dimension;
- square roots encode observable fluctuation scale.

The Lean theorem formalizes only the square-error stability step and the finite
coordinate-counting penalty. That is the honest bridge from synthesis prose to
a compiled mathematical artifact.
