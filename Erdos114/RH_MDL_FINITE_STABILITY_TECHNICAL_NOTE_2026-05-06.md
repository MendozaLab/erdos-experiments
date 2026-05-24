# Finite MDL Stability for Hilbert-Space Beurling-Nyman Approximants

Internal draft v0.1
Date: 2026-05-06
Lean anchor: `ErdosLean4V427.RH.BeurlingNymanMDL`
Claim ceiling: finite Hilbert-space MDL stability for frozen finite
dictionaries; not a Riemann Hypothesis result

## Abstract

The Nyman-Beurling and Baez-Duarte criteria place the Riemann Hypothesis in a
Hilbert-space approximation setting. This note isolates a smaller unconditional
problem that appears in finite MDL experiments built from that setting:
coefficient quantization should not be counted as signal unless its induced
residual penalty is explicit.

For a bounded linear operator `A`, an unquantized coefficient vector `a`, a
perturbed or quantized vector `q`, and a target `y`, we prove

```text
||A q - y|| <= ||A a - y|| + ||A|| ||q - a||.
```

For `k` active Hilbert coordinates, if every coordinate error is bounded by
`Delta / 2`, then

```text
||q - a|| <= sqrt(k) * Delta / 2.
```

Combining the two gives the finite MDL certificate

```text
||A q - y||
  <= ||A a - y|| + ||A|| * sqrt(k) * Delta / 2.
```

The result is formalized in Lean 4. It is finite stability infrastructure for
description-length accounting; it does not prove a zero-free region, an
asymptotic approximation rate, or any statement about zeta zeros.

## 1. Context

The Nyman-Beurling criterion and Baez-Duarte's strengthening make the Riemann
Hypothesis equivalent to an `L^2` closure problem involving fractional-part
dilates. That literature already contains serious approximation, coefficient,
Gram, Hankel, and probabilistic-dictionary structure. The contribution here is
therefore intentionally narrower: define the finite quantization-stability
bound needed before a numerical MDL curve can be treated as an auditable finite
object.

In the finite experiments motivating this note, a fixed dictionary and a fixed
quadrature rule give a weighted least-squares problem

```text
A_ij = sqrt(w_i) rho_j(x_i)
y_i  = sqrt(w_i)
min_a ||A a - y||_2.
```

The MDL issue is not merely the fitted residual. The stored coefficient vector
has finite description length, and finite precision can be amplified by the
conditioned design operator. The certificate below makes that amplification
explicit.

## 2. Main Theorem

Let `A : E ->L[Real] F` be a bounded linear map between normed real spaces.
For all `a q : E` and `y : F`,

```text
||A q - y|| <= ||A a - y|| + ||A|| ||q - a||.
```

Now specialize the coefficient space to `EuclideanSpace Real (Fin k)`. If
`Delta >= 0` and

```text
for every coordinate i,  ||(q - a)_i|| <= Delta / 2,
```

then

```text
||q - a|| <= sqrt(k) * Delta / 2.
```

Consequently,

```text
||A q - y||
  <= ||A a - y|| + ||A|| * sqrt(k) * Delta / 2.
```

This is the finite residual penalty used by the MDL certificate. In experiment
language, `k` is the active coefficient count and `Delta` is the certified
coordinatewise quantization width.

## 3. Proof Sketch

The residual identity is

```text
A q - y = (A a - y) + A(q - a).
```

The triangle inequality gives

```text
||A q - y|| <= ||A a - y|| + ||A(q - a)||.
```

The operator norm gives

```text
||A(q - a)|| <= ||A|| ||q - a||.
```

For the coordinate estimate, the Hilbert norm on `EuclideanSpace Real (Fin k)`
is

```text
||x|| = sqrt(sum_i ||x_i||^2).
```

If each coordinate is at most `Delta / 2`, the sum of squared coordinate errors
is at most

```text
k * (Delta / 2)^2,
```

so taking square roots gives

```text
||x|| <= sqrt(k) * Delta / 2.
```

## 4. Lean Verification

The theorem is formalized in:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/Lean4/erdos-lean4-v427/ErdosLean4V427/RH/BeurlingNymanMDL.lean
```

Compiled theorem and interface names:

```lean
Erdos.RH_MDL.residual_bound
Erdos.RH_MDL.residual_bound_of_coeff_error
Erdos.RH_MDL.finite_coordinate_error_bound
Erdos.RH_MDL.CoordinatewiseQuantized
Erdos.RH_MDL.residual_bound_of_coordinate_error
Erdos.RH_MDL.residual_bound_of_coordinatewise_quantized
```

Verification commands:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/Lean4/erdos-lean4-v427
lake env lean ErdosLean4V427/RH/BeurlingNymanMDL.lean
lake build ErdosLean4V427.RH.BeurlingNymanMDL
rg -n "sorry|admit|axiom" ErdosLean4V427/RH/BeurlingNymanMDL.lean
```

The Lean formalization deliberately stops at a predicate named
`CoordinatewiseQuantized`. That predicate is the certified output of a
quantizer: each coordinate error is bounded by `Delta / 2`. A future file may
define a nearest-grid quantizer, but this note does not need that machinery.

## 5. Relationship To The Finite MDL Artifacts

The existing finite artifacts compute residuals, singular-value data,
condition numbers, coefficient bit-depths, and certified upper residuals:

- `EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-03_RESULTS.json`
- `EXP-MATH-RH-BEURLING-NYMAN-MDL-SENSITIVITY-20260506-01_RESULTS.json`
- `EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_RESULTS.json`

The compiled theorem justifies the deterministic penalty term

```text
operator_norm * sqrt(active_count) * quantization_step / 2.
```

If the quadrature weights are normalized so `||y|| = 1`, that penalty is on
the same scale as the certified relative residual reported by the finite MDL
artifacts.

## 6. Prior Art Positioning

This note does not claim that Nyman-Beurling is newly an approximation problem.
That is the original point of the criterion and its later refinements.

Required citation line:

- Nyman (1950) and Beurling (1955) introduce the Hilbert-space approximation
  criterion.
- Baez-Duarte (2002/2003, 2005) strengthens and generalizes the criterion.
- Burnol (2001/2002) gives quantitative lower-bound structure involving zeta
  zeros.
- Alouges-Darses-Hillion (2022) and Darses-Hillion (2021) cover polynomial,
  probabilistic, coefficient-control, and Gram/Hankel variants.
- Lagarias (2004/2007) is relevant for Li-coefficient and positivity
  context, but this note does not identify Li coefficients with MDL slack.
- Recent gray literature on Preprints.org describes dictionary-approximation
  and compressed-sensing language for the Baez-Duarte problem. That preempts
  broad novelty claims around "dictionary approximation" itself.

The remaining defensible novelty is only the finite MDL accounting layer:
given a frozen finite dictionary and cost model, coefficient quantization has a
compiled residual penalty.

## 7. Claim Filter

Safe claims:

- A finite Hilbert-space residual perturbation theorem has been formalized in
  Lean.
- A finite coordinatewise coefficient-error bound has been formalized in Lean.
- Together, these prove the deterministic quantization penalty used in the
  finite MDL certificate for a frozen dictionary.

Unsafe claims:

- The theorem proves, supports, or gives evidence for the Riemann Hypothesis.
- The theorem proves any asymptotic Nyman-Beurling approximation rate.
- The theorem proves invariance of the MDL curve across dictionary choices.
- The theorem formalizes zeta, prime distribution, fractional-part functions,
  or closure of the Nyman-Beurling span.

## 8. Bibliography

1. B. Nyman, "On the One-Dimensional Translation Group and Semi-Group in
   Certain Function Spaces," thesis, Uppsala, 1950.
2. A. Beurling, "A closure problem related to the Riemann zeta-function,"
   Proceedings of the National Academy of Sciences 41 (1955), 312-314.
   DOI: `10.1073/pnas.41.5.312`.
3. L. Baez-Duarte, "A strengthening of the Nyman-Beurling criterion for the
   Riemann Hypothesis," arXiv: `math/0202141`.
4. L. Baez-Duarte, "A general strong Nyman-Beurling Criterion for the Riemann
   Hypothesis," arXiv: `math/0505453`.
5. J.-F. Burnol, "A lower bound in an approximation problem involving the zeros
   of the Riemann zeta function," arXiv: `math/0103058`.
6. F. Alouges, S. Darses, and M. Hillion, "Polynomial approximations in a
   generalized Nyman-Beurling criterion," Journal de Theorie des Nombres de
   Bordeaux 34 (2022), 767-785.
7. S. Darses and M. Hillion, "On probabilistic generalizations of the
   Nyman-Beurling criterion for the zeta function," Confluentes Mathematici 13
   (2021), 43-59.
8. J. C. Lagarias, "Li Coefficients for Automorphic L-Functions," arXiv:
   `math/0404394`.
9. "Spectral and Analytic Structure of the Nyman-Beurling-Baez-Duarte
   Approximation," Preprints.org manuscript `202506.0772`; gray literature,
   cited only as dictionary/compressed-sensing framing risk.

## 9. Next Work

1. Define a concrete nearest-grid quantizer and prove it satisfies
   `CoordinatewiseQuantized`.
2. Specialize the abstract finite theorem to the stored Beurling-Nyman design
   matrix schema.
3. Run a composite-first residual experiment that separates composite bulk
   compression from prime-indexed marginal gain.
4. Only after those checks, prepare an external version with DOI/arXiv-style
   formatting and a fresh prior-art pass.
