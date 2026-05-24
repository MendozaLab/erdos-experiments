# EXP-MATH-EHP114-HESSIAN-FANTASY-20260505-01 Report

## Status

Diagnostic aggregate, not a proof of EHP #114 and not a publication packet.

The ordinary smooth-Hessian fantasy is dead as a radial local model: the radial
boundary layer follows a Puiseux exponent close to 1/n, not an eps^2 law.

The useful survivor is stratified. Read the certificate as three pieces:
outer-domain coercivity, radial Puiseux lower bound, and nonradial shape
Hessian/remainder control. This is a shadow signature, not universal law.

## Verdict

- Overall decision: `STRATIFIED`
- Smooth Hessian layer: `DEAD_AS_ORDINARY_RADIAL_MODEL`
- Stratified certificate layer: `ALIVE_AS_DIAGNOSTIC_ANSATZ`

Meaning: do not try to make a single ordinary Hessian at `z^n - 1` carry the
boundary. The radial axis is singular. The finite-difference shape probes still
show positive Hessian-like curvature on quotient Fourier modes, so the right
attack is a decomposed certificate rather than a smooth Taylor certificate.

## Radial Puiseux Evidence

| n | expected 1/n slope | tail-window 4 slope | abs error | relative error | deficit / eps^(1/n) tail proxy |
|---:|---:|---:|---:|---:|---:|
| 14 | 0.0714285714286 | 0.0714266524662 | 1.91896239662e-06 | 2.68654735527e-05 | 29.2164468516 |
| 15 | 0.0666666666667 | 0.0666648756369 | 1.79102980891e-06 | 2.68654471336e-05 | 31.2911143052 |

The n = 14 and n = 15 radial probes both land on the predicted 1/n exponent to
small relative error in the tail fit. This kills the radial quadratic model and
supports a fixed-n Puiseux lower-bound theorem as the radial component.

## Tensor-Cone Shape Probe

- Source: `EXP-MATH-EHP114-N10-TENSOR-CONE-SCAFFOLD-20260502-02`
- Degree: 10
- Rigorous: False
- Singular radial slope: 0.099828404964
- Shape slope range: 0.110969377353 to 0.200796273453
- Shape slope mean: 0.157571184867
- Shape slopes >= 1.5: 0
- Shape slopes < 1.5: 16
- All tested shape symmetric deficits positive: True
- Mixed shape proxy eigenvalues: min 160875.053745, max 498598.792372, positive 16/16
- Condition proxy: 3.09929215728

Interpretation: the tensor-cone data does not rescue a smooth eps^2 scale law,
because every measured shape slope is sub-Hessian. It does preserve a positive
finite-difference shape cone at the tested scale. The shape layer is therefore
Hessian-like as a cone positivity probe, not as a full smooth local model.

## Fourier-Hessian Shape Probe

- Source: `EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02`
- Degree: 15
- Rigorous: False
- Basis rank: 27 = expected reduced dimension 27
- Positive modes at all eps: 27/27
- Rows not positive at all eps: 0
- Curvature range: 68179.5415698 to 1569895.84277
- Large reference-error flag: True

Interpretation: all 27 n = 15 Fourier modes carry positive deficit curvature
under this floating estimator, but the singularity sanity flag matters. This is
evidence for the nonradial cone component, not an interval certificate.

## Certificate Ansatz

Let `D_n(p) = L(z^n - 1) - L(p)` and decompose a local perturbation as
`u = r e_0 + s`, where `e_0` is the singular radial direction and `s` is the
quotient nonradial Fourier/tensor shape component.

The surviving stratified certificate should have this form:

```text
Outer domain:
  D_n(p) >= B_n(rho) > 0
  for p outside the singular local cone / middle annulus.

Radial axis:
  D_n(r e_0) >= C_n |r|^(1/n)
  with fixed-n interval constants, starting at n = 14.

Shape cone:
  D_n(r e_0 + s) >= C_n |r|^(1/n) + lambda_n ||s||^2 - R_n(r, s)
  where lambda_n > 0 comes from the nonradial cone and R_n is bounded small
  enough that mixed radial/shape terms cannot cancel the Puiseux deficit.
```

That is the honest target. The current evidence supports the architecture but
does not establish the analytic remainder bound.

## Next Run Target

Recommended next run: n = 20 triad, because it is the cleanest falsification
test for the stratified fantasy.

Inputs:

- Exact radial hypergeometric family for n = 20 with eps grid
  `(0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001, 1e-05, 1e-06, 1e-07, 1e-08)`.
- Fourier-Hessian basis rank `2n - 3 = 37`, eps `(0.02, 0.01, 0.005)`, grid resolution at least matching n = 15.
- Tensor-cone mixed matrix for the n = 20 quotient shape basis, with the singular radial direction excluded from the mixed shape block.

Outputs:

- `EXP-MATH-EHP114-N20-RADIAL-HYPERGEOMETRIC-20260505-02_RESULTS.json`
- `EXP-MATH-EHP114-N20-FOURIER-HESSIAN-20260505-02_RESULTS.json`
- `EXP-MATH-EHP114-N20-TENSOR-CONE-SCAFFOLD-20260505-02_RESULTS.json`

Pass signal: radial slope within 5 percent of 1/20, all or nearly all Fourier
shape modes positive, and a positive mixed shape proxy with a bounded condition
proxy. Fail signal: radial slope drifts from 1/20 or the shape cone loses broad
positive curvature.

Proof-oriented alternative: n = 14 interval hardening for the fixed
`gauss_2F1_puiseux_lower_bound_n14_interval` theorem, then a separate n = 14 or
n = 15 Fourier/tensor interval cone for nonradial and mixed terms.

## Source Boundary

All numbers in this aggregate are computed from saved local artifacts listed in
the companion results JSON. No D1, scorecard, public document, email, push, or
publication surface was changed.
