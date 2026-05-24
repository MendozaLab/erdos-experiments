# EXP-MATH-EHP114-N10-TENSOR-CONE-SCAFFOLD-20260502-02 Report

## Status

Diagnostic tensor-cone scaffold. Not a proof.

## Claim Tested

Does the `n = 10` quotient tangent cone look like an ordinary smooth Hessian
problem, or does it already show stratified/Puiseux behavior?

## Reference

- Degree: 10
- Exact L*: `22.886060328165432870217912684368174106539132975253`
- Reference artifact: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n10-inari_RESULTS.json`
- Artifact verdict: `EHP_N10_PROVEN`
- Exact value inside artifact interval: True

## Basis

- Quotient basis rank: 17
- Expected rank `2n - 3`: 17
- Shape basis rank after removing singular radial: 16
- eps values: [0.04, 0.02, 0.01, 0.005]
- Grid resolution: 340

## Summary

- Radial slope: 0.099828
- Shape slope min: 0.110969
- Shape slope mean: 0.157571
- Shape slope max: 0.200796
- Shape directions with slope >= 1.5: 0
- Shape directions with slope < 1.5: 16
- All directions have positive symmetric deficit: True

## Mixed Shape Hessian Proxy

- eps: 0.005
- min eigenvalue: 1.608751e+05
- max eigenvalue: 4.985988e+05
- condition proxy abs(max)/abs(min): 3.099292e+00

This matrix is a finite-difference proxy, not an interval Hessian.

## Direction Table

| label | role | mode | kind | phase | fitted slope | smallest-eps D/eps^2 | positive? |
|---|---|---:|---|---|---:|---:|---|
| radial_singular_m0 | singular_radial | 0 | radial | constant | 0.099828 | 4.912128e+05 | yes |
| m1_cos_radial | shape_candidate | 1 | radial | cos | 0.111287 | 4.024963e+05 | yes |
| m1_cos_tangent | shape_candidate | 1 | tangent | cos | 0.110969 | 4.024761e+05 | yes |
| m2_cos_radial | shape_candidate | 2 | radial | cos | 0.135517 | 3.189279e+05 | yes |
| m2_cos_tangent | shape_candidate | 2 | tangent | cos | 0.130453 | 3.302646e+05 | yes |
| m2_sin_radial | shape_candidate | 2 | radial | sin | 0.130585 | 3.302649e+05 | yes |
| m2_sin_tangent | shape_candidate | 2 | tangent | sin | 0.134989 | 3.187952e+05 | yes |
| m3_cos_radial | shape_candidate | 3 | radial | cos | 0.155120 | 2.813387e+05 | yes |
| m3_cos_tangent | shape_candidate | 3 | tangent | cos | 0.155174 | 2.813368e+05 | yes |
| m3_sin_radial | shape_candidate | 3 | radial | sin | 0.159924 | 2.823702e+05 | yes |
| m3_sin_tangent | shape_candidate | 3 | tangent | sin | 0.159604 | 2.823505e+05 | yes |
| m4_cos_radial | shape_candidate | 4 | radial | cos | 0.179922 | 2.514793e+05 | yes |
| m4_cos_tangent | shape_candidate | 4 | tangent | cos | 0.188179 | 2.501655e+05 | yes |
| m4_sin_radial | shape_candidate | 4 | radial | sin | 0.188510 | 2.501655e+05 | yes |
| m4_sin_tangent | shape_candidate | 4 | tangent | sin | 0.179664 | 2.514909e+05 | yes |
| m5_cos_radial | shape_candidate | 5 | radial | cos | 0.200796 | 2.397896e+05 | yes |
| m5_cos_tangent | shape_candidate | 5 | tangent | cos | 0.200444 | 2.397872e+05 | yes |

## Interpretation

If shape slopes are near 2, a smooth tensor/Hessian shape cone may be viable. If slopes are well below 2, the local proof must be stratified rather than an ordinary tensor Hessian.

This scaffold uses floating marching-squares estimates for perturbed lengths and is not an interval certificate.

If n=10 shows stable positive deficit and a coherent scaling class, repeat for n=11-14, then interval-harden the first case whose structure is representative.

## Source Boundary

- Reference artifact SHA-256: `2b72e052aa7200f7ac5d40992843601988de234093c44860fde99a7871e19581`
