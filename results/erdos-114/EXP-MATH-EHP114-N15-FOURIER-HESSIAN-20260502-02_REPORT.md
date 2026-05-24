# EXP-MATH-EHP114-N15-FOURIER-HESSIAN-20260502-02 Report

## Claim Tested

Can the EHP114 deficit near `z^n - 1` be seen in Hilbert/Fourier coordinates
instead of raw coefficient coordinates?

This is a diagnostic numerical probe, not a proof.

## Result

- Degree: 15
- Basis rank: 27 / expected `2n-3` = 27
- Grid resolution: 360
- eps values: [0.02, 0.01, 0.005]
- L0 marching-squares estimate: 30.465107922357
- L0 relative error vs shortcut L*: 0.072527
- Rows positive for all eps: 27/27
- Worst mode: `m7_sin_tangent` with min curvature 6.817954e+04
- Best mode: `m0_radial` with max curvature 1.569896e+06
- All tested modes positive: True

## Interpretation

Positive deficit curvature means the regular root polygon locally beats the tested Fourier perturbation direction under this floating-point estimator.

This is not an interval certificate and does not prove local maximality. Because z^n - 1 is a singular lemniscate, the measured effect is better read as a stratified deficit/unfolding signal than as an ordinary smooth Hessian certificate.

If all or most directions are positive, upgrade this to an interval or analytic Fourier-mode Hessian certificate and then bound tensor remainders.

## Singularity Sanity Check

z^n - 1 has a multiple critical point at the origin lying on the lemniscate. Small radial/constant perturbations unfold this singular level set, so an ordinary smooth Hessian interpretation is unsafe.

- Radial family checked: `p_R(z) = z^n - R^n, with R = 1 +/- eps/sqrt(n)`
- Large reference-error flag: True

| eps | radius | length | deficit vs L0 marching |
|---:|---:|---:|---:|
| 0.02 | 0.994836022205 | 9.032478457880 | 21.432629464477 |
| 0.02 | 1.005163977795 | 8.398302239560 | 22.066805682796 |
| 0.01 | 0.997418011103 | 9.980619849603 | 20.484488072754 |
| 0.01 | 1.002581988897 | 9.622609453938 | 20.842498468419 |
| 0.005 | 0.998709005551 | 10.941013303612 | 19.524094618745 |
| 0.005 | 1.001290994449 | 10.741806471789 | 19.723301450567 |

## Mode Table

| label | mode | kind | phase | min curvature | max curvature | all eps positive? |
|---|---:|---|---|---:|---:|---|
| m0_radial | 0 | radial | constant | 1.087486e+05 | 1.569896e+06 | yes |
| m1_cos_radial | 1 | radial | cos | 9.489883e+04 | 1.356412e+06 | yes |
| m1_cos_tangent | 1 | tangent | cos | 9.550020e+04 | 1.366121e+06 | yes |
| m2_cos_radial | 2 | radial | cos | 8.487327e+04 | 1.196694e+06 | yes |
| m2_cos_tangent | 2 | tangent | cos | 8.454434e+04 | 1.191081e+06 | yes |
| m2_sin_radial | 2 | radial | sin | 8.456766e+04 | 1.191187e+06 | yes |
| m2_sin_tangent | 2 | tangent | sin | 8.486964e+04 | 1.196706e+06 | yes |
| m3_cos_radial | 3 | radial | cos | 7.769194e+04 | 1.066773e+06 | yes |
| m3_cos_tangent | 3 | tangent | cos | 7.921683e+04 | 1.093848e+06 | yes |
| m3_sin_radial | 3 | radial | sin | 7.921560e+04 | 1.093848e+06 | yes |
| m3_sin_tangent | 3 | tangent | sin | 7.764386e+04 | 1.066587e+06 | yes |
| m4_cos_radial | 4 | radial | cos | 7.441607e+04 | 1.005258e+06 | yes |
| m4_cos_tangent | 4 | tangent | cos | 7.439161e+04 | 1.002736e+06 | yes |
| m4_sin_radial | 4 | radial | sin | 7.441140e+04 | 1.002727e+06 | yes |
| m4_sin_tangent | 4 | tangent | sin | 7.441017e+04 | 1.005239e+06 | yes |
| m5_cos_radial | 5 | radial | cos | 7.060412e+04 | 9.233093e+05 | yes |
| m5_cos_tangent | 5 | tangent | cos | 7.220841e+04 | 9.568219e+05 | yes |
| m5_sin_radial | 5 | radial | sin | 7.222848e+04 | 9.569175e+05 | yes |
| m5_sin_tangent | 5 | tangent | sin | 7.058561e+04 | 9.232999e+05 | yes |
| m6_cos_radial | 6 | radial | cos | 6.917565e+04 | 9.049020e+05 | yes |
| m6_cos_tangent | 6 | tangent | cos | 6.971870e+04 | 9.037934e+05 | yes |
| m6_sin_radial | 6 | radial | sin | 6.973683e+04 | 9.038387e+05 | yes |
| m6_sin_tangent | 6 | tangent | sin | 6.915884e+04 | 9.048956e+05 | yes |
| m7_cos_radial | 7 | radial | cos | 6.819814e+04 | 8.739847e+05 | yes |
| m7_cos_tangent | 7 | tangent | cos | 6.818233e+04 | 8.778573e+05 | yes |
| m7_sin_radial | 7 | radial | sin | 6.820636e+04 | 8.778016e+05 | yes |
| m7_sin_tangent | 7 | tangent | sin | 6.817954e+04 | 8.739316e+05 | yes |

## Source Boundary

- Shortcut L* artifact consulted only for reference: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/EXP-MM-EHP-007-n15-inari_RESULTS.json`
- Shortcut artifact SHA-256: `1d2cd344c258f72e1408b05c56e61407689486f4694533bd8bf5476a15782b2a`

The shortcut artifact verdict is not used as proof in this probe.
