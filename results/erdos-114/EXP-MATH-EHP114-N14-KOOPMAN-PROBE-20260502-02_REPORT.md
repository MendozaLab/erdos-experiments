# EXP-MATH-EHP114-N14-KOOPMAN-PROBE-20260502-02 Report

## Claim Tested

Bind the n=14 EHP interval artifact to the Lean-side MDL regime rule:

- `spectral_gap > exp(1)` -> `QuantumShadow`
- otherwise -> `Classical`

## Result

- Degree: 14
- Spectral gap: 0.869067352436
- Threshold exp(1): 2.718281828459
- gap / exp(1): 0.319712
- MDL regime under this probe: `Classical`
- Spectral entropy: 0.358325452335
- Effective unitary dimension: 1
- Sampled perimeter: 54.354928372685 (`diagnostic_only`)
- alpha_star: 0.506329113924
- P1 floor collapse: True
- P3 entropy curvature: True

## Interval Artifact Binding

- Proof verdict: `EHP_N14_PROVEN`
- Interval proof complete: True
- Rigor: `ieee_1788_interval_arithmetic_inari`
- L* lower: 30.852910841548532
- L* upper: 30.852910841548546
- L* midpoint: 30.85291084154854
- Sampled-perimeter relative error vs L* midpoint: 0.761744

## Interpretation

The one-degree n=14 Koopman-kernel probe does not cross exp(1); under this operator, n=14 classifies as Classical.

The n=14 run preserves the Leg-4 floor signature: P1 and P3 are satisfied for z^14-1.

This is a sampled kernel measurement, not a theorem about the true Koopman spectrum. It binds n=14 to the existing probe rule and the existing interval perimeter artifact. The reused Leg-4 perimeter estimator is diagnostic only for n=14; the interval artifact is the perimeter truth surface.

## Source Artifacts

- Leg-4 script: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdosatlas-workbench/experiments/leg4_unitarity_114.py`
- Leg-4 script SHA-256: `eebf630301d4846333844c199d9d3471fa81364ee90456dfbbbea39f24778e0c`
- n=14 interval source: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json`
- n=14 interval SHA-256: `50b1c965c842ced25b2930c2b71ffb6e2da693872aa464a19fbd9d5d9efa0ca7`

## Status Boundary

This supports `N14_CLASSICAL_BY_CURRENT_KOOPMAN_GAP` and `N14_PHASE_FLOOR_PRESENT`.
It does not prove the true Koopman spectrum is Classical. The sampled perimeter is not
used as the perimeter truth surface.
