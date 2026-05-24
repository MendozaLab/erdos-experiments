# EXP-MATH-EHP114-MDL-PROBE-20260502-01 Report

## Claim Tested

Apply the Lean-side MDL regime rule to the existing Erdős #114 Koopman-gap artifact:

- `spectral_gap > exp(1)` -> `QuantumShadow`
- otherwise -> `Classical`

This is an evidence-binding probe over existing artifacts, not a new optimizer run.

## Result

- Degrees with measured gaps: 11
- QuantumShadow crossings: 0
- Classical classifications: 11
- Max measured spectral gap: 0.974558224107 at n=10
- Threshold exp(1): 2.718281828459
- alpha_star values observed: [0.5063291139240507]
- P1 floor-collapse support: 11/11
- P3 entropy-curvature support: 11/11
- z^n-1 perimeter win fraction in Leg-4 random comparison: 0.727273

## Interpretation

No measured z^n-1 Koopman gap crosses exp(1); under the current operator/probe, EHP114 stays Classical for n=3..13.

A stable alpha_star floor is present in the prior Leg-4 artifact, with P1 and P3 satisfied across measured degrees.

This does not prove EHP114 is globally Classical. It only says the current sampled Koopman-kernel probe does not observe a QuantumShadow gap crossing for z^n-1 over n=3..13. n=14 has interval perimeter support but no bound spectral-gap measurement in this artifact set.

## Per-Degree Probe Table

| n | spectral_gap | gap/exp(1) | regime | alpha_star | z^n-1 win? | proof verdict |
|---:|---:|---:|---|---:|---|---|
| 3 | 0.651860149394 | 0.239806 | Classical | 0.506329113924 | no | EHP_N3_PROVEN |
| 4 | 0.708441732894 | 0.260621 | Classical | 0.506329113924 | yes | EHP_N4_PROVEN |
| 5 | 0.743543858081 | 0.273534 | Classical | 0.506329113924 | no | EHP_N5_PROVEN |
| 6 | 0.897455346050 | 0.330155 | Classical | 0.506329113924 | yes | EHP_N6_PROVEN |
| 7 | 0.876438375746 | 0.322424 | Classical | 0.506329113924 | no | EHP_N7_PROVEN |
| 8 | 0.876587600088 | 0.322479 | Classical | 0.506329113924 | yes | EHP_N8_PROVEN |
| 9 | 0.880554922924 | 0.323938 | Classical | 0.506329113924 | yes | EHP_N9_PROVEN |
| 10 | 0.974558224107 | 0.358520 | Classical | 0.506329113924 | yes | EHP_N10_PROVEN |
| 11 | 0.848074722786 | 0.311989 | Classical | 0.506329113924 | yes | EHP_N11_PROVEN |
| 12 | 0.929945260802 | 0.342108 | Classical | 0.506329113924 | yes | EHP_N12_PROVEN |
| 13 | 0.885752603340 | 0.325850 | Classical | 0.506329113924 | yes | EHP_N13_PROVEN |

## Source Artifacts

- Leg-4 source: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdosatlas-workbench/experiments/LEG4_UNITARITY_114_RESULTS.json`
- Leg-4 SHA-256: `933a43711fbe6b90c783d5efa22d80ae28364ba24af6dc455f4494f14d64433c`
- n=14 interval source: `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json`
- n=14 SHA-256: `50b1c965c842ced25b2930c2b71ffb6e2da693872aa464a19fbd9d5d9efa0ca7`

## Status Boundary

This supports a `CLASSICAL_BY_CURRENT_KOOPMAN_GAP` finding for n=3..13, plus a separate
`PHASE_FLOOR_PRESENT` finding. It does not support a public claim that EHP114 is globally
Classical or globally QuantumShadow.
