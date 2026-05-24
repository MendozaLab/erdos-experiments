# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_REGULARITY_INTERVAL_UNRESOLVED`
- Worst subcell: `(6, 4)`
- Candidate level boxes: `20935`
- z-subdivision per marching cell: `16`
- Derivative mode: `recurrent`
- Regularity strategy: `level_set_collar`
- Uncertain corner cells: `6099`
- Interval-rejected boxes: `17057`
- Collar-rejected boxes: `74136`
- Regularity unresolved cells: `4443`
- Recurrent unresolved cells: `null`
- Root-factored unresolved cells: `null`
- Regularity resolution delta: `null`
- Minimum candidate gradient lower: `0.0`
- Maximum candidate Hessian upper: `8137.030701453701`
- Sum candidate normal-drift error: `29.5787588492358`
- Sum candidate relative-length error: `83.66055545944704`
- Available exact-length bridge budget: `2.5620009612530126`

- Normal error / budget: `11.545178669554264`
- Relative error / budget: `32.654380979830194`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

Interval boxes for p' still include zero on some active cells. Subdivide z-cells, use Bernstein/affine arithmetic, or restrict to a validated level-set collar before any length bridge can be proof-grade.
