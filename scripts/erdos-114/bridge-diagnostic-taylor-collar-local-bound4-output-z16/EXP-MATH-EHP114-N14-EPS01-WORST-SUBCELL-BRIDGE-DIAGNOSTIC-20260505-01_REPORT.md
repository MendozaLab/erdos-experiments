# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_CANDIDATE_ERROR_OVER_BUDGET`
- Worst subcell: `(6, 4)`
- Candidate level boxes: `16270`
- z-subdivision per marching cell: `16`
- Derivative mode: `recurrent`
- Regularity strategy: `taylor_collar`
- Local bound subdivision: `4`
- Uncertain corner cells: `4740`
- Interval-rejected boxes: `17057`
- Collar-rejected boxes: `78801`
- Regularity unresolved cells: `0`
- Recurrent unresolved cells: `null`
- Root-factored unresolved cells: `null`
- Regularity resolution delta: `null`
- Minimum candidate gradient lower: `2.6179810866449595`
- Maximum candidate Hessian upper: `4925.7474439382095`
- Sum candidate normal-drift error: `2.6981605337895673`
- Sum candidate relative-length error: `7.631001636111534`
- Available exact-length bridge budget: `2.5620009612530126`

- Normal error / budget: `1.0531457929156913`
- Relative error / budget: `2.9785319176380773`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

The simple condition-ratio contour-error candidates are too coarse. Need sharper local charting, smaller z-cells, or a coarea/implicit-function bound with better constants.
