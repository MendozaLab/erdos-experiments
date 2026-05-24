# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_CANDIDATE_ERROR_OVER_BUDGET`
- Worst subcell: `(6, 4)`
- Candidate level boxes: `20935`
- z-subdivision per marching cell: `16`
- Derivative mode: `recurrent`
- Regularity strategy: `taylor_collar`
- Uncertain corner cells: `6099`
- Interval-rejected boxes: `17057`
- Collar-rejected boxes: `74136`
- Regularity unresolved cells: `0`
- Recurrent unresolved cells: `null`
- Root-factored unresolved cells: `null`
- Regularity resolution delta: `null`
- Minimum candidate gradient lower: `2.4845471318753294`
- Maximum candidate Hessian upper: `8137.030701453701`
- Sum candidate normal-drift error: `5.606367282955769`
- Sum candidate relative-length error: `15.856460244499416`
- Available exact-length bridge budget: `2.5620009612530126`

- Normal error / budget: `2.1882768069742764`
- Relative error / budget: `6.189092230763413`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

The simple condition-ratio contour-error candidates are too coarse. Need sharper local charting, smaller z-cells, or a coarea/implicit-function bound with better constants.
