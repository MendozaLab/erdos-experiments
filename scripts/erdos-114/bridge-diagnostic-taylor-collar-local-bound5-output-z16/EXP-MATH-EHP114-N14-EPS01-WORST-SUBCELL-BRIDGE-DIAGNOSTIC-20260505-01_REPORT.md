# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_CANDIDATE_ERROR_UNDER_BUDGET_NOT_THEOREM`
- Worst subcell: `(6, 4)`
- Candidate level boxes: `15976`
- z-subdivision per marching cell: `16`
- Derivative mode: `recurrent`
- Regularity strategy: `taylor_collar`
- Local bound subdivision: `5`
- Uncertain corner cells: `4653`
- Interval-rejected boxes: `17057`
- Collar-rejected boxes: `79095`
- Regularity unresolved cells: `0`
- Recurrent unresolved cells: `null`
- Root-factored unresolved cells: `null`
- Regularity resolution delta: `null`
- Minimum candidate gradient lower: `2.6250820810208424`
- Maximum candidate Hessian upper: `4758.03929434784`
- Sum candidate normal-drift error: `2.5603409719563692`
- Sum candidate relative-length error: `7.241193257048661`
- Available exact-length bridge budget: `2.5620009612530126`

- Normal error / budget: `0.9993520731171657`
- Relative error / budget: `2.8263819438643645`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

Turn the candidate contour-error estimate into a theorem and preserve the same constants in a formal row record.
