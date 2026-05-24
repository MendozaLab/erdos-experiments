# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_CANDIDATE_ERROR_UNDER_BUDGET_NOT_THEOREM`
- Worst subcell: `(6, 4)`
- Candidate level boxes: `15793`
- z-subdivision per marching cell: `16`
- Derivative mode: `recurrent`
- Regularity strategy: `taylor_collar`
- Local bound subdivision: `6`
- Uncertain corner cells: `4614`
- Interval-rejected boxes: `17057`
- Collar-rejected boxes: `79278`
- Regularity unresolved cells: `0`
- Recurrent unresolved cells: `null`
- Root-factored unresolved cells: `null`
- Regularity resolution delta: `null`
- Minimum candidate gradient lower: `2.627337117546352`
- Maximum candidate Hessian upper: `4648.132314117793`
- Sum candidate normal-drift error: `2.4736484172357147`
- Sum candidate relative-length error: `6.995992478971961`
- Available exact-length bridge budget: `2.5620009612530126`

- Normal error / budget: `0.9655142424403749`
- Relative error / budget: `2.730675196761203`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

Turn the candidate contour-error estimate into a theorem and preserve the same constants in a formal row record.
