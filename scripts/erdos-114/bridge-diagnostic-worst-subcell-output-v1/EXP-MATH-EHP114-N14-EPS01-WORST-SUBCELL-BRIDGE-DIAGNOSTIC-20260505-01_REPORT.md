# EHP114 n=14 Worst-Subcell Bridge Diagnostic

Experiment: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

Source: `EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01`

## Verdict

- Status: `BRIDGE_REGULARITY_INTERVAL_UNRESOLVED`
- Worst subcell: `(6, 4)`
- Active marching cells: `438`
- Uncertain corner cells: `138`
- Regularity unresolved cells: `437`
- Minimum candidate gradient lower: `0.0`
- Maximum candidate Hessian upper: `344393.5146641399`
- Sum candidate normal-drift error: `65.54909491880817`
- Sum candidate relative-length error: `12.080735767210504`
- Available exact-length bridge budget: `2.5620009612530126`

## Claim Ceiling

This is a bridge diagnostic only. It does not validate exact lemniscate length, does not prove the marching-squares-to-exact comparison theorem, and does not prove Erdős #114.

## Next Blocker

Interval boxes for p' still include zero on some active cells. Subdivide z-cells, use Bernstein/affine arithmetic, or restrict to a validated level-set collar before any length bridge can be proof-grade.
