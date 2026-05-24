# EHP114 n=14 eps=0.1 One-Cell Interval-Coefficient Oracle Prototype

Experiment: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-INTERVAL-COEFF-ORACLE-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01`

## Meaning

This run tries the first continuous-cell enclosure for the selected eps=0.1
coefficient rectangle. It propagates the full coefficient box through a
conservative interval marching-squares upper bound.

## Verdict

- Status: `ONE_CELL_INTERVAL_COEFF_ORACLE_FAIL`
- Length upper bound: `2847.9689848799253`
- Deficit lower bound: `-2817.1160740383766`
- Target: `10.180114778928864`
- Margin lower bound: `-2827.2961888173054`
- Active cells: `36920`
- Definite-case cells: `0`
- Uncertain-corner cells: `36920`
- Ambiguous-average cells: `0`

## Certificate Scope

Conservative interval-coefficient upper bound for the current marching-squares oracle functional over one selected eps=0.1 coefficient rectangle.

## Claim Ceiling

This does not certify exact lemniscate length and does not prove Erdős #114. It is a prototype continuous-cell enclosure for the existing grid oracle.

## Next Blocker

Replace the grid-oracle functional by an exact/validated lemniscate-length enclosure, or formalize why this marching-squares interval oracle is an accepted certificate layer.
