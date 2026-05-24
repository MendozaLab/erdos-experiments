# EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01 Report

## Meaning

This is a budget prototype for the exact-length lift after the Python SUBDIV8 pass.
It does not compute exact lemniscate length and does not certify the marching-squares oracle as proof-grade length.

## Source

- Source experiment: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01`
- Source path: `erdos-experiments/Erdos114/low_dim_cone/EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.json`
- Source failures: `0`

## Budget

- `eps`: `0.1`
- `lstar_lower`: `30.852910841548532`
- `target_reserve`: `10.180114778928864`
- uniform exact length cap: `20.672796062619668`
- source max marching length upper: `18.110795101362747`
- minimum additive exact-length error budget: `2.5620009612569206`
- mean additive exact-length error budget: `3.3689587217741392`
- minimum relative error budget vs marching upper: `0.14146264407044962`

Worst additive-budget cell:

- subcell: `(6, 4)`
- `u0_interval`: `[-0.00043749999999999995, -0.00021874999999999998]`
- `u1_interval`: `[0.0008749999999999999, 0.0010937499999999999]`
- marching length upper: `18.110795101362747`
- allowed additive exact-length error: `2.5620009612569206`

## Sufficient Bridge Condition

For every SUBDIV8 subcell C, prove a validated exact lemniscate length upper enclosure L_exact(C) <= L_ms(C) + E(C), with E(C) <= margin_lower(C). A uniform additive bridge with E <= 2.56200096125692 would preserve the current scalar reserve on the whole selected cell.

## Missing For Proof Grade

- a lower bound on |grad(|p|)| or |p'| on the relevant level-set collar
- a validated collar width tau around |p| = 1
- an interval enclosure for the area of the collar or a direct implicit-curve length enclosure
- a proof that the numerical contour extraction overcounts exact length by at most E(C)

Claim ceiling: budget-only prototype; not exact lemniscate certification, not a Lean theorem, and not a proof of Erdos #114.
