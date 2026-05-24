# EHP114 n=14 eps=0.1 One-Cell Root-Affine Oracle Prototype

Experiment: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-ORACLE-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01`

## Meaning

This run avoids coefficient-box dependency blowup by evaluating the polynomial
directly as a product over root-affine intervals.

## Verdict

- Status: `ONE_CELL_ROOT_AFFINE_ORACLE_FAIL`
- Length upper bound: `44.5971758610202`
- Deficit lower bound: `-13.744265019471666`
- Target: `10.180114778928864`
- Margin lower bound: `-23.92437979840053`
- Active cells: `640`
- Definite-case cells: `104`
- Uncertain-corner cells: `536`
- Ambiguous-average cells: `0`

## Certificate Scope

Root-affine interval upper bound for the current marching-squares oracle functional over one selected eps=0.1 coefficient rectangle.

## Claim Ceiling

This does not certify exact lemniscate length and does not prove Erdős #114. It is a prototype continuous-cell enclosure for the existing grid oracle.

## Next Blocker

If this passes, formalize the root-affine interval oracle and then replace the grid functional with a validated exact-length enclosure. If it fails, the remaining route is analytic Cauchy/Lipschitz.
