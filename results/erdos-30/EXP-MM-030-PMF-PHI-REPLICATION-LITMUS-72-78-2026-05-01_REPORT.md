# EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-78-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_PHI_REPLICATION_LITMUS / EXACT_PACKET_BACKED / INTERPRETIVE

## Question

Does the post-72 branch fan show phi/Fibonacci-mediated self-replication, or only a more general inherited-plus-translated replication mechanism?

## Result

The finite window supports inherited-plus-translated self-replication, but it does not support phi/Fibonacci mediation at the current tolerance.

## Checks

- every post-72 face contains previous face union previous face shifted by `+1`: `True`
- finite ratios are within tolerance of phi: `False`
- Fibonacci residuals are small relative to exact face counts: `False`
- skeleton ratios are within tolerance of phi: `False`

## Replication Table

| n | face | inherited union | inherited overlap | new face | skeletons | new skeletons | face ratio | skeleton ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 72 | 4 | 0 | 0 | 4 | 2 | 2 | NA | NA |
| 73 | 8 | 8 | 0 | 0 | 2 | 0 | 2.000000 | 1.000000 |
| 74 | 34 | 12 | 4 | 22 | 13 | 11 | 4.250000 | 6.500000 |
| 75 | 84 | 60 | 8 | 24 | 25 | 12 | 2.470588 | 1.923077 |
| 76 | 214 | 134 | 34 | 80 | 65 | 40 | 2.547619 | 2.600000 |
| 77 | 482 | 344 | 84 | 138 | 134 | 69 | 2.252336 | 2.061538 |
| 78 | 970 | 750 | 214 | 220 | 244 | 110 | 2.012448 | 1.820896 |

## Claim Boundary

This packet can support self-replication language in the finite exact-face sense. It does not support phi/Fibonacci mediation unless the explicit ratio and recurrence checks pass. No asymptotic claim and no Erdős #30 proof claim are made.
