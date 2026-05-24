# EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-64-2026-05-01

Date: 2026-05-01
Problem: Erdos #30
Status: DERIVED_PHI_REPLICATION_LITMUS / EXACT_PACKET_BACKED / INTERPRETIVE

## Question

Does the post-58 branch fan show phi/Fibonacci-mediated self-replication, or only a more general inherited-plus-translated replication mechanism?

## Result

The finite window supports inherited-plus-translated self-replication, but it does not support phi/Fibonacci mediation at the current tolerance.

## Checks

- every post-58 face contains previous face union previous face shifted by `+1`: `True`
- finite ratios are within tolerance of phi: `False`
- Fibonacci residuals are small relative to exact face counts: `False`
- skeleton ratios are within tolerance of phi: `False`

## Replication Table

| n | face | inherited union | inherited overlap | new face | skeletons | new skeletons | face ratio | skeleton ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 58 | 10 | 0 | 0 | 10 | 2 | 2 | NA | NA |
| 59 | 18 | 14 | 6 | 4 | 4 | 2 | 1.800000 | 2.000000 |
| 60 | 54 | 26 | 10 | 28 | 18 | 14 | 3.000000 | 4.500000 |
| 61 | 152 | 90 | 18 | 62 | 49 | 31 | 2.814815 | 2.722222 |
| 62 | 398 | 250 | 54 | 148 | 123 | 74 | 2.618421 | 2.510204 |
| 63 | 1022 | 644 | 152 | 378 | 312 | 189 | 2.567839 | 2.536585 |
| 64 | 2360 | 1646 | 398 | 714 | 669 | 357 | 2.309198 | 2.144231 |

## Claim Boundary

This packet can support self-replication language in the finite exact-face sense. It does not support phi/Fibonacci mediation unless the explicit ratio and recurrence checks pass. No asymptotic claim and no Erdős #30 proof claim are made.
