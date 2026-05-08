# Erdos #242 Salez Basic-Formula Alignment Audit

Experiment: `EXP-MATH-ERDOS242-SALEZ-BASIC-ALIGNMENT-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `BASIC_FORMULAS_RECOVERED`

## Meaning

The local residue operator recovers the opening Salez reduction exactly: the families for `3t - 1`, `4t - 1`, and `8t - 3` correspond to our `n mod 3 = 2`, `n mod 4 = 3`, and `n mod 8 = 5` operator families.

This is useful, but it is not enough to claim SOTA. Salez's real frontier is the later complete set of seven modular equations and the `10^17` sieve; Mihnea-Bogdan report `10^18` computational verification.

## Audit Results

- Checked Salez basic-family instances up to n = 100000: 58332
- Mismatches: 0
- Divisor-descent coverage: 98769 / 99998 = 0.987710
- Remaining hard-strip count: 1229
- First hard-strip values: `[73, 97, 193, 241, 313, 337, 409, 433, 457, 577, 601, 673, 769, 937, 1009, 1033, 1129, 1153, 1201, 1249]`
- First hard-strip primes: `[73, 97, 193, 241, 313, 337, 409, 433, 457, 577, 601, 673, 769, 937, 1009, 1033, 1129, 1153, 1201, 1249]`

## Claim Boundary

We can now say the operator is aligned with the first Salez reduction and has compiled scaling infrastructure. We cannot yet say it recovers the seven modular equations, improves the sieve frontier, or resolves the `p = 1 mod 24` core.

## Next Gate

Extract the seven Salez modular equations from the paper/program, encode each as a named operator family, and rerun this audit with a family-by-family recovery table.
