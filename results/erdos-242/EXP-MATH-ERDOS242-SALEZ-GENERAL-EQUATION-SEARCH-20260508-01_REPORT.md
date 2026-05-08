# Erdos #242 Salez General-Equation Bounded Search

Experiment: `EXP-MATH-ERDOS242-SALEZ-GENERAL-EQUATION-SEARCH-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `SALEZ_GENERAL_BOUNDED_SEARCH_READY`

## Meaning

This audit reconstructs Salez's seven reference equations as bounded constant-parameter search families. For each hit it rebuilds Rosati variables `A,B,C,D`, reconstructs the denominators, and verifies the identity `4/p = 1/x + 1/y + 1/z` directly.

It is a stronger alignment than the Example 1 check, but it is not Salez's optimized sieve and it does not move the computational frontier.

## Run Shape

- Max prime checked: `1000000`
- Constant bound: `21`
- Target: `prime hard strip p = 1 mod 24`
- Odd primes: `78497`
- Basic local-family primes: `68765`
- Hard-strip primes `p = 1 mod 24`: `9732`

## Bounded Search Result

- Covered target primes: `9732 / 9732 = 1.000000`
- Equation hit counts: `{'eqmod1a': 136506, 'eqmod1b': 37780, 'eqmod1c': 46867, 'eqmod2a': 37690, 'eqmod2b': 53353, 'eqmod2c': 22330, 'eqmod2d': 22523}`
- First covered target primes: `[73, 97, 193, 241, 313, 337, 409, 433, 457, 577, 601, 673, 769, 937, 1009, 1033, 1129, 1153, 1201, 1249]`
- First uncovered target primes: `[]`

## Boundary

This makes the #242 lane genuinely Salez-facing: the seven equations are now executable search families, not just labels in a report. The remaining gap is to match Salez's optimized sieve behavior and parameter strategy, then decide whether our operator view simplifies or only repackages that known machinery.
