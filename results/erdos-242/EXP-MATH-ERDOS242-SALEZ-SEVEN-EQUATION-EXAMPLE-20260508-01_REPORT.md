# Erdos #242 Salez Seven-Equation Example Audit

Experiment: `EXP-MATH-ERDOS242-SALEZ-SEVEN-EQUATION-EXAMPLE-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `SALEZ_SEVEN_EXAMPLE_VERIFIED`

## Meaning

This is the first direct alignment with Salez's seven reference-equation layer. It does not implement the full sieve, but it encodes the seven equation labels and verifies Salez's Example 1 family where all seven are present.

## Result

- Family checked: `p = 120*t - 23; hence p = 1 mod 24 and p = 2 mod 5`
- t range: `[1, 1000]`
- Failure count: 0
- Verified counts by equation: `{'eqmod1a': 1000, 'eqmod1b': 1000, 'eqmod1c': 1000, 'eqmod2a': 1000, 'eqmod2b': 1000, 'eqmod2c': 1000, 'eqmod2d': 1000}`

## Boundary

This advances the #242 packet from basic-formula alignment to seven-equation example alignment. The next step is still larger: implement each reference equation as a general constant-coefficient search family and compare its coverage with the local residue operator.
