# Erdos #242 Mihnea/Bogdan Certificate Comparison

Experiment: `EXP-MATH-ERDOS242-MIHNEA-BOGDAN-CERT-COMPARISON-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `CERTIFIED_PARITY_LOCAL`

## Meaning

This is a local comparison against the Salez/Mihnea-Bogdan filter lineage on the certificate axis. It does not compete with the `10^18` computational frontier. It asks whether the filter layer can be made witness-bearing at a modest scale.

## Main Result

- Target: hard-strip primes `p = 1 mod 24` up to `1000000`
- Constant bound: `21`
- Certified targets: `9732 / 9732 = 1.000000`
- Filter-only/no-witness targets: `0`
- Invalid certificates: `0`
- Certificate count: `357049`
- Filter moduli used: `970`
- Runtime seconds: `71.003845`

## Regression And Controls

- 100k regression: `1181 / 1181` with verdict `CERTIFIED_PARITY_LOCAL`
- Near-miss lower bound control: bound `20` covers `9731 / 9732` with verdict `NEEDS_FILTER_COVERAGE`
- Restricted-modulus control: moduli `[3]` covers `0 / 1181` with verdict `NEEDS_FILTER_COVERAGE`

## Boundary

Mihnea/Bogdan remain the computational-frontier comparison point. This packet's contribution is different: each local certified target carries a concrete Rosati witness and verified denominators.
