# Erdos #242 Salez Filter Theorem Packet

Experiment: `EXP-MATH-ERDOS242-SALEZ-FILTER-THEOREM-PACKET-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `SALEZ_FILTER_LEMMAS_FORMALIZED_LOCAL`

## Meaning

This packet turns the local #242 certificate layer into formal certificate infrastructure. The Lean file proves seven Salez-style sufficient-condition lemmas that produce cleared Erdős-Straus witnesses.

## Lean Manifest

- Lean file: `erdos-experiments/Erdos242/lean/Erdos242SalezFilters.lean`
- Lean SHA-256: `718b5afb84ec13553aab3ba7bd279f93c54d20d0315b1d4822daaf1023084a80`
- Expected theorem labels present: `7 / 7`
- Forbidden Lean tokens: `0`

## Certificate Coverage

- Target: `prime p = 1 mod 24`, `p <= 1000000`
- Constant bound: `21`
- Certified targets: `9732 / 9732`
- Missing targets: `0`
- Certificates: `357049`
- Invalid certificates: `0`
- Unknown equation labels: `{}`
- Equation counts: `{'eqmod1a': 136506, 'eqmod1b': 37780, 'eqmod1c': 46867, 'eqmod2a': 37690, 'eqmod2b': 53353, 'eqmod2c': 22330, 'eqmod2d': 22523}`

## Boundary

This is formal certificate infrastructure, not a new computational bound. It is not `10^18` replication, not exact upstream pipeline equality, and not a proof of Erdős-Straus.
