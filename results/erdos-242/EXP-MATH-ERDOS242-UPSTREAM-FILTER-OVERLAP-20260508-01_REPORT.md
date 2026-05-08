# Erdos #242 Upstream Filter-Table Overlap Gate

Experiment: `EXP-MATH-ERDOS242-UPSTREAM-FILTER-OVERLAP-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `UPSTREAM_FILTER_OVERLAP_QUANTIFIED`

## Meaning

This gate compares our certificate-bearing filter surface with Mihnea/Bogdan's exact pinned `Filters.txt` table. It quantifies overlap and non-overlap; it does not assert exact equality or computational-frontier replication.

## Upstream Table

- Pinned commit: `e36eef1815d339701b9f168fe7fa504ccfa401e8`
- Filter rows: `148923`
- Upstream moduli: `148923`
- Upstream `(modulus,residue)` pairs: `12199540`
- Filters SHA-256: `4aafb1e2cd909439be5a32c84bff5d6cd27fdce95a91edaee219de237f2c0ffd`
- Residues SHA-256: `efee180138f7f25738d6b4c6b9e311da0be17a5e2d791e1df134c9e09ab361e1`
- Residues count: `2101514`

## Local Surface

- Max n: `1000000`
- Constant bound: `21`
- Certificates: `357049`
- Local moduli: `970`
- Local `(modulus,residue)` pairs: `7626`
- Invalid certificates: `0`

## Overlap

- Pair intersection: `2823`
- Fraction of upstream pairs: `0.00023140`
- Fraction of local pairs: `0.37018096`
- Upstream-only pair sample: `[{'modulus': 31, 'residue': 0}, {'modulus': 31, 'residue': 22}, {'modulus': 31, 'residue': 24}, {'modulus': 37, 'residue': 0}, {'modulus': 37, 'residue': 5}, {'modulus': 37, 'residue': 15}, {'modulus': 37, 'residue': 18}, {'modulus': 37, 'residue': 22}, {'modulus': 37, 'residue': 23}, {'modulus': 37, 'residue': 29}, {'modulus': 37, 'residue': 32}, {'modulus': 37, 'residue': 35}]`
- Local-only pair sample: `[{'modulus': 7, 'residue': 3}, {'modulus': 7, 'residue': 5}, {'modulus': 7, 'residue': 6}, {'modulus': 11, 'residue': 7}, {'modulus': 11, 'residue': 8}, {'modulus': 11, 'residue': 10}, {'modulus': 15, 'residue': 7}, {'modulus': 15, 'residue': 13}, {'modulus': 19, 'residue': 14}, {'modulus': 19, 'residue': 15}, {'modulus': 19, 'residue': 18}, {'modulus': 20, 'residue': 13}]`

## Boundary

Local-only pairs are classified as certificate-enriched natural moduli unless later proven equivalent to upstream filters. Upstream `Residues.txt` is treated as a residual progression/search surface, not a witness table.
