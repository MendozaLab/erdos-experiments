# Erdos #242 Source-Concordance Gate

Experiment: `EXP-MATH-ERDOS242-SOURCE-CONCORDANCE-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `SOURCE_CONCORDANT_CERTIFICATE_ENRICHMENT`

## Meaning

This gate checks whether the local certificate-parity runner is aligned with the Salez / Mihnea-Bogdan filter lineage. It is not a proof, not a `10^18` replication, and not a public SOTA claim.

## Concordance Classification

- Direct Salez-equation reconstruction: `True`
- Mihnea/Bogdan-style filter analogue: `True`
- Local certificate enrichment: `True`

## Local Packet Interpreted

- Certificates at `1e6`: `357049`
- Filter moduli: `970`
- Filter residue classes: `7626`
- Certified hard-strip primes: `9732 / 9732`
- Invalid certificates: `0`
- Filter-only/no-witness targets: `0`

Interpretation: the certificates are direct seven-equation Salez reconstructions with added Rosati variables and denominators. The filter moduli/residue classes are Mihnea/Bogdan-style filter analogues, not asserted as exact equality to the upstream filter tables.

## External Reference Snapshot

- Mihnea/Bogdan repo status: `FETCHED_REFERENCE_METADATA_ONLY`
- Mihnea/Bogdan HEAD: `e36eef1815d339701b9f168fe7fa504ccfa401e8`
- Expected repo paths present: `True`
- Third-party code committed locally: `False`

## Intentional Differences

- local run stops at n <= 1,000,000; Mihnea/Bogdan report a 10^18-scale computational frontier
- local output attaches Rosati variables and denominators for certified targets
- local code reconstructs source formulas; it does not vendor or execute third-party code
- local filter moduli/residue classes are certificate-bearing natural moduli, not a byte-for-byte copy of Mihnea/Bogdan Filters.txt/Residues.txt

## Boundary

This packet supports the claim that our local runner is source-concordant certificate enrichment. It does not support claiming computational-frontier parity with Mihnea/Bogdan.
