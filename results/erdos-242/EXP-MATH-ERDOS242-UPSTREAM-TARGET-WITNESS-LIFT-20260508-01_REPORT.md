# Erdos #242 Upstream Target Witness Lift

Experiment: `EXP-MATH-ERDOS242-UPSTREAM-TARGET-WITNESS-LIFT-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `UPSTREAM_TARGETS_WITNESS_LIFTED`

## Meaning

This gate moves from residue-pair overlap to target-level meaning. It asks whether a hard-strip prime hit by Mihnea/Bogdan's pinned upstream filter table can be lifted locally to an explicit Rosati witness with verified denominators.

## Upstream Target Scan

- Pinned commit: `e36eef1815d339701b9f168fe7fa504ccfa401e8`
- Filter rows: `148923`
- Upstream moduli: `148923`
- Upstream `(modulus,residue)` pairs: `12199540`
- Filters SHA-256: `4aafb1e2cd909439be5a32c84bff5d6cd27fdce95a91edaee219de237f2c0ffd`
- SHA matches prior overlap gate: `True`
- Upstream-filtered hard-strip targets: `9732`
- Total upstream target-hit records: `334567`
- Upstream hit validation errors: `0`

## Local Witness Lift

- Max n: `1000000`
- Constant bound: `21`
- Hard-strip targets: `9732`
- Locally witnessed targets: `9732`
- Upstream-filtered and locally witnessed: `9732`
- Upstream-filtered but local witness missing: `0`
- Local witness but no observed upstream hit: `0`
- Invalid emitted witnesses: `0`
- Local certificates: `357049`
- Local filter pairs: `7626`

## Supporting Pair Context

- Prior pair intersection: `2823`
- Prior local-pair overlap fraction: `0.37018095987411487`
- Prior upstream-pair overlap fraction: `0.00023140216762271363`

## Boundary

This is not `10^18` replication. It is not exact pipeline equality. It is not a proof of Erdős-Straus. It is a local certificate-enrichment layer over upstream filter hits.
