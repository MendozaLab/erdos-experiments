# EHP114 n=14 Local Hard-Cell Certificate Packet

Experiment: `EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-01-06-20260506-01`

## Local Theorem Statement

For the requested n=14 root-affine subcell at `eps = 0.1`, the certified local hard-cell length upper assembled from the accepted slab length and the L24/L27/L28/L29 residual-chain closure is at most `15.932766113819524`. This is below the local exact cap `20.672796062619668` by margin `4.740029948800144`.

## Verdict

- Status: `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF`
- Local certificate pass: `true`
- Source SHA fail count: `0`
- L31 residual-chain failures: `0`
- Residual closures verified by L31: `64` of `64`
- Ownership duplicate count: `0`
- Source-filter mismatch count: `0`
- Candidate-count drift count: `0`
- First failed condition: `none`

## Dependency Table

| Role | Experiment |
|---|---|
| L24 regular monotone closure | `EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01` |
| L27 p-prime partition closure | `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01` |
| L28 p-prime root-location closure | `EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01` |
| L29 sharp p-prime root-location closure | `EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01` |
| L31 residual-chain integration | `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01` |

## Interpretation

This packet is a local theorem-facing wrapper around verified artifacts. It does not recompute interval geometry. It checks that the L31 integration result and all child source checksums are stable, then states the exact local claim that follows. The next proof-facing dependency is global n=14 coverage over all root-affine cells using the same residual-chain discipline.

## Claim Ceiling

Local n=14 cell certificate packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. Operator-grammar and MDL sidecars are not evidence for this packet.
