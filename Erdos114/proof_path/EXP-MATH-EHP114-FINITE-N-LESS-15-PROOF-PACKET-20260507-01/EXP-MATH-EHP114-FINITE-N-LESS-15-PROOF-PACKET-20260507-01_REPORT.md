# EHP114 Finite Proof Packet For n < 15

Experiment: `EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01`

## Verdict

Status: `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`.

This packet assembles the finite proof frontier for positive degrees below 15. It adds the direct `n=1` analytic row, imports `n=2` as a literature row, and revalidates the DOI-backed Rust/inari rows for `n=3..14` by SHA.

## Counts

- Degree range: `1 <= n < 15`
- Covered degrees: `14`
- Missing degrees: `0`
- Trivial analytic rows: `1`
- Literature-only rows: `1`
- Local interval-certified rows: `12`
- Source SHA failures: `0`

## Important Caveat

The `n=13` row remains accepted as part of the DOI-backed local interval surface but carries a route exception because its stored branch-and-bound profile has zero evaluations and an empty level list.

## Claim Ceiling

This packet proves only the finite range n < 15 under the stated literature and interval-certificate dependencies. It does not prove Erdos #114 for all n and does not bridge Tao's opaque high-degree threshold.
