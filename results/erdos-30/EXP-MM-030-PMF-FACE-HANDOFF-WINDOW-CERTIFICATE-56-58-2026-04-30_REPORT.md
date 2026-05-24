# EXP-MM-030-PMF-FACE-HANDOFF-WINDOW-CERTIFICATE-56-58-2026-04-30

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_WINDOW_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Derived from `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`. The source packet SHA was verified from the packet directory before derivation.

## Window Result

The exported exact faces for `56 <= n <= 58` are all complete under the source packet cap. Every exported witness passes a direct Sidon check. Prefix, mass, joint, and Pareto roles recomputed from coordinates match the source packet for all three rows.

| n | h(n) | exact maximizers | exported | field split | prefix winners | mass winners | joint winners | Pareto minima |
|---:|---:|---:|---:|---|---|---|---|---|
| 56 | 10 | 4 | 4 | false | [3] | [3] | [3] | [3] |
| 57 | 10 | 6 | 6 | true | [4,5] | [3] | [5] | [3,5] |
| 58 | 10 | 10 | 10 | true | [2,8,9] | [7] | [7] | [7,9] |

## Handoff Pattern

- `n = 56`: split absent; prefix, mass, joint, and Pareto all select index `3`.
- `n = 57`: split appears; four maximizers persist from `n = 56`, and two new endpoint-shift states appear. Mass selects remnant index `3`; joint selects endpoint-shift index `5`.
- `n = 58`: split persists; all six `n = 57` maximizers persist, and four new endpoint-shift states enter.

## Certificate Checks

- all rows exported all witnesses: `true`
- all exported witnesses are Sidon: `true`
- all recomputed field/Pareto roles match packet: `true`
- n=57 exposed mass/joint symmetric-difference distance: `18`

## Claim Boundary

This packet certifies a finite exact-face window phenomenon. It does not prove the Sidon asymptotic, does not make Maxwell's Litmus public theorem language, and does not upgrade PMF to a proof.
