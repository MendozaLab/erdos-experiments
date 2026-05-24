# EXP-MM-030-PMF-N57-FACE-CERTIFICATE-2026-04-30

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Derived from `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`. The source packet SHA was verified from the packet directory before derivation.

## Result

The six exported exact `h = 10` maximizers at `n = 57` all pass a direct Sidon witness check. The field-winner and Pareto roles recomputed from the exported coordinates match the source packet.

## Face Classification

| index | class | inherited from n=56 | persists in n=58 | contains 57 | prefix | mass | joint | Pareto |
|---:|---|---:|---:|---|---|---|---|---|
| 0 | remnant | 0 | 0 | no | no | no | no | no |
| 1 | remnant | 1 | 1 | no | no | no | no | no |
| 2 | remnant | 2 | 4 | no | no | no | no | no |
| 3 | remnant | 3 | 5 | no | no | yes | no | yes |
| 4 | endpoint-shift | none | 6 | yes | yes | no | no | no |
| 5 | endpoint-shift | none | 7 | yes | yes | no | yes | yes |

## Certificate Checks

- all six exported witnesses are Sidon: `true`
- prefix winners match packet: `true`
- mass winners match packet: `true`
- joint winners match packet: `true`
- Pareto minima match packet: `true`
- remnant count: `4`
- endpoint-shift count: `2`
- symmetric-difference distance between exposed indices 3 and 5: `18`

## Interpretation

At `n = 57`, the exact face decomposes into four inherited remnant maximizers and two endpoint-shift maximizers. The mass field selects the remnant exposed candidate `3`; the joint field selects the endpoint-shift exposed candidate `5`. This is a finite exact-face handoff certificate, not a theorem-level Sidon result.

## Claim Boundary

Safe: finite exact-face handoff, remnant-to-endpoint field switch, first-hit accounting loses face geometry.

Unsafe: PMF proves Sidon, physics solves Erdos #30, Mendoza Limit proves the handoff, SOTA theorem result.
