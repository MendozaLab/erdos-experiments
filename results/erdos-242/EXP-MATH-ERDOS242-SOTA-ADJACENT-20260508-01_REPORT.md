# Erdos #242 SOTA-Adjacent Review Packet

Experiment: `EXP-MATH-ERDOS242-SOTA-ADJACENT-20260508-01`
Status: `REVIEW_ONLY`
Verdict: `SOTA_ADJACENT_REVIEW_READY`

## Meaning

This packet does not compete with the current global verification frontier. It makes the #242 modular-coverage framing SOTA-adjacent by giving it a literature-aware boundary, a compiled scaling lemma, and a stronger local coverage audit.

## Local Signal

- Structured operator coverage up to n = 100000: 98769 / 99998 = 0.987710
- Matched random divisor-descent baseline: p50 = 0.946149, p95 = 0.981540
- Bounded split recovery with max_a = 100: 99997 / 99998 = 0.999990
- First unresolved sampled values after bounded split: `[2521]`

## Formal Infrastructure

- Lean file: `erdos-experiments/Erdos242/lean/Erdos242Scaling.lean`
- Checked theorem: `Erdos242Scout.erdosStrausClearedSolves_mul`
- Compile command: `lake env lean erdos-experiments/Erdos242/lean/Erdos242Scaling.lean`
- Status: individual Lean compile PASS, 0 sorry in this file.

## SOTA Boundary

- Global computational frontier is not ours: Mihnea-Bogdan report `10^18` verification.
- Modular-equation frontier is not ours: Salez reports seven modular equations and `10^17` checking.
- Our adjacent niche is operator/formalization: residue coverage + divisor descent + Lean scaling infrastructure.

## Next Gate

To move beyond adjacent, the packet needs a reviewed reconstruction of Salez/Mordell/Terzi-style modular equations and a direct comparison showing whether the operator recovers, simplifies, or formalizes any of those equations.
