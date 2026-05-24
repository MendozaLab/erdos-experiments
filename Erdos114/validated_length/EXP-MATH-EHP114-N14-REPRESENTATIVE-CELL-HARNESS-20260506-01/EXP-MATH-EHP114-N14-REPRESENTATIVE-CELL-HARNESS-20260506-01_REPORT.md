# EHP114 n=14 Representative-Cell Harness

Experiment: `EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-20260506-01`

## Verdict

- Status: `REPRESENTATIVE_CELL_HARNESS_BLOCKED_REGULAR_RESIDUAL_DRIFT`
- Processed cells: `CELL-00-02, CELL-03-03`
- Pass count: `1`
- Fail count: `1`
- First failed cell: `CELL-03-03`
- First failed condition: `CELL-03-03 L24 produced 12 regular regions, but only 8 were monotone-excluded; 4 regular wall-separation failures remain, and L31 correctly reported count drift rather than promoting a packet.`
- Source SHA fail count: `0`

## Interpretation

The harness confirms that the local packet path transfers to `CELL-00-02`, including the derivative-root closure chain. It then stops on `CELL-03-03`, where the source length remains under cap after one high-slope repair, but the regular residual slice has twelve regular regions rather than the prior eight and leaves four wall-separation failures. This is a proof-pipeline generality blocker, not a length-budget failure.

## Next Dependency

L33K should generalize L24/L31 to variable regular-region counts and add a proof-facing repair for the four CELL-03-03 regular wall-separation failures before returning to representative-cell batching.

## Claim Ceiling

Representative-cell harness diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
