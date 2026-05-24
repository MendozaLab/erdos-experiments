# EHP114 n=14 Bernstein Tube Resolver

Experiment: `EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01`

## Verdict

- Status: `"BERNSTEIN_TUBE_FAIL_DEPTH_LIMIT"`
- Source accepted length: `20.316451752723314`
- Resolved tube length upper: `-0.0`
- Total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source unresolved tubes: `73080`
- Certified pieces: `0`
- Bernstein edge certified pieces: `0`
- Excluded pieces: `87216`
- Remaining unresolved pieces: `165782`
- Recursion nodes: `432916`
- Max depth reached: `2`

## Interpretation

This run only targets the remaining pieces from the targeted tube resolver. It uses Bernstein-form edge ranges to try to certify chart boundary signs. It does not rerun global slab length and it does not claim a proof of Erdős #114. A pass requires zero remaining unresolved pieces and total length under the cap.

## Next Blocker

Bernstein edge certificates did not isolate every remaining piece. Next method should use interval Newton/Krawczyk contraction in the active derivative direction, or promote to full bivariate Bernstein subdivision.
