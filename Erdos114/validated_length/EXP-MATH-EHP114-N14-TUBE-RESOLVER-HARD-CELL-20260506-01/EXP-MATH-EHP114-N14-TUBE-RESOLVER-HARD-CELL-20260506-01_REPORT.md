# EHP114 n=14 Targeted Tube Resolver

Experiment: `EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03`

## Verdict

- Status: `"TUBE_RESOLVER_FAIL_DEPTH_LIMIT"`
- Source accepted length: `20.316451752723314`
- Resolved tube length upper: `-0.0`
- Total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source unresolved tubes: `4398`
- Certified pieces: `0`
- Excluded pieces: `93132`
- Remaining unresolved pieces: `73080`
- Recursion nodes: `328026`
- Max depth reached: `6`

## Interpretation

This run only targets the unresolved tubes from the z32 slab artifact. It does not rerun global slab length and it does not claim a proof of Erdős #114. A pass requires zero remaining unresolved pieces and total length under the cap.

## Next Blocker

Recursive subdivision did not isolate every tube. Next method should use interval Newton/Krawczyk contraction or Bernstein form certificates on the remaining pieces.
