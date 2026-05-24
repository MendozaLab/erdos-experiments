# EHP114 n=14 Krawczyk Tube Resolver

Experiment: `EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01`

## Verdict

- Status: `"KRAWCZYK_TUBE_FAIL_DEPTH_LIMIT"`
- Source accepted length: `20.316451752723314`
- Resolved tube length upper: `-0.0`
- Total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source unresolved tubes: `165782`
- Certified pieces: `0`
- Krawczyk certified pieces: `0`
- Excluded pieces: `60796`
- Remaining unresolved pieces: `270768`
- Recursion nodes: `497346`
- Max depth reached: `1`

## Interpretation

This run only targets the remaining pieces from the Bernstein tube artifact. It uses a one-dimensional interval Krawczyk inclusion in the active chart direction. It does not rerun global slab length and it does not claim a proof of Erdős #114. A pass requires zero remaining unresolved pieces and total length under the cap.

## Next Blocker

Krawczyk contraction did not isolate every remaining piece. Next method should use tighter root-parameter intervals, per-piece root-affine narrowing, or full bivariate Bernstein/Krawczyk subdivision.
