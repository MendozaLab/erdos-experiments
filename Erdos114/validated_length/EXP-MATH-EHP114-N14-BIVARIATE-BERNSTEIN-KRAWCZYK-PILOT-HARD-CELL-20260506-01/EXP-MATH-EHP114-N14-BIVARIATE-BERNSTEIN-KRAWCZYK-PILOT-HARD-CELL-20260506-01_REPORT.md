# EHP114 n=14 Bivariate Bernstein/Krawczyk Pilot

Experiment: `EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01`

## Verdict

- Status: `"BIVARIATE_BERNSTEIN_KRAWCZYK_PILOT_FAILS"`
- Processed pieces: `4096`
- Source unprocessed pieces: `266672`
- Source accepted length: `20.316451752723314`
- Resolved length upper: `-0.0`
- Total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source unresolved pieces: `270768`
- Certified pieces: `0`
- Krawczyk certified pieces: `0`
- Bivariate-excluded pieces: `0`
- Remaining unresolved pieces: `4096`
- Recursion nodes: `4096`
- Max depth reached: `0`

## Interpretation

This run only targets the remaining pieces from the last full Krawczyk artifact. It replaces the old edge/one-dimensional Bernstein checks with a two-dimensional Bernstein hull for `F`, `Fx`, and `Fy` on each piece, then attempts strict interval Krawczyk contraction in the active chart direction. It does not rerun global slab length and it does not claim a proof of Erdős #114. A pass requires zero remaining unresolved pieces and total length under the cap.

## Next Blocker

Bivariate Bernstein/Krawczyk did not certify or exclude material pieces. Retire this route and move to affine/Taylor model arithmetic or an analytic root-collar theorem.
