# EHP114 n=14 Parameter-Sliced Krawczyk Resolver

Experiment: `EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01`

## Verdict

- Status: `"PARAM_SLICED_KRAWCZYK_PILOT_FAILS"`
- Source accepted length: `20.316451752723314`
- Max tile resolved length upper: `0.0`
- Max tile total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source unresolved pieces: `270768`
- Processed pieces per tile: `4096`
- Parameter subdivision: `4`
- Parameter tiles: `16`
- Total certified piece outcomes: `0`
- Total excluded piece outcomes: `19130`
- Total remaining unresolved outcomes: `46406`

## Interpretation

This run targets the remaining pieces from the prior Krawczyk artifact after slicing the root-affine parameter cell. Parameter tiles are alternatives, so length is aggregated by maximum tile total, not by summing over tiles. This is not a proof of Erdős #114 and not a global n=14 proof.

## Next Blocker

Pilot did not certify any piece after root-parameter slicing. Retire this cheap route and move to full bivariate Bernstein/Krawczyk subdivision.
