# EHP114 n=14 Sharp p-Prime Root-Location Diagnostic

Experiment: `EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-CELL-00-04-20260506-01`

Source: `EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-CELL-00-04-20260506-01`

## Verdict

- Status: `"SHARP_PPRIME_ROOT_LOCATION_PASS_NOT_GLOBAL_PROOF"`
- Processed remaining candidates: `28`
- p-prime root-excluded candidates: `28`
- p-prime root-near candidates: `0`
- Still unresolved candidates: `0`
- Minimum root-box distance lower: `-0.018206156016638105`
- Param subdivision: `8`
- Max local split depth: `4`
- Total validated length upper: `15.557494057270386`
- Exact length cap: `20.672796062619668`
- Margin to cap: `5.115302005349282`
- First failed condition: `"none"`

## Interpretation

This diagnostic uses the invariant that a critical point on `|p| = 1` must have `p'(z) = 0`. It preserves the L28 disk test and then tests the sharper Taylor/Rouche variation bound `r |p''(z0)| + 0.5 r^2 sup |p'''|` on targeted local splits inside the 17 unresolved boxes. It does not promote branch length.

## Claim Ceiling

Local n=14 hard-cell sharp p-prime root-location diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
