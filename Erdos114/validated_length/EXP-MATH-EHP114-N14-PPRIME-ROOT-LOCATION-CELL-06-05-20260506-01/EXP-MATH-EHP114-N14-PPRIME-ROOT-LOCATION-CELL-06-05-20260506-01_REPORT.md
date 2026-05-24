# EHP114 n=14 p-Prime Root-Location Diagnostic

Experiment: `EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-CELL-06-05-20260506-01`

Source: `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-CELL-06-05-20260506-01`

## Verdict

- Status: `"PPRIME_ROOT_LOCATION_PARTIAL"`
- Processed remaining candidates: `48`
- p-prime root-excluded candidates: `23`
- p-prime root-near candidates: `25`
- Still unresolved candidates: `25`
- Minimum root-box distance lower: `-0.022968605598378003`
- Param subdivision: `8`
- Total validated length upper: `18.696441846009524`
- Exact length cap: `20.672796062619668`
- Margin to cap: `1.9763542166101438`
- First failed condition: `"2588:2:root/y0/y0 remains near a possible p-prime root under the root-free disk test"`

## Interpretation

This diagnostic uses the invariant that a critical point on `|p| = 1` must have `p'(z) = 0`. For each remaining L27 candidate and root-parameter tile, it tests whether `|p'(z0)|` dominates the spatial radius times a bound for `|p''|`; when it does, that tile is certified free of p-prime roots. It does not promote branch length.

## Claim Ceiling

Local n=14 hard-cell p-prime root-location diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
