# EHP114 n=14 Validated Taylor Root-Collar Pilot

Experiment: `EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01`

## Verdict

- Status: `"VALIDATED_TAYLOR_ROOT_COLLAR_DERIVATIVE_SIGNED_NO_COLLAR"`
- Processed pieces: `10`
- Validated derivative signed pieces: `10`
- Strict collar candidates: `0`
- Remaining non-candidate pieces: `10`
- Source accepted length: `20.316451752723314`
- Candidate length upper sum: `-0.0`
- Exact length cap: `20.672796062619668`
- Margin to cap if candidates promoted: `0.3563443098963539`

## Interpretation

This run carries the root-affine uncertainty explicitly through interval roots and interval Taylor/Hessian enclosures. It is stronger than the midpoint Taylor/affine diagnostic, but it is still a local pilot over a small prefix of pieces. A strict collar candidate requires an interval Newton/Krawczyk image to be a strict subset of the original dependent interval.

## Next Blocker

Root-affine Taylor enclosures keep derivative signs but Krawczyk is not a strict subset. The next proof-facing step is a sharper analytic collar lemma or higher-order Taylor model with explicit third-derivative remainder.

## Claim Ceiling

Validated Taylor root-collar pilot only. Not a proof of Erdos #114, not a global n=14 proof, and not an exact lemniscate-length certificate.
