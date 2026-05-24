# EHP114 n=14 Root-Collar Taylor/Affine Pilot

Experiment: `EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01`

Source: `EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02`

## Verdict

- Status: `"ROOT_COLLAR_TAYLOR_AFFINE_PILOT_DEPENDENCY_RECOVERED_NO_COLLAR"`
- Processed pieces: `10`
- Source unprocessed pieces: `4086`
- Source accepted length: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- Source derivative signed count: `10`
- Taylor derivative signed count: `10`
- Taylor tighter than Bernstein count: `10`
- Root-collar candidate count: `0`
- Remaining non-candidate pieces: `10`

## Interpretation

This pilot tests whether a centered Taylor/affine local model recovers dependency information that the bivariate Bernstein hull lost. It uses midpoint root-affine roots and sampled Hessian inflation as a diagnostic. It is therefore proof-facing triage, not a validated proof certificate. A useful outcome is either a strict root-collar candidate on at least one worst piece, or a clean failure that sends the route toward a hand analytic collar theorem.

## Next Blocker

Taylor/affine coordinates recover dependency but do not yet produce a strict collar. Next step is a validated Taylor model with root-affine uncertainty or a hand analytic root-collar lemma.

## Claim Ceiling

Root-collar Taylor/affine pilot only. Not a proof of Erdős #114, not a global n=14 proof, and not an exact lemniscate-length certificate.
