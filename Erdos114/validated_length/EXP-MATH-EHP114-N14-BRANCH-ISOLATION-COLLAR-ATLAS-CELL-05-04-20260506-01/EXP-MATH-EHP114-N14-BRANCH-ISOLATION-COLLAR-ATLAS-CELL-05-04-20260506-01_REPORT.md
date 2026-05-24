# EHP114 n=14 Branch Isolation Collar Atlas

Experiment: `EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-CELL-05-04-20260506-01`

Source: `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03`

## Verdict

- Status: `"BRANCH_ISOLATION_FAIL_COLLAR"`
- Processed unresolved branches: `16`
- Certified branches: `0`
- Excluded regions: `0`
- Remaining unresolved branches: `64`
- Source accepted length: `19.71332226297289`
- Resolved unresolved-branch length upper: `-0.0`
- Total validated length upper: `19.71332226297289`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.9594737996467764`
- Ownership duplicate count: `0`
- First failed condition: `"collar_wall_sign_separation_failed"`

## Interpretation

This run targets only unresolved branch tubes from the z32 slab artifact. It uses exclusion first, then a monotone collar lemma: a fixed nonzero derivative component plus opposite signed collar walls certifies one owned graph branch. It does not rerun global slab length.

## Next Blocker

Collar wall sign separation did not close every processed tube. Need sharper Taylor remainder bounds or analytic root-collar inequalities.

## Claim Ceiling

Local n=14 hard-cell branch-isolation diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
