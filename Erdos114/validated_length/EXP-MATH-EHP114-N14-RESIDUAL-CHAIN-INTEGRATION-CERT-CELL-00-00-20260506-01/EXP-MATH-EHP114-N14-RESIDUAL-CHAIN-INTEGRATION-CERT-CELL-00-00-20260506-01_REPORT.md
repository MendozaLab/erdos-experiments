# EHP114 n=14 Residual-Chain Integration Certificate

Experiment: `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-00-00-20260506-01`

## Verdict

- Status: `RESIDUAL_CHAIN_FAIL_COUNT_DRIFT`
- Regular regions closed by L24: `8`
- Critical regions closed by L27 partition closure: `14`
- Critical regions closed by L28 root-location: `21`
- Critical regions closed by L29 sharp root-location: `21`
- Final critical closed count: `56` of `56`
- Total residual closure count: `64`
- Ownership duplicate count: `0`
- Source-filter mismatch count: `2`
- Candidate-count drift count: `4`
- Total validated length upper: `20.316451752723314`
- Exact length cap: `20.672796062619668`
- Margin to cap: `0.3563443098963539`
- First failed condition: `L24 total_validated_length_upper drifted: 15.410486279252329`

## Interpretation

This certificate does not recompute the analytic inequalities. It verifies that the already emitted L24, L27, L28, and L29 local certificates compose into a single residual-chain closure for the n=14 hard cell `(6,4)`. The pass condition is bookkeeping integrity: source checksums, source filters, candidate counts, ownership uniqueness, and length-cap consistency must all agree.

## Claim Ceiling

Local n=14 hard-cell residual-chain integration certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
