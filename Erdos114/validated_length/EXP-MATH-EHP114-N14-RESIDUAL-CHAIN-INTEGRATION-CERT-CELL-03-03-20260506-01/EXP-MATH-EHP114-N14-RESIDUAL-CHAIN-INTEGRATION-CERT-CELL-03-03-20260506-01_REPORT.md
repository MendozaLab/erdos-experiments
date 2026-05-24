# EHP114 n=14 Residual-Chain Integration Certificate

Experiment: `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-03-03-20260506-01`

## Verdict

- Status: `RESIDUAL_CHAIN_FAIL_COUNT_DRIFT`
- Regular regions closed by L24: `8`
- Critical regions closed by L27 partition closure: `9`
- Critical regions closed by L28 root-location: `19`
- Critical regions closed by L29 sharp root-location: `24`
- Final critical closed count: `52` of `52`
- Total residual closure count: `60`
- Ownership duplicate count: `0`
- Source-filter mismatch count: `0`
- Candidate-count drift count: `1`
- Total validated length upper: `18.01426357052187`
- Exact length cap: `20.672796062619668`
- Margin to cap: `2.658532492097798`
- First failed condition: `L24 processed_regular_region_count is not 8`

## Interpretation

This certificate does not recompute the analytic inequalities. It verifies that the already emitted L24, L27, L28, and L29 local certificates compose into a single residual-chain closure for the n=14 hard cell `(6,4)`. The pass condition is bookkeeping integrity: source checksums, source filters, candidate counts, ownership uniqueness, and length-cap consistency must all agree.

## Claim Ceiling

Local n=14 hard-cell residual-chain integration certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
