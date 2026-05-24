# EHP114 n=14 Residual-Chain Integration Certificate

Experiment: `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-00-03-20260506-01`

## Verdict

- Status: `RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF`
- Regular regions closed by L24: `8`
- Critical regions closed by L27 partition closure: `7`
- Critical regions closed by L28 root-location: `33`
- Critical regions closed by L29 sharp root-location: `16`
- Final critical closed count: `56` of `56`
- Total residual closure count: `64`
- Ownership duplicate count: `0`
- Source-filter mismatch count: `0`
- Candidate-count drift count: `0`
- Total validated length upper: `15.415136574295897`
- Exact length cap: `20.672796062619668`
- Margin to cap: `5.257659488323771`
- First failed condition: `none`

## Interpretation

This certificate does not recompute the analytic inequalities. It verifies that the already emitted L24, L27, L28, and L29 local certificates compose into a single residual-chain closure for the n=14 hard cell `(6,4)`. The pass condition is bookkeeping integrity: source checksums, source filters, candidate counts, ownership uniqueness, and length-cap consistency must all agree.

## Claim Ceiling

Local n=14 hard-cell residual-chain integration certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
