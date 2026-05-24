# EHP114 n=14 eps=0.1 One-Cell Variation Probe

Experiment: `EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01`

## Meaning

This run selects the strongest sampled eps=0.1 coefficient cell and tests the
first bridge from point evidence toward a cell certificate.

## Verdict

- Status: `ONE_CELL_ROOT_BOX_PASS_EMPIRICAL_VARIATION_PASS`
- Root box admissible by affine radius bound: `True`
- Max affine root-radius upper bound: `0.9936164637647286`
- Oracle grid points inside selected cell: `289`
- Deficit failures on sampled points: `0`
- Minimum sampled margin: `12.101568027485637`
- Maximum sampled margin: `12.103239066384457`
- Max neighbor margin slope: `1.4657307798578183`
- Required Lipschitz upper bound for this margin: `9779.543788829036`
- Empirical margin after one cell-radius variation: `12.099754278181432`
- Continuous deficit certificate: `False`

## Why This Is Not Yet A Deficit Certificate

The affine root-radius bound encloses root admissibility for the selected coefficient cell. The deficit lower bound is still sampled at grid points. Certifying it continuously requires a rigorous Lipschitz/Cauchy variation bound or an interval-coefficient length oracle.

## Next Blocker

Turn the empirical margin-variation diagnostic into a rigorous Lipschitz/Cauchy bound over the selected cell.

## Claim Ceiling

One-cell root-box enclosure plus dense deficit variation probe at eps=0.1. This is not a continuous deficit certificate, not a Lean proof, and not a proof of Erdős #114.
