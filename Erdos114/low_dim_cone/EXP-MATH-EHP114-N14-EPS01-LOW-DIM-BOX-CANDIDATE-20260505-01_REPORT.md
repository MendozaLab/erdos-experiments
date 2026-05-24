# EHP114 n=14 eps=0.1 Two-Dimensional Box-Candidate Probe

Experiment: `EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01`

Parent packet: `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01`

## Meaning

This run densifies the eps=0.1 two-dimensional spectral span where the previous
low-dimensional probe had full root-admissibility. It is a bridge toward an
interval-box proof, but it is still point-evaluation evidence.

## Verdict

- Status: `EPS01_LOW_DIM_BOX_CANDIDATE_PASS`
- Evaluated oracle points: `361`
- Failure count: `0`
- Minimum point margin: `12.05951541635811`
- Candidate cells whose four corners and center all pass: `164`
- Minimum sample margin among candidate cells: `12.062197693283327`
- Continuous certificate: `False`

## Why This Is Not Yet A Box Certificate

The current oracle evaluates concrete coefficient points. It does not enclose the full (u0,u1) coefficient boxes. A rigorous box certificate still needs an interval lower bound or a Lipschitz/Cauchy variation bound over each cell.

## Next Blocker

Add a rigorous cell-wise variation bound for D14 over (u0,u1) coefficient boxes, or upgrade the Rust oracle to accept interval coefficient boxes.

## Claim Ceiling

Dense eps=0.1 2D spectral-span box-candidate evidence. This is not a continuous interval-box certificate, not a Lean proof, and not a proof of Erdős #114.
