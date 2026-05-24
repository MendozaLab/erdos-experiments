# EHP114 n=14 Slab Root-Isolation Validated-Length Hard-Cell Diagnostic

Experiment: `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-07-03-20260506-01`

Source: `EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01`

## Verdict

- Status: `SLAB_FAIL_ROOT_ISOLATION`
- Degree: `14`
- eps: `0.1`
- Root-affine subcell: `(7, 3)`
- Exact length cap: `20.672796062619668`
- Outside extent excluded: `true` (margin `1.0072307982716286`)
- Slab branch count: `2413`
- Candidate y-cells: `550664`
- Excluded y-cells: `49010936`
- Unresolved branch count: `4402`
- Ownership duplicate count: `0`
- Sum branch length upper: `20.50697022838823`
- Endpoint/overlap tax: `0.0`
- Total validated length upper: `20.50697022838823`
- Margin to cap: `0.1658258342314376`

## Interpretation

This run demotes marching squares to a diagnostic and attempts a slab-wise implicit root-isolation enclosure. Each accepted branch has endpoint sign separation and a fixed nonzero `Fy` interval, so it is counted as one owned graph segment over its x-slab. The prior `relative error / budget = 2.730675196761203` is retained only as historical context and is not load-bearing for the pass/fail decision.

## Claim Ceiling

This is a local hard-cell diagnostic. It is not a proof of Erdős #114 and not a global n=14 proof.

## Next Blocker

Some candidate y-runs failed endpoint sign separation or fixed nonzero Fy. Need narrower slabs, vertical refinement, rotated charts, or interval Newton isolation.
