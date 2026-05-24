# EHP114 n=14 CELL-05-06 High-Slope Slab Repair Diagnostic

Experiment: `EXP-MATH-EHP114-N14-CELL-05-06-HIGH-SLOPE-SLAB-REPAIR-20260506-01`

## Verdict

- Status: `HIGH_SLOPE_SLAB_REPAIR_ADJUSTED_UNDER_CAP_DIAGNOSTIC`
- Source total upper: `38.96678653155184`
- Exact cap: `20.672796062619668`
- Source margin: `-18.293990468932172`
- Target count: `1`
- Target original length sum: `21.26085781145849`
- Repaired target length upper: `0.0012648804210201544`
- Target length delta: `-21.25959293103747`
- Adjusted total if replaced: `17.70719360051437`
- Adjusted margin to cap: `2.9656024621052985`
- Remaining unresolved source branches: `4394`
- First failed condition: `none`

## Interpretation

This diagnostic only repairs the selected high-slope slab branches exposed by L33E. It splits each target slab in x, re-isolates candidate y-tubes inside the original branch tube, and recomputes graph length on the repaired collars. A good result here says the source overage is repairable locally; it does not certify the whole cell while other unresolved branches remain.

## Claim Ceiling

CELL-05-06 high-slope slab repair diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
