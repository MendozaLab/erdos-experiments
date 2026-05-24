# EHP114 n=14 Controlled Remaining-Cell Batch

Experiment: `EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-04`

## Verdict

- Status: `CONTROLLED_REMAINING_CELL_BATCH_PASS_NOT_GLOBAL_PROOF`
- Processed cell count: `14`
- Passing cell count: `14`
- Stopped on blocker: `False`
- First failed condition: `none`
- Remaining cell count after this run: `0`
- Exact length cap: `20.672796062619668`

| Cell | Status | Final length upper | Margin to cap | Repair targets |
|---|---|---:|---:|---|
| `CELL-05-07` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `16.357819331297453` | `4.314976731322215` | `none` |
| `CELL-06-00` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.048560443161602` | `3.624235619458066` | `none` |
| `CELL-06-01` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.959595359865958` | `2.71320070275371` | `none` |
| `CELL-06-02` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.893133643859095` | `0.779662418760573` | `3656:0` |
| `CELL-06-03` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `20.65662501573826` | `0.016171046881407136` | `none` |
| `CELL-06-05` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `18.696441846009524` | `1.9763542166101438` | `3713:1` |
| `CELL-06-06` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `18.137586413849377` | `2.535209648770291` | `3713:1,2740:3` |
| `CELL-06-07` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.627330610449697` | `3.045465452169971` | `3713:1` |
| `CELL-07-01` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.86151012994801` | `0.8112859326716588` | `3656:0,2997:0` |
| `CELL-07-02` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `18.358947122360235` | `2.313848940259433` | `3656:0` |
| `CELL-07-03` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `20.66460519288245` | `0.008190869737216389` | `none` |
| `CELL-07-04` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.0475036927646` | `1.6252923698550674` | `none` |
| `CELL-07-05` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `20.50343400901042` | `0.16936205360924816` | `none` |
| `CELL-07-06` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.81811347662989` | `0.854682585989778` | `none` |

## Interpretation

This batch is a proof-track control gate. It runs the same per-cell source,
repair, residual-chain, and local-packet machinery used by the representative
cells, and it stops on the first blocker. A pass here means only that the
processed local cells received theorem-packet-shaped artifacts under the local
n=14 contract.

## Claim Ceiling

local n=14 cell-batch diagnostic only; not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate
