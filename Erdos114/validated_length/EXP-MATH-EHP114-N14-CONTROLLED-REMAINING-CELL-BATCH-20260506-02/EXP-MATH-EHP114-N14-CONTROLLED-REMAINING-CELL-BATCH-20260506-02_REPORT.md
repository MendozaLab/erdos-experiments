# EHP114 n=14 Controlled Remaining-Cell Batch

Experiment: `EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-02`

## Verdict

- Status: `CONTROLLED_REMAINING_CELL_BATCH_BLOCKED`
- Processed cell count: `17`
- Passing cell count: `16`
- Stopped on blocker: `True`
- First failed condition: `CELL-02-07: 2716:0:root/y0/y1 failed because x-wall F intervals do not have strict opposite signs; adaptive y-subdivision left 21 of 26 leaf boxes unresolved`
- Remaining cell count after this run: `37`
- Exact length cap: `20.672796062619668`

| Cell | Status | Final length upper | Margin to cap | Repair targets |
|---|---|---:|---:|---|
| `CELL-00-07` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.627882036790496` | `5.044914025829172` | `3655:0` |
| `CELL-01-00` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.945524819071935` | `0.7272712435477331` | `4063:0` |
| `CELL-01-01` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `18.308490790487795` | `2.3643052721318725` | `3713:1` |
| `CELL-01-02` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.431146086004166` | `3.241649976615502` | `none` |
| `CELL-01-03` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `16.813915919651034` | `3.8588801429686335` | `none` |
| `CELL-01-04` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.53773753688207` | `5.135058525737598` | `none` |
| `CELL-01-05` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.719816162745309` | `4.952979899874359` | `2423:3` |
| `CELL-01-06` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.932766113819524` | `4.740029948800144` | `3655:0,2423:3` |
| `CELL-01-07` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `16.033405318919826` | `4.639390743699842` | `2423:3` |
| `CELL-02-00` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `16.59830584503669` | `4.0744902175829765` | `4063:0` |
| `CELL-02-01` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.604237870652188` | `3.0685581919674796` | `4063:0` |
| `CELL-02-02` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.787525654612743` | `2.8852704080069245` | `3713:1` |
| `CELL-02-03` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.64690787614151` | `3.0258881864781593` | `none` |
| `CELL-02-04` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `20.469469503184456` | `0.2033265594352116` | `none` |
| `CELL-02-05` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `19.71223661397008` | `0.9605594486495868` | `none` |
| `CELL-02-06` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `17.4243076793094` | `3.248488383310267` | `2423:3` |
| `CELL-02-07` | `CONTROLLED_REMAINING_CELL_BATCH_BLOCKED` | `None` | `None` | `none` |

## Interpretation

This batch is a proof-track control gate. It runs the same per-cell source,
repair, residual-chain, and local-packet machinery used by the representative
cells, and it stops on the first blocker. A pass here means only that the
processed local cells received theorem-packet-shaped artifacts under the local
n=14 contract.

## Claim Ceiling

local n=14 cell-batch diagnostic only; not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate
