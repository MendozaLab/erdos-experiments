# EHP114 n=14 Controlled Remaining-Cell Batch

Experiment: `EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-01`

## Verdict

- Status: `CONTROLLED_REMAINING_CELL_BATCH_BLOCKED`
- Processed cell count: `5`
- Passing cell count: `4`
- Stopped on blocker: `True`
- First failed condition: `CELL-00-07: source_over_cap_but_no_repair_targets`
- Remaining cell count after this run: `53`
- Exact length cap: `20.672796062619668`

| Cell | Status | Final length upper | Margin to cap | Repair targets |
|---|---|---:|---:|---|
| `CELL-00-03` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.415136574295897` | `5.257659488323771` | `none` |
| `CELL-00-04` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.557494057270386` | `5.115302005349282` | `none` |
| `CELL-00-05` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.393007841356184` | `5.279788221263484` | `none` |
| `CELL-00-06` | `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF` | `15.119804427110331` | `5.552991635509336` | `none` |
| `CELL-00-07` | `CONTROLLED_REMAINING_CELL_BATCH_BLOCKED` | `None` | `None` | `none` |

## Interpretation

This batch is a proof-track control gate. It runs the same per-cell source,
repair, residual-chain, and local-packet machinery used by the representative
cells, and it stops on the first blocker. A pass here means only that the
processed local cells received theorem-packet-shaped artifacts under the local
n=14 contract.

## Claim Ceiling

local n=14 cell-batch diagnostic only; not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate
