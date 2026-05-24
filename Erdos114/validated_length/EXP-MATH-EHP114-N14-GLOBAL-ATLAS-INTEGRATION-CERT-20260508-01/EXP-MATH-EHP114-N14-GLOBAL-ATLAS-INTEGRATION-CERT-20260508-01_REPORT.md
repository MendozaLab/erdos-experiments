# EHP114 n=14 Global Atlas Integration Certificate

Experiment: `EXP-MATH-EHP114-N14-GLOBAL-ATLAS-INTEGRATION-CERT-20260508-01`

## Verdict

- Status: `GLOBAL_ATLAS_INTEGRATION_PASS_NOT_FULL_PROOF`
- Discovered cell count: `64` of `64`
- Passing cell count: `64` of `64`
- Blocker count: `0`
- First failed condition: `none`
- Cap consistent across cells: `True`
- Global max length upper: `20.66460519288245`
- Global min margin to cap: `0.008190869737216389`
- Exact length cap: `20.672796062619668`
- Cross-cell unique source IDs: `321`
- Cross-cell source ID collisions: `0`

## Tightest Cells (smallest margin to cap)

| Cell | Length upper | Margin to cap |
|---|---:|---:|
| `CELL-07-03` | `20.66460519288245` | `0.008190869737216389` |
| `CELL-06-03` | `20.65662501573826` | `0.016171046881407136` |
| `CELL-07-05` | `20.50343400901042` | `0.16936205360924816` |
| `CELL-02-04` | `20.469469503184456` | `0.2033265594352116` |
| `CELL-06-04` | `20.316451752723314` | `0.3563443098963539` |

## Loosest Cells (largest margin to cap)

| Cell | Length upper | Margin to cap |
|---|---:|---:|
| `CELL-00-01` | `15.078165474813495` | `5.594630587806172` |
| `CELL-00-06` | `15.119804427110331` | `5.552991635509336` |
| `CELL-00-05` | `15.393007841356184` | `5.279788221263484` |
| `CELL-00-00` | `15.410486279252329` | `5.262309783367339` |
| `CELL-00-03` | `15.415136574295897` | `5.257659488323771` |

## Interpretation

This certificate composes 64 already-emitted local-cell certificate packets into one global integration artifact. It does not recompute interval geometry. The pass condition is bookkeeping integrity: every per-cell packet must carry a passing local certificate, every cell artifact's own RESULTS.json SHA-256 must match its sidecar, every cell must use the same exact length cap, no source experiment ID may be claimed by two cells, and the maximum local length upper across the atlas must remain below the exact length cap.

## Claim Ceiling

Global n=14 atlas integration certificate only. Not a proof of Erdos #114, not a result for n != 14, and not a Tao-bridge or all-degree statement. The all-degree EHP conjecture remains open.

## Next Dependencies

- n=15 boundary slice atlas (Track B in plan yes-wondrous-blum.md), Tao threshold effectivization (Track A), and finite-degree certificate index for n=3..14. This integration cert is necessary but not sufficient for an all-degree EHP proof.
