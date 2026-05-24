# EHP114 n=14 Per-Cell Source Generation Smoke

Experiment: `EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-01`

## Verdict

- Status: `PER_CELL_SOURCE_BLOCKED_SOURCE_ARTIFACT_MISSING`
- Smoke subcell: `CELL-00-00`
- Upstream binaries audited: `4`
- Hardcoded blocker count: `0`
- Source contract status: `SOURCE_ARTIFACT_MISSING`
- First missing source artifact: `../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json`
- First failed condition: `required first slab source artifact is missing for CELL-00-00: ../../Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json`

## Meaning

This L33C smoke gate asks whether the first missing cell can start the real source chain. It does not generate a certificate and does not rename hard-cell data as cell-local evidence. If the `CELL-00-00` slab source is absent or a source declares the wrong subcell, the proof pipeline stops there.

## Claim Ceiling

Per-cell source-generation smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
