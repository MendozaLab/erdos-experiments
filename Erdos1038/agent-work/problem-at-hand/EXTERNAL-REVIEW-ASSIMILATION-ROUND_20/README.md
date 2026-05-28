# EXTERNAL-REVIEW-ASSIMILATION-ROUND_20

Local assimilation packet for Perplexity Computer's Round 20 work product on Erdős #1038.

## Files

- `MANIFEST.json` — packet metadata, bundle provenance, per-artifact verdicts, headline results, confound status, route state, scope tags, forbidden-claims confirmation.
- `VERDICT_LEDGER.md` — per-artifact accept/reject reasoning. Aggregate verdict.
- `REPRODUCIBILITY_CHECK.md` — local rerun of `genus_growth_sweep.py` (~5s) and `transform_diagnostic.py` (~64s); per-genus / per-block numerical agreement.
- `METHODOLOGY_NOTES.md` — positive assessment of PC's choices (Fraction-arithmetic basis-change construction, iterated Van der Sluis equilibration, structurally + empirically non-tautological diagnostic). Two carve-outs PC named for future work. Confound score-card.
- `NEXT_LOCAL_GATE.md` — G2.5 (Chebyshev-rescaled period matrix conditioning, can be PC or local), G3 (interval re-implementation, local-only), G4 (endpoint-limit source vector, local-only). Round 21 dispatch question (Option α extend Round 20 with G2.5 vs Option β redirect PC to summit-level).

## Status

`ASSIMILATED` at claim level 0. PC's Round 20 output was the strongest #1038 PC return so far — methodologically rigorous (Fraction-arithmetic, Van der Sluis equilibration), reproducible at f64 precision (G1) and to displayed precision (G2 raw cond exact via Fraction), honestly scoped, every dispatch constraint respected including the corrective-dispatch branch checkout (no main-fallback).

## Confound score-card

- **C1** (tautological QR diagnostic): **RESOLVED**
- **C2** (unprobed g=24): **RESOLVED**
- **C3** (f64-only sampling): still open — G3 local-only closes it
- **C4** (silent main-fallback): RESOLVED at template level (Research-Hub commit `a4b5bc6`)

## Route status

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL — **working basis pivoted to Chebyshev-rescaled** within the canonical family (within-family refinement, not a re-demotion)
- `weighted_qr_basis`: DEMOTED to diagnostic (unchanged)
- `dependent_vieta_image_consumer`: BLOCKED on six absent receipts (unchanged)
- `independent_coefficient_box_theorem`: DANGEROUS (unchanged)
- `endpoint_limit_source_kernel`, `kkt_strict_slack`, `global_reduction`: SUMMIT-LEVEL OPEN (unchanged)

## Bundle source

PC bundle at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_20/` (with `BUNDLE.sha256`). 8 files: MANIFEST.json, GENUS_GROWTH_CONDITIONING_REPORT.json (36KB), genus_growth_sweep.py, CONDITIONING_LITERATURE.md, TRANSFORM_DIAGNOSTIC_REPORT.json (478KB), transform_diagnostic.py, TRANSFORM_DIAGNOSTIC_RATIONALE.md, PROPOSED_PLAYBACK_ROWS.jsonl.

## Playback

`evt-20260528-round-20-complete` event logged in `Erdos1038/agent-work/EVEREST_ROUTE_PLAYBACK.jsonl` with full route-state snapshot + confound status after assimilation.

## Linear

- Initial Round 20 dispatch comment: `44d36001` (2026-05-28T03:49:08Z)
- Corrective Round 20 dispatch comment: `a5553f4e` (2026-05-28T04:22:50Z) — after C4 first-attempt freshness defect
- Round 20 receipt comment: posted at assimilation time, links to PR #3 + this packet
