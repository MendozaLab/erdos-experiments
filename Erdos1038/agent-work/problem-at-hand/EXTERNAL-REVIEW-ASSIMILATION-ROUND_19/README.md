# EXTERNAL-REVIEW-ASSIMILATION-ROUND_19

Local assimilation packet for Perplexity Computer's Round 19 work product on Erdős #1038.

## Files

- `MANIFEST.json` — packet metadata, bundle provenance, per-artifact verdicts, status, scope tags, forbidden-claims confirmation.
- `VERDICT_LEDGER.md` — per-artifact accept/reject reasoning. Aggregate verdict.
- `REPRODUCIBILITY_CHECK.md` — local rerun of `round19_canonical_falsifier_sweep.py`; structural identity + f64 last-bit numerical agreement.
- `METHODOLOGY_NOTES.md` — two flags surfaced during assimilation: (1) QR-transform diagnostic is tautological with cond(M), (2) genus-24 conditioning unprobed.
- `NEXT_LOCAL_GATE.md` — four-gate ladder (G1–G4) for promoting canonical-basis route from "starting proposal" to "interval-certified seed." Round 20 dispatch question addressed.

## Status

`ASSIMILATED` at claim level 0. PC's Round 19 output was honest, well-scoped, and respected every constraint in the dispatch. Two methodology gaps named for the next gate to close.

## Bundle source

PC bundle lives at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_19/` (with `BUNDLE.sha256` for provenance). The bundle itself is referenced by SHA from `MANIFEST.json` in this packet — no copies of PC's files are duplicated here.

## Playback

`evt-20260528-round-19-complete` event logged in `Erdos1038/agent-work/EVEREST_ROUTE_PLAYBACK.jsonl`.
