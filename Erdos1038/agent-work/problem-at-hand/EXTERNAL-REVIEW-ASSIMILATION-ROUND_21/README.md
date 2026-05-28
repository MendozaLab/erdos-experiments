# EXTERNAL-REVIEW-ASSIMILATION-ROUND_21

Local assimilation packet for Perplexity Computer's Round 21 work product on Erdős #1038.

## Files

- `MANIFEST.json` — packet metadata, bundle provenance, per-artifact verdicts, headline results, confound status, route state, scope tags, forbidden-claims confirmation.
- `VERDICT_LEDGER.md` — per-artifact accept/reject reasoning. Aggregate verdict: Round 21 walks back Round 20's working-basis pivot in the strongest possible form — M_T is 8–12 orders of magnitude worse than M_monomial at g=24, not merely comparable.
- `REPRODUCIBILITY_CHECK.md` — two reruns: (1) PC's sweep verbatim (reproduces to within f64 ULP amplification at ~1.4e-7 relative at g=24); (2) sanity follow-up with Round 20's exact branch_points, confirming cond(M_monomial) = 2.119e+11 matches Round 20 G1 anchor and cond(M_T) = 3.119e+20 at g=24 — qualitative finding construction-robust.
- `METHODOLOGY_NOTES.md` — assessment of PC's methods. Notably strong: honest construction caveat with explicit sanity-follow-up recommendation; side-by-side M_monomial recomputation enabling ratio-based decision; Trefethen-2019 structural explanation for the M_T failure mode. Minor flag: frozen dataclass import compatibility in Python 3.9.
- `NEXT_LOCAL_GATE.md` — Round 20 G2.5 carve-out now CLOSED with negative verdict. Working basis UNCERTAIN. Next: Option (c) row-space/cycle-space equivalence (local-only). G3 interval re-implementation remains the gate for any certified seed. Further PC dispatch on the conditioning question is not warranted — the sweep phase is complete.

## Status

`ASSIMILATED` at claim level 0. PC's Round 21 output closes the carve-out that Round 20 explicitly named. The result is a genuine setback for the Chebyshev-rescaled pivot, but PC handled it with the same methodological honesty that made Round 20 the strongest prior return. The construction caveat was surfaced proactively, the side-by-side comparison was the right design, and the Trefethen-2019 structural explanation is correct and non-trivial.

## Confound score-card

- **C1** (tautological QR diagnostic): RESOLVED (Round 20 G2)
- **C2** (unprobed g=24 conditioning): RESOLVED (Round 20 G1)
- **C3** (f64-only sampling): still open — G3 local-only closes it
- **C4** (silent main-fallback): RESOLVED at template level (commit `a4b5bc6`); Round 21 branch-verification preamble passed at `41d8a0f3`
- **G2.5** (M_T own conditioning): **RESOLVED** (Round 21 — negative verdict)

## Route status

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL — **working basis UNCERTAIN** (reverted from Chebyshev-rescaled pivot; neither monomial nor Chebyshev-rescaled viable in f64 at g=24)
- `weighted_qr_basis`: DEMOTED to diagnostic (unchanged)
- `dependent_vieta_image_consumer`: BLOCKED on six absent receipts (unchanged)
- `independent_coefficient_box_theorem`: DANGEROUS (unchanged)
- `endpoint_limit_source_kernel`, `kkt_strict_slack`, `global_reduction`: SUMMIT-LEVEL OPEN (unchanged)

## Bundle source

PC bundle at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_21/`. 5 files: MANIFEST.json, CHEBYSHEV_M_T_CONDITIONING_REPORT.json (53KB, 147 trials), CHEBYSHEV_M_T_RATIONALE.md, chebyshev_m_t_sweep.py, PROPOSED_PLAYBACK_ROWS.jsonl.

## Playback

`evt-20260528-round-21-complete` event to be logged in `Erdos1038/agent-work/EVEREST_ROUTE_PLAYBACK.jsonl` with full route-state snapshot and updated confound status.

## Linear

- Round 21 dispatch comment: `6eec5a97` (2026-05-28, MODE_2_INTENSE_SOLVE_END_TO_END)
