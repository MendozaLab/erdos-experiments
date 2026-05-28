# Round 20 — Verdict Ledger

External agent: Perplexity Computer
Dispatch: Linear KEN-5 (initial `44d36001` 2026-05-28T03:49:08Z; corrective `a5553f4e` 2026-05-28T04:22:50Z after C4 freshness defect)
Bundle return: 2026-05-28 ~21:21 PDT
Assimilation status: **ASSIMILATED** — claim level remains 0; route state pivots within the canonical hyperelliptic family

## Bundle

PC returned a substrate bundle (output rule preference #2) of 8 files, all SHA256-verified byte-for-byte against the Downloads source and against PC's own self-manifest. Bundle landed at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_20/`. PC's MANIFEST.json declares the branch read (`agent/claude/round-19-assimilation-20260527` @ `94606179`, ahead of the anchor `4d5c42b27` per the standing accept-later-head rule) and confirms the corrective dispatch was honored — no main-fallback this round.

## Per-artifact verdicts

### MANIFEST.json — ACCEPT

PC's self-declared manifest is internally consistent: claim level 0 throughout, scope tags (F64_SAMPLED_ONLY, PUBLIC_SAFE_SCAFFOLD, TRIAGE_PROBE_NOT_CERTIFICATE), five non-evidence clauses explicit. SHA-aggregate matches when recomputed.

### GENUS_GROWTH_CONDITIONING_REPORT.json — ACCEPT

147 trials, full sweep grid as specified (`g ∈ {2, 4, 8, 12, 16, 20, 24}`, `cluster_eps ∈ {1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4}`, `jitter ∈ {0, ±1e-4}`). Per-genus max cond(M):

| genus | max cond(M)     | order      |
|-------|-----------------|------------|
| 2     | 2.748           | O(1)       |
| 4     | 1.867e1         | O(10)      |
| 8     | 1.440e3         | O(1e3)     |
| 12    | 1.444e5         | O(1e5)     |
| 16    | 1.581e7         | O(1e7)     |
| 20    | 1.807e9         | O(1e9)     |
| 24    | 2.119e11        | O(1e11)    |

Growth is essentially geometric in genus (~2 orders per Δg = 4), fully consistent with Gautschi 1990 / Pan 2016 Vandermonde-conditioning theory.

Decisive FAIL at g=24 (cond crosses `1e10` threshold uniformly across all 21 (eps, jitter) configurations; first-crossing-genus is 24 in every row). The unrescaled monomial canonical basis is unsuitable at g=24 even before interval-arithmetic considerations.

### genus_growth_sweep.py — ACCEPT

Deterministic, numpy-only. 320 quadrature nodes per gap (Round 19 used 160; PC doubled "for higher genus headroom" — sound choice given amplified ill-conditioning). Endpoint-safe Chebyshev quadrature unchanged from Round 19 (the same scheme cited in Frauendiener-Klein arXiv:1408.2201). Branch-point construction documented and defensively coded against degenerate jitter values. Code reads cleanly cold; reproducibility verified — see REPRODUCIBILITY_CHECK.md.

### CONDITIONING_LITERATURE.md — ACCEPT

Cites Gautschi 1990 (Vandermonde conditioning), Pan 2016 (sharper rates + rescaling motivation), Trefethen 2019 (Vandermonde with Arnoldi + Chebyshev substitute), Frauendiener-Klein 2014 (period-matrix Chebyshev quadrature), Mumford and Cantor (motivational only, explicitly), Higham 2002 (the `cond ≥ 1e10` half-precision-loss threshold rationale). Each citation correctly classified (literature-established vs. motivational vs. threshold-justifying); none claimed to certify #1038.

### TRANSFORM_DIAGNOSTIC_REPORT.json — ACCEPT

147 (genus, eps, jitter) rows, each containing 21 (or g-many) gap-block conditions in both raw and equilibrated form. Aggregate at g=24: max block equilibrated cond = 1.176e7; max block raw cond = 4.631e118. Verdict `DEFENSIVE_PIVOT_TO_CHEBYSHEV_RESCALED_IS_ACTIONABLE_UNDER_EQUILIBRATION` is correctly licensed by the data (1.176e7 < 1e10 threshold; 9 orders below weighted-QR demote cond 4.7e16).

The four-orders-of-magnitude separation between cond(M) (~2.1e11) and equilibrated cond(P) (~1.18e7) at g=24 is the empirical non-tautology proof.

### transform_diagnostic.py — ACCEPT_WITH_NOTES

Acceptable; the methodological choices PC made here are notably strong (see METHODOLOGY_NOTES.md §1). Two design choices worth flagging in writing — neither blocks acceptance:

1. **Fraction-arithmetic basis-change construction.** P is built in exact rational (`Fraction`) arithmetic to avoid f64 catastrophic cancellation in the binomial expansion of `(mid + half·t)^n` at small `half` and high `n`. Down-cast to f64 only at the end before SVD. This is the right call — building in f64 directly would have produced `inf` in the SVD at g=24 / eps=1e-4 (as PC notes in the rationale). The cost is ~64 seconds on this machine for the full sweep, which is fine for triage.
2. **Iterated Van der Sluis equilibration (6 iterations).** Standard 1-pass equilibration converges to within `√min(m,n)` of the optimal diagonal-scaled cond; iteration converges fully. PC's choice of 6 iterations is well within the convergence regime for this problem size. Equilibrated cond is the scale-invariant diagnostic per Van der Sluis 1969 — the right number to report when the raw cond is dominated by trivial column-scale gaps (here ~118 orders of magnitude).

### TRANSFORM_DIAGNOSTIC_RATIONALE.md — ACCEPT

Excellent self-documentation. Explicitly addresses (a) why Option (b) over (a)/(c), (b) why equilibrated cond is the correct number, (c) what the verdict does NOT say (no quadrature error certificate, no certification of Chebyshev-rescaled period matrix's own conditioning, no relaxation of the six-receipt blocker, no altitude movement). The defer-with-reasoning for Options (a) and (c) is sound: (a) would re-exhibit weighted-QR's already-encoded failure; (c) requires private receipt data.

One implicit claim worth surfacing (not a flaw, just a follow-up for the local gate): the diagnostic establishes that the **basis-change** from monomial to Chebyshev-rescaled is well-conditioned, NOT that the **Chebyshev-rescaled period matrix M_T itself** has good condition. M_T conditioning is a separate G2-follow-up (or part of G3) that the local agent must compute. PC explicitly names this carve-out in the rationale.

### PROPOSED_PLAYBACK_ROWS.jsonl — ACCEPT_AND_INCORPORATE_VARIATION

PC proposes one `round_substrate_return` event with full evidence + blockers + next_actions + route-state-aware playback_note. Local agent variation: this assimilation packet's playback row instead uses `event_type: round_complete` (PC's `round_substrate_return` is the dispatch-side framing; `round_complete` is the assimilation-side completion event), `actor: claude-on-ken-machine`, and incorporates C1 + C2 RESOLVED status, C3 STILL_OPEN, C4 RESOLVED_BY_TEMPLATE_UPDATE annotations per the standing real-time-with-route-snapshot pattern.

## Aggregate verdict

**Round 20 is the strongest PC return so far for #1038.** PC closed both of Round 19's methodology gaps with rigorous, reproducible, well-grounded work:

- **C1 (tautological QR diagnostic) is fully resolved.** PC didn't just pick a different diagnostic — they proved non-tautology empirically (4-order-of-magnitude separation between cond(M) and cond(P) at g=24) and structurally (P depends only on gap geometry, not on M's cycle integrals; this is verifiable by reading the code).
- **C2 (unprobed genus-24 conditioning) is fully resolved.** The sweep is complete through g=24 with all 21 (eps, jitter) configurations covered. The FAIL at g=24 is decisive and quantitatively consistent with Vandermonde-conditioning theory.
- **C3 (f64-only sampling) remains open** — both sweeps are still F64_SAMPLED_ONLY. This is expected; G3 (interval-arithmetic re-implementation) is local-only and out of PC's scope.

PC respected every dispatch constraint: explicit branch checkout (no main-fallback this round), no push, no PR, no claim promotion, no receipt invention, no chaining into Round 21. Honest scope statements throughout.

**The route status changes:**

- `canonical_hyperelliptic_basis` stays PRIMARY PARALLEL — but the working basis pivots from unrescaled monomial to Chebyshev-rescaled within the canonical family.
- This is a meaningful refinement, not an altitude move. The route is now narrower (one working basis specified) but better-grounded.
- Six receipts still absent. Six-receipt blocker on dependent-Vieta consumer unchanged.
- Summit-level open gates (endpoint-limit kernel, KKT/strict-slack, global reduction) unchanged.

No #1038 claim made. Claim level remains 0.
