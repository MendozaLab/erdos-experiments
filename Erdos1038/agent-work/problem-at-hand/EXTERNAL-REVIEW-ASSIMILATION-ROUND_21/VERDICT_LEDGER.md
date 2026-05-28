# Round 21 — Verdict Ledger

External agent: Perplexity Computer
Dispatch: Linear KEN-5 comment `6eec5a97` (2026-05-28, MODE_2_INTENSE_SOLVE_END_TO_END)
Bundle return: 2026-05-28 ~05:00 UTC
Assimilation status: **ASSIMILATED** — claim level remains 0; Round 20's working-basis pivot to Chebyshev-rescaled is walked back

## Bundle

PC returned a substrate bundle of 5 files, all SHA256-verified against the downloaded source and against PC's own self-manifest. Bundle landed at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_21/`. PC's MANIFEST.json declares the branch read (`agent/claude/round-19-assimilation-20260527` @ `41d8a0f3`, later than dispatch anchor `33109f0`; accepted per the standing later-head rule) and confirms no main-fallback was used this round. Branch verification preamble explicitly passed.

## Per-artifact verdicts

### MANIFEST.json — ACCEPT

PC's self-declared manifest is internally consistent: claim level 0 throughout, scope tags correct (F64_SAMPLED_ONLY, PUBLIC_SAFE_SCAFFOLD, TRIAGE_PROBE_NOT_CERTIFICATE), ten non-evidence clauses explicit. Construction caveats are honestly surfaced in a dedicated field — PC names the public-safe-scaffold limitation and recommends the sanity follow-up before the local agent needs to ask. SHA-aggregate matches when recomputed.

### CHEBYSHEV_M_T_CONDITIONING_REPORT.json — ACCEPT

147 trials, full sweep grid as specified (`g ∈ {2, 4, 8, 12, 16, 20, 24}`, `cluster_eps ∈ {1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4}`, `jitter ∈ {0, ±1e-4}`). Per-genus maxima, side-by-side for both M_T and M_monomial on the same branch-point construction:

| genus | max cond(M_T)  | max cond(M_monomial) | ratio M_T/M_mono |
|-------|----------------|----------------------|-----------------|
| 2     | 4.03e+01       | 4.25e+00             | 9.5             |
| 4     | 1.30e+04       | 2.82e+01             | 5.0e+02         |
| 8     | 3.61e+09       | 1.03e+03             | 3.5e+06         |
| 12    | 1.03e+17       | 3.75e+04             | 2.8e+12         |
| 16    | 4.46e+18       | 1.35e+06             | 3.3e+12         |
| 20    | 3.10e+19       | 4.81e+07             | 6.4e+11         |
| 24    | **2.28e+19**   | **1.70e+09**         | **1.3e+10**     |

The M_T growth pattern is super-geometric and essentially singular well before g=24. cond(M_T) first crosses the Higham 1e10 threshold at g=8 — sixteen genus steps earlier than M_monomial's first crossing (g=24 in Round 20 G1). The ratio at intermediate genus (g=12-16) peaks near 3e12, worse than at g=24 where f64 saturation effects compress the spread. The decisive word is in the side-by-side: cond(M_T) / cond(M_monomial) is between 10^9 and 10^12 across every genus ≥ 8. The Chebyshev-rescaled period matrix is not better-conditioned than the monomial one — it is orders of magnitude worse.

Decision: `RELABELING_ONLY_PIVOT_TO_ROW_SPACE_CYCLE_SPACE_EQUIVALENCE`, strengthened to "M_T is empirically worse than M_monomial, not merely comparable."

### chebyshev_m_t_sweep.py — ACCEPT_WITH_NOTES

Acceptable; the core logic is sound. One design flag worth documenting, not a blocker:

The script uses `@dataclass(frozen=True)` for the `SweepResult` container. When the script is used as an importable module (rather than run as `__main__`), Python 3.9's `dataclasses` module can fail to reconstruct the frozen dataclass from a dynamically imported module scope — this is a known Python 3.9 compatibility issue with `frozen=True` and reload/import cycles. The local sanity rerun worked around this by extracting the `branch_points()` function manually rather than importing the module whole (see REPRODUCIBILITY_CHECK.md for the workaround). Recommendation: either drop `frozen=True` or provide a function-only entry point (e.g., a `run_sweep()` function that returns a plain dict) for downstream consumers that need to import the module programmatically. As a standalone script run via `python3 chebyshev_m_t_sweep.py`, it works without modification.

The quadrature implementation, the side-by-side M_monomial recomputation, and the decision-rule execution are all correct.

### CHEBYSHEV_M_T_RATIONALE.md — ACCEPT

Excellent self-documentation, on par with Round 20's `TRANSFORM_DIAGNOSTIC_RATIONALE.md`. PC's structural explanation for why M_T fails where the basis-change P succeeds is the strongest piece of writing in this bundle — it closes a conceptual gap that Round 20 could not address without the direct computation.

The Trefethen 2019 framing is genuinely insightful. The key observation is that Trefethen's Vandermonde-with-Arnoldi result applies to *evaluation* matrices `V[i,j-1] = T_{j-1}(x_i)` over a single interval, where the Chebyshev structural advantage (orthogonality, controlled growth) holds globally across all nodes. M_T here is structurally different: per-gap normalization breaks global Chebyshev orthogonality across rows, and integration against `1/y(x)` introduces a smooth per-gap Jacobian whose Chebyshev expansion decays geometrically. That geometric decay is the source of ill-conditioning — each row of M_T becomes nearly rank-1 in its Chebyshev expansion, pushing small singular values toward zero as genus grows. This is not a hand-wave; it is correct mathematical reasoning that would survive scrutiny from a numerical analyst.

PC's statement of what the verdict does NOT establish is equally strong: no #1038 claim, no altitude movement, no implementation of Option (c), no certification that the construction isn't a pathological outlier, and an explicit recommendation to rerun with the Round 20 byte-identical branch_points once Research-Hub is accessible. That recommendation is exactly what the local sanity follow-up executed — see REPRODUCIBILITY_CHECK.md.

### PROPOSED_PLAYBACK_ROWS.jsonl — ACCEPT_AND_INCORPORATE_VARIATION

PC proposes one `round_substrate_return` event with full evidence field, blockers (six receipts, G3, G4, construction caveat), next_actions (sanity follow-up, Option (c) as next local gate, route status revert to UNCERTAIN), and a playback_note that correctly distinguishes the Trefethen 2019 evaluation-matrix context from the period-matrix failure mode here. All blockers and next_actions are correctly captured.

Local agent variation, identical to Round 20 practice: the assimilation-side event uses `event_type: round_complete`, `actor: claude-on-ken-machine`, and incorporates the G2.5 RESOLVED annotation alongside C1/C2 RESOLVED, C3 STILL_OPEN, C4 RESOLVED status per the standing real-time-with-route-snapshot pattern. Route-state change (working basis UNCERTAIN) and construction caveat are incorporated in the playback_note.

## Aggregate verdict

Round 21 walks back Round 20's working-basis pivot in the strongest possible form. Where Round 20 ended with cautious optimism — "the basis-change to Chebyshev-rescaled is well-conditioned; M_T's own conditioning is a named carve-out for follow-up" — Round 21 closes that carve-out with a decisive negative result. M_T is not just comparable to the monomial basis at g=24; it is 8–12 orders of magnitude worse, depending on construction. The first crossing of the Higham 1e10 threshold arrives at g=8 for M_T, versus g=24 for M_monomial.

PC's transparency about the construction caveat is exemplary. Rather than presenting cond(M_T) in isolation against Round 20's published cond(M_monomial) anchor, PC recomputed cond(M_monomial) on the same branch-point construction so the ratio is controlled. That's the right call: the cond(M_monomial) anchor from Round 20 (2.12e11) vs. this sweep's same-construction value (1.70e9) differ by ~100×, consistent with the public-safe-scaffold being a representative but not byte-identical reconstruction. The ratio cond(M_T)/cond(M_monomial) is 1.34e10 in this sweep and 1.47e9 in the sanity rerun with Round 20's exact branch_points — both are deeply in the "M_T fails" regime, regardless of which branch_points construction is used.

The route state changes are real:

- `canonical_hyperelliptic_basis` stays PRIMARY PARALLEL — but the working basis reverts from "Chebyshev-rescaled" to UNCERTAIN. This is not a re-demotion of the route itself; it is an honest acknowledgment that Round 20's within-family refinement was premature.
- The next working-basis candidate is Option (c) row-space/cycle-space equivalence, which is local-only because it requires private receipt data to fix alternate gap-row choices.
- Six receipts still absent. Six-receipt blocker on the dependent-Vieta consumer unchanged.
- Summit-level open gates (endpoint-limit kernel, KKT/strict-slack, global reduction) unchanged.

No #1038 claim made. Claim level remains 0. The path forward for altitude movement is the local-Codex receipt sprint, not additional PC rounds on the working-basis question.
