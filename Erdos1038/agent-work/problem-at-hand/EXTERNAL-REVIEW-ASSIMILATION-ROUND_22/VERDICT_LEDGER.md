# Round 22 — Verdict Ledger

External agent: Perplexity Computer
Dispatch: Linear KEN-5 comment `f13c87ca`, 2026-05-28T05:00:00Z
Mode: MODE_1_REVIEW_RESEARCH_OPINE
Bundle return: 2026-05-28 ~06:30 UTC
Assimilation status: **ASSIMILATED** — claim level remains 0; route gains literature foundation and schema target; one new confound (C5) added

---

## Bundle

PC returned a substrate bundle of 5 named files (plus BUNDLE.sha256) under the standing output-rule preference. Bundle landed at `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_22/`. PC's MANIFEST.json declares branch verified (`agent/claude/round-19-assimilation-20260527` @ `33109f0`), no main-fallback. One access anomaly noted: PC returned HTTP 404 when attempting to read `EXTERNAL-REVIEW-ASSIMILATION-ROUND_19/NEXT_LOCAL_GATE.md` — the file exists locally (confirmed). PC correctly pivoted to the playback-row tail for Task C deferred status rather than halting. See METHODOLOGY_NOTES §4 for investigation note.

---

## Per-artifact verdicts

### MANIFEST.json — ACCEPT

PC's self-manifest is internally consistent. Claim level 0 throughout, scope tags (LITERATURE_SCOUT_OUTPUT, SCHEMA_PROPOSAL_ONLY, PUTATIVE_MORPHISMS_ONLY, PUBLIC_SAFE_SCAFFOLD, MODE_1_NO_THEOREM_ADVANCE) accurate, all eight forbidden-claim checks explicit. PC correctly identifies the three deliverables and their scope boundaries. Branch verification confirmed; anchor accepted.

---

### ENDPOINT_LIMIT_LITERATURE_SCOUT.md — ACCEPT

The bibliography is well-organized and substantively correct across all six sections (A–F). Twenty-four entries spanning Deift-Zhou multi-cut outer parametrices, Lubinsky endpoint universality, Kuijlaars-McLaughlin equilibrium theory, Eynard-Orantin recursion, and hyperelliptic period-fiber geometry. Each entry carries a DOI or arXiv identifier, a classification (established theorem / dual I+III), and an explicit relevance note connecting it to the #1038 endpoint-limit gate.

The classification scheme is honest. Three entries carry dual (I)+(III) labels where a result both establishes a theorem and excludes a naive morphism: D3 (Claeys-Wang: universality without CD formula — excludes naive CD morphism for biorthogonal structure), E4 (Eynard et al. 2024: quantization with smooth ramification — excludes topological recursion at hard-wall), and C3 (Kuijlaars-McLaughlin: critical endpoint classification adds an exclusion condition on top of the established theorem). PC did not inflate these to "established morphisms" — they are classified correctly as informing exclusion conditions.

The most directly relevant citation is F2 (Kandaian 2026: hard-soft two-edge Bessel parametrix in a 3×3 RHP). This is the closest structural analog to the #1038 endpoint situation in the literature and is correctly flagged as the primary local model.

No citation is treated as proof of anything. None is drawn from Erdős Atlas or Collider internals. All 24 entries are verifiable against their stated DOIs or arXiv identifiers. ACCEPT unconditionally.

---

### ENDPOINT_LIMIT_KERNEL_MORPHISMS.md — ACCEPT_WITH_NOTES

Five morphisms, each with a falsifiable local test runnable on toy configurations without private receipt access. Two exclusion results. The overall quality is good — better than most Mode 1 morphism proposals. Two notes worth carrying forward, neither blocking acceptance.

**M1 (Lubinsky-Bessel → finite period contribution):** Well-designed. The density-vanishing-exponent fit (T1) on a toy polynomial with a root at x=1 is clean and executable in minutes. The falsification threshold (alpha_fit < −1.5 implies potential period-integral divergence) is correctly stated. ACCEPT without note.

**M2 (Kuijlaars-Mo Cauchy kernel → source kernel residue):** Well-designed. The residue convergence test as ε→0 (T2) is the right mathematical object to probe — the Cauchy kernel on the hyperelliptic surface is exactly what Kuijlaars-Mo identify as the source kernel. Falsification condition (residue diverges) is correct. ACCEPT without note.

**M3 (Bogatyrev period fiber degeneracy → rank drop):** Well-designed, and important for a reason PC explicitly notes: this morphism surfaces C5 (boundary-component rank-drop risk). Round 21's G2.5 result (in flight) probes cond(M_T) directly, which will bound whether C5 materializes at g=24. It is worth noting that Round 20 established cond(M_monomial) ~ 2.1e11 at g=24 *regardless* of component width — the rank-drop effect from M3 is a weaker, more specific mechanism (narrow boundary component specifically). The relationship is: Round 20 showed the monomial basis fails uniformly at g=24; M3 says the Chebyshev-rescaled basis could *also* fail if a component is narrow. These are compatible but distinct failure modes. ACCEPT; C5 tracked as open confound.

**M4 (Deift et al. theta-function parametrix → normal theta-row):** The falsification test is reasonable, but there is a conceptual subtlety in the test's framing that the local agent should note when executing T4. PC describes checking whether the "period row = normal theta row" — but the relevant period matrix is the full symplectic period matrix Ω (both a-cycles and b-cycles), not the bare real a-period matrix M that Round 20 built. The Deift et al. theta function depends on the Abel map, which uses both a- and b-period integrals. T4's numerical test builds only the a-period matrix (integrating holomorphic differentials over support components). This is the same conceptual layer issue as a separate concern in Round 23's FT-02A diagnostic. For the purpose of T4, the test can still give useful information if interpreted correctly — it checks whether the boundary component's a-period row is linearly independent from the normalization row, which is a necessary (not sufficient) condition for the full Deift et al. theta-function to have a normal row for that component. ACCEPT with note to check against symplectic Ω when executing T4 locally.

**M5 (Eichinger-Lukić Robin constant → FIXED_CLOUD_BOUND via Vieta bound):** Well-designed and the most immediately executable. For the toy case f(x)=x^{24}, the true Robin constant is log 24 ≈ 3.178 (the equilibrium measure of the disk `{|x| < 1}^{1/24}` is classical). The Vieta bound is a closed-form upper bound from the leading coefficient of f. T5 is runnable locally in under a minute and gives a direct sanity check on the M5 morphism. ACCEPT without note.

**X1 (CD formula fails for biorthogonal structure, Claeys-Wang):** Correctly identified. PC flags this as a prerequisite check before applying M1–M4 — if the boundary component has non-generic (critical) density vanishing, the source kernel is not the standard CD kernel and the morphisms require revision. This is good scientific hygiene. ACCEPT.

**X2 (Eynard-Orantin quantization requires smooth ramification):** Correctly identified. The hard-wall boundary at x=±1 breaks the smooth-ramification hypothesis of Eynard et al. 2024. Any attempt to derive the endpoint-limit source kernel via topological recursion must address this gap explicitly. ACCEPT.

---

### FIXED_CLOUD_BOUND_CERTIFICATE_SCHEMA_PROPOSAL.json — ACCEPT

The schema is well-structured and covers the right mathematical territory. Seven top-level sections: cloud_descriptor, bound_statement, potential_theory_data, endpoint_limit_data, period_row_compatibility, validation_witnesses, certificate_status. The field names and types are sensible; the `additionalProperties: false` constraints are appropriate for a schema that will be used programmatically.

The 10 cross-field validation rules are the most important design decision here, and they are sound:

- VLD-01 through VLD-03 enforce matrix-shape consistency (gap_period_row_count = component_count − 1; matrix is square of the right dimension).
- VLD-04 through VLD-05 enforce that FULL_INTERVAL_CERTIFIED status requires all six receipts present and all arithmetic modes set to INTERVAL_ARITHMETIC. This is the right gating logic — the schema enforces honest accounting.
- VLD-06 through VLD-07 enforce that ADMISSIBLE admissibility requires FINITE_ROW_COMPATIBLE period contribution and no normalization leak detected. Correct.
- VLD-08 through VLD-09 enforce the condition-number threshold gate (σ_min > 0 and cond < threshold for FULL_INTERVAL_CERTIFIED). Correct.
- VLD-10 enforces that bessel_order_alpha is non-null only when the local_universality_class is BESSEL_HARD_EDGE_ALPHA. Clean.

The enum for `equilibrium_measure_density_type` (GENERIC_SQRT_VANISHING / CRITICAL_HIGHER_ORDER / HARD_EDGE / MIXED) maps correctly to the Kuijlaars-McLaughlin classification that X1 references. The `local_universality_class` enum covers the expected cases (AIRY_SOFT_EDGE, BESSEL_HARD_EDGE, BESSEL_HARD_EDGE_ALPHA, PAINLEVE_CRITICAL, SINE_BULK, UNDETERMINED). The `certificate_status` enum's terminal states include both success (FULL_INTERVAL_CERTIFIED) and three failure modes (REJECTED_NORMALIZATION_LEAK, REJECTED_SINGULAR_MATRIX, REJECTED_CONDITION_THRESHOLD_EXCEEDED). These failure modes correspond exactly to the gates the local agent must pass.

One observation for local adoption: the schema's `arithmetic_mode` enum includes CHEBYSHEV_RESCALED_F64 for the period_matrix section, which correctly tags the Round 20 working basis. A future schema version may want to distinguish CHEBYSHEV_RESCALED_F64 from CHEBYSHEV_RESCALED_INTERVAL when interval arithmetic is added. For now the schema handles this via the outer `certificate_status` field — only FULL_INTERVAL_CERTIFIED requires INTERVAL_ARITHMETIC throughout.

Certificate_status is explicitly set to UNCERTIFIED_SCHEMA_ONLY in the proposal, matching scope. No receipt values are populated. ACCEPT.

---

### PROPOSED_PLAYBACK_ROWS.jsonl — ACCEPT_AND_INCORPORATE_VARIATION

PC proposes two rows: `round_dispatch` (evt-20260528-round-22-dispatch) and `round_complete` (evt-20260528-round-22-complete). Both are claim_level 0, no altitude movement, WORK_PRODUCT_INTENT_ONLY status. The `playback_note` on the round_complete row correctly carries the full route-state snapshot including the new C5 risk and the in-flight R21/R23 context.

Local agent variation: the assimilation packet uses `actor: claude-on-ken-machine` on the assimilation-side completion event (PC's `round_complete` framing is the dispatch-side completion; the assimilation packet's event is the local registration). The C5 confound is annotated as OPEN with falsifiable_via: T3. C1, C2, C4 carry-forward annotations are added per standing pattern.

---

## Aggregate verdict

**Round 22 is a productive Mode 1 scout.** PC delivered what a Mode 1 dispatch should deliver: a solid literature foundation, a clean schema proposal, and five morphisms with mostly-clean falsifiable tests — all without overstating, populating certificates, or inventing receipts.

The 24-citation bibliography covers the right territory. The endpoint-limit regime for real hyperelliptic period integrals was genuinely under-researched in the assimilation record up to this point; the Kandaian 2026 hard-soft two-edge result (F2) and the Lubinsky endpoint universality family (B1–B3) are the most directly useful additions.

The schema proposal is the other significant contribution. Having a well-specified schema for FIXED_CLOUD_BOUND_CERTIFICATE.json means that when the local agent eventually produces private receipts, the target structure is already defined and validated. The 10 cross-field rules will catch consistency errors early.

The one new confound (C5: boundary-component rank-drop risk at g=24) is a real concern that PC deserves credit for surfacing explicitly. The Bogatyrev M3 morphism raises the question of whether the endpoint-limited component in the 25-component cloud is narrow enough to destabilize M_T conditioning, even after the Round 20 Chebyshev-rescaled pivot.

Two methodology notes (T4's bare a-period matrix vs. symplectic Ω framing, and the HTTP 404 access anomaly) are flagged in METHODOLOGY_NOTES without blocking acceptance of any artifact.

**Route state changes:** none to status labels. endpoint_limit_source_kernel remains SUMMIT-LEVEL OPEN; canonical_hyperelliptic_basis remains PRIMARY PARALLEL. The route gains a literature foundation for the endpoint-limit gate and a schema target for the certificate. C5 is added to the open confound list.

No #1038 claim made. Claim level remains 0.
