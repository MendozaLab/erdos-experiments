# Round 23 — Verdict Ledger

External agent: Perplexity Computer
Dispatch: Linear KEN-5 comment `9638bdc3`, 2026-05-28
Mode: MODE_1_REVIEW_RESEARCH_OPINE
Assimilation status: **ASSIMILATED** — claim level remains 0; global_reduction route advances from 'named open thread' to 'three candidate morphisms with explicit maps'

---

## Bundle

PC returned a substrate bundle of 5 files under `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_23/`. This was a Mode 1 (review / research / opine) round, not an end-to-end solve, so no computational sweep scripts or result JSON files are present — the deliverables are a literature survey, a structured theorem-family inventory, and three morphism proposals with falsifiable tests. That is appropriate and exactly what was requested.

PC declared `claim_ceiling: 0` throughout, no altitude movement, no period legitimacy claim, no private receipts. The round is honest.

---

## Per-Artifact Verdicts

### MANIFEST.json — ACCEPT

PC's self-declared manifest is internally consistent. Claim ceiling 0. Forbidden-claims list is complete and accurate (twelve items, covering the expected exclusions: #1038 solved, altitude upgrade, period legitimacy claimed, KKT/global reduction closed, any morphism verified, main-fallback used, ErdosAtlas/Collider treated as external literature). SHA hashes are listed as PENDING_LOCAL_HASH, which is correct for a Mode 1 bundle (no computational outputs to hash-verify beyond the prose documents). Morphism summary section correctly classifies all six families.

### GLOBAL_REDUCTION_SCOUT.md — ACCEPT

A well-structured annotated bibliography covering six theorem families with clear, principled methodology. Each family gets: overview, structural assessment, key references with DOIs, and a classification. The exclusion arguments are specific rather than vague:

- TF-01 (modular forms): exclusion grounded in the correct observation that the gap structure in #1038 is topological (connected components of a real sublevel set), not spectral (eigenvalue gaps in a Hecke module). The "topological obstruction" paragraph naming the absence of a natural Galois action compatible with Hecke commutation is exactly the right diagnosis.

- TF-03 (Weil bounds): exclusion grounded in the correct observation that Weil/Deligne applies to characteristic-q fields, while #1038 lives entirely over ℝ. PC notes the structural analogy (both are "saving" estimates) but correctly declines to treat analogy as morphism.

- TF-04 (Sato-Tate): exclusion correctly identifies that Sato-Tate is an ensemble-average statement (over primes or over families of curves), whereas #1038's gap surface is a single fixed object. The note on indirect connection via TF-02 (the Sato-Tate group of the associated genus-24 Jacobian would be a structural invariant) is a legitimate carve-out — PC doesn't over-claim it but doesn't throw it away either.

**Praise: the exclusion documentation is unusually careful.** PC did not fall into the common trap of calling these theorem families "related" and leaving them in limbo. NEGATIVE_EXCLUSION with a specific stated reason is the right disposition.

### THEOREM_FAMILIES_INVENTORY.json — ACCEPT

Well-formed structured catalog. Each family has a `classification`, `morphism_status`, and (for putative morphisms) `morphism_strength` and a `falsifiable_test` block with `test_id`, `description`, `requires_private_receipts`, and `test_type`. This structure matches what a local agent needs to triage and execute.

The `requires_private_receipts: false` flag on FT-02A, FT-05A, and FT-06A is correct — all three tests work with synthetic toy configurations. One of the three tests has a methodological issue (FT-02A — see METHODOLOGY_NOTES.md for details), but the structural intent of the test is sound; only the specific formula for what to verify needs correction.

### PUTATIVE_MORPHISMS_TO_1038.md — ACCEPT_WITH_METHODOLOGY_FLAG

**The overall quality of this document is high.** Each morphism is presented with: the theorem-family context, an explicit component-to-component map, a statement of what the morphism would reduce #1038 to, and a falsifiable test. The ranking A > B > C is defensible.

**Morphism A (hyperelliptic Jacobian, strongest).** The structural alignment is genuine. Mapping gap-period matrix rows to Abelian integrals on the genus-24 hyperelliptic curve `y² = ∏(x-a_j)(x-b_j)` is the right identification. The reduction question ("does M satisfy the Jacobian constraints?") is correctly framed. The composability observation with Morphism C (A + C => Siegel upper half-space constraint) at the end of the document is highlighted as a research direction, not a claim.

**Flag on FT-02A.** The test as written asks to verify `Omega_{12} = Omega_{21}` (symmetry) and `det(Im(Omega)) > 0` (positive-definiteness) for the matrix `Omega_{ij} = integral_{G_j} x^{i-1} / sqrt(|prod_k (x-a_k)(x-b_k)|) dx`. The issue is that this integral computes the **real a-period matrix** `M_a`, not the full symplectic period matrix `Omega`. The real a-period matrix `M_a` is generally **not** symmetric (and has no imaginary part). Running FT-02A as written would produce a false FAIL on symmetry, incorrectly suggesting morphism A is invalid.

A local toy check at g=2 (curve `y² = (x+1)(x+0.5)(x-0.2)(x-0.5)(x-0.7)(x-1)`) confirms:
- Real a-period matrix M_a = [[1.92, -1.35], [4.53, 3.70]]; cond = 2.625; det = 13.2
- M_a[0,1] = -1.35 ≠ M_a[1,0] = 4.53 — correctly **asymmetric**
- This is expected and mathematically correct behavior

The corrected FT-02A requires computing both M_a (real a-period integrals) and M_b (b-period integrals, which involve complex contours), then forming `Omega = M_a^{-1} M_b`, and checking symmetry + positive-definiteness of *that* matrix. See METHODOLOGY_NOTES.md §2 for the corrected protocol.

**Morphism B (NPS theorem, second strongest).** The explicit map is well-constructed. The identification of `nu(T, h) ≤ 2 × gap_count ≤ 50` is correct given the cloud structure's 25-component maximum. The derivation of the degree-independent lower bound is sound. The honest caveat — that the component-count bound applies to a specific candidate configuration, not to all monic real-rooted polynomials — is correctly included. FT-05A is a sensible test: component-count sweep for degrees {10, 20, 50, 100} is cheap to run inline and directly falsifies the bounded-component-count assumption.

**Morphism C (KKT optimality, third).** The co-area formula derivation connecting `d/d(a_i) m` to the gap-period matrix entries is explicit and structurally sound. The reduction to "extremal polynomial's KKT conditions equivalent to period Jacobian being zero" is correctly framed as a research observation, not a proof. FT-06A (finite-difference Jacobian vs. analytic formula at degree-10) is a clean, cheap, inline test.

**Composability note.** The observation that Morphisms A and C may compose — if M = Omega (A) and M = KKT Jacobian (C) both hold, then KKT stationarity ⇔ Omega ∈ H_{24} — is the most intellectually interesting product of this round. PC presents it as a research-level question, not a claim. That's the right disposition. If formalized, this would be a meaningful global reduction.

### PROPOSED_PLAYBACK_ROWS.jsonl — ACCEPT_AND_INCORPORATE_VARIATION

PC proposes two rows: `round_dispatch` (evt-20260528-round-23-dispatch) and `round_complete` (evt-20260528-round-23-complete). Both are internally consistent, claim-level 0, with accurate blocker lists and next-action items. Local agent variation: this assimilation packet's playback row uses `actor: claude-on-ken-machine`, incorporates the FT-02A methodology gap as a new confound entry, and updates the `global_reduction` route status from "SUMMIT-LEVEL OPEN (unnamed)" to "SUMMIT-LEVEL OPEN (three candidate morphisms with explicit maps)" per the standing real-time-with-route-snapshot pattern.

---

## Aggregate Verdict

Round 23 is a productive Mode 1 literature scout. PC delivered exactly what was asked: a careful survey of six theorem families, three morphism proposals with explicit component-to-component maps and falsifiable tests requiring no private payload, and three principled exclusions with specific stated reasons.

The strongest result — the hyperelliptic Jacobian / period matrix alignment (Morphism A) — is a structurally serious research direction. The period-legitimacy problem (does M lie in the Siegel upper half-space H_{24}?) is directly connected to the existing route's open question, and the morphism would be meaningful if it survives local testing.

The methodological gap in FT-02A (testing M_a for symmetry, which is the wrong matrix) would have produced a false FAIL if executed naively. The local toy check surfaced this at minimal cost. The morphism A's structural claim is sound; only the test protocol needed correction.

No #1038 claim is made. No altitude movement. Claim level remains 0. The global_reduction route is now richer with named, structured morphisms, but remains SUMMIT-LEVEL OPEN until tests pass.
