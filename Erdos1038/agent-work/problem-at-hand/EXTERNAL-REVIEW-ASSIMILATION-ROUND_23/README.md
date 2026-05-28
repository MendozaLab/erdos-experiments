# EXTERNAL-REVIEW-ASSIMILATION-ROUND_23

Local assimilation of Perplexity Computer's Round 23 return for Erdős #1038.

**PC dispatch:** Linear KEN-5 comment `9638bdc3`, 2026-05-28
**Mode:** MODE_1_REVIEW_RESEARCH_OPINE (global reduction scout — theorem-family morphism hunt)
**Claim level:** 0 throughout

---

## What PC did

PC surveyed six theorem families for structural alignment with #1038's φ-k-gap surface:

- Three received **PUTATIVE_MORPHISM** classification with explicit component-to-component maps and falsifiable local tests requiring no private payload.
- Three received **NEGATIVE_EXCLUSION** with specific stated reasons (not loose analogies).

The strongest putative morphism maps the 24-row gap-period matrix to the period matrix of a genus-24 hyperelliptic curve, with period legitimacy equivalent to the Riemann bilinear relations. The second maps the sublevel-set measure bound to the NPS lower bound via bounded gap-component count. The third identifies gap-period matrix entries with KKT Jacobian components at the extremal configuration.

PC also noted that Morphisms A and C may compose: if both hold, the extremal polynomial's KKT stationarity conditions translate into the Siegel upper half-space constraint on the period matrix. This is the most intellectually interesting direction to formalize next.

---

## What this assimilation adds

**Methodology critique.** The local agent ran a toy g=2 check and found that FT-02A as written is conceptually mismatched: it tests the real a-period matrix M_a for symmetry, but symmetry is a property of the full symplectic period matrix Omega = M_a^{-1} M_b. Running the test naively would produce a false FAIL. The corrected FT-02A protocol (computing both M_a and M_b, then forming Omega) is documented in METHODOLOGY_NOTES.md §2 and NEXT_LOCAL_GATE.md Action 3.

The toy g=2 sanity check confirmed that the route infrastructure is sound: M_a is well-conditioned (cond = 2.625, det = 13.2, non-degenerate) and correctly asymmetric — the asymmetry is expected mathematical behavior.

---

## Files in this packet

| File | Contents |
|------|----------|
| `MANIFEST.json` | Bundle provenance, per-artifact verdicts, route state after assimilation |
| `VERDICT_LEDGER.md` | Per-artifact verdicts with narrative assessment |
| `METHODOLOGY_NOTES.md` | Detailed methodology assessment: exclusion quality, FT-02A gap and corrected protocol, FT-05A/FT-06A quality, composability note |
| `NEXT_LOCAL_GATE.md` | Recommended local actions: FT-05A (inline, cheap), FT-06A (inline, cheap), corrected FT-02A (substantial), A+C composability (future PC round) |
| `README.md` | This orientation file |

---

## Route state one-liner

Global reduction: SUMMIT-LEVEL OPEN — three putative morphisms named with explicit maps (hyperelliptic Jacobian, NPS theorem, KKT optimality), none locally verified. FT-05A and FT-06A are runnable inline. Corrected FT-02A requires complex b-period contour integration. No #1038 altitude movement.

---

## Prior round continuity

| Round | Mode | Key outcome |
|-------|------|------------|
| 19 | MODE_2 | Canonical hyperelliptic basis established; conditioning blocker at g=24 surfaced (C1, C2) |
| 20 | MODE_2 | C1 (tautological QR) and C2 (unprobed g=24) resolved; working basis pivots to Chebyshev-rescaled |
| 21/22 | MODE_2 | G2.5 Chebyshev period matrix conditioning resolved |
| **23** | **MODE_1** | **Global reduction scout: 6 families, 3 putative morphisms, 3 exclusions; FT-02A methodology gap surfaced** |
