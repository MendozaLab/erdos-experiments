# Erdős Problem #30 — $1,000 Prize Eligibility Assessment

**Date:** 2026-04-15
**Assessor:** MendozaLab agent (honest evaluation per Prime Directive)
**Score (current):** A0 / B5 / C1 (93% confidence)

---

## What the Prize Requires

Erdős Problem #30 ($1,000, posted by Erdős): Prove or disprove that for every ε > 0,

    h(N) = √N + O_ε(N^ε)

where h(N) is the maximum size of a Sidon set in {1, …, N}.

The prize is for **resolving the conjecture** — a proof or disproof. Partial results, formalizations, and reformulations do not qualify.

## Current SOTA (Not Ours)

- **Upper bound:** h(N) ≤ √N + 0.98183·N^{1/4} + O(1) (Carter–Hunter–O'Bryant, 2025)
- **Lower bound:** h(N) ≥ √N + ½ + O(N^{-1/4}) (Singer 1938, via perfect difference sets)
- **The gap:** The upper bound has error O(N^{1/4}); the conjecture asks for O_ε(N^ε). Closing this gap has been open for 85 years.

## What We Have — Honest Inventory

### Track 1: Classical Formalization (A0/B5/C1)

8 Lean 4 files, 29 theorems, 5 cited axioms, 0 sorry, lake build PASS. Formalizes:
- Erdős–Turán upper bound (1941)
- Lindström parametric/weak/full bounds (1969)
- BFR core inequality (2023)
- Singer construction lower bound (1938)
- Additive energy + dispersion machinery

**Prize relevance: NONE.** This formalizes *known* results. No new mathematics. No improvement to any constant. This is correctly scored A0.

### Track 2: Spectral Sidon / Choreography Note

Finite computational verification that D_n = {2sin(πm/n)} is Sidon in ℝ for n=3..29, with collision at n=30.

**Prize relevance: NONE.** This is a finite verification over 27 cases. It does not establish anything asymptotic and has no bearing on h(N) for integer Sidon sets in {1,…,N}. The connection to choreographies is interesting but does not advance the conjecture.

### Track 3: Dispersion Paper Draft (sidon_dispersion_paper_draft.tex)

Claims four theorems:
1. **Thm 1 (zero unordered dispersion):** Correct and elementary — follows directly from Sidon injectivity.
2. **Thm 2 (ordered dispersion = (log 2)²/k):** Correct — clean calculation.
3. **Thm 3 (PPV upper bound: h(N) ≤ √N + C_ε·N^{1/4}/√(log N)):** UNVERIFIED — see analysis below.
4. **Thm 4 (Erdős equivalence via dispersion):** ESSENTIALLY TAUTOLOGICAL — see analysis below.

#### Critical Analysis of Theorem 3

The paper claims to improve the Lindström bound by a factor of √(log N) using the PPV finite-blocklength converse. Issues:

1. **Self-referential dispersion.** The author acknowledges this: "But this is not quite right." The dispersion V = (log 2)²/k depends on the Sidon set size k, which is the quantity being bounded. The fixed-point argument on lines 278–303 is heuristic, not rigorous.

2. **PPV framework applicability.** The PPV converse (their Theorem 54) is designed for discrete memoryless channels with n independent uses at blocklength n. Here n = ⌈log₂(2N)⌉ and there is effectively one channel use (one sum). The second-order term √(nV)·Q⁻¹(ε) is not designed for this regime.

3. **Even if correct, does not resolve the conjecture.** The claimed bound is h(N) ≤ √N + O(N^{1/4}/√(log N)). The conjecture requires O_ε(N^ε). The gap between N^{1/4}/√(log N) and N^ε remains enormous. This would be a modest improvement to the upper bound constant, not a resolution.

#### Critical Analysis of Theorem 4

Proposition 4 reveals that V_ord depends ONLY on |A|, not on the structure of A. This means the dispersion reformulation reduces to: "Erdős conjecture ⟺ there exist Sidon sets with |A| > √N + f(N) for appropriate f." This is tautologically equivalent to the original conjecture, not a new characterization. The author writes "This might seem circular, but it is not" — but the argument for non-circularity is unconvincing.

## Gap Analysis: What Would Be Needed

To claim the $1,000 prize, one of the following would be required:

### Path A — Resolve the Conjecture (Proof)
- Close the 85-year gap between O(N^{1/4}) and O_ε(N^ε)
- This requires a fundamentally new technique; no existing approach (Cauchy-Schwarz, PPV, or Fourier L⁴) is known to bridge this gap
- Estimated difficulty: **Extreme** (85 years open, top-tier combinatorialists have tried)

### Path B — Resolve the Conjecture (Disproof)
- Construct a family of Sidon sets achieving h(N) = √N + Ω(N^δ) for some fixed δ > 0, OR
- Prove that the O(N^{1/4}) error is tight (i.e., the conjecture is false)
- Both would be major results in additive combinatorics

### Path C — Sufficient Partial Progress (Prize Committee Discretion)
- Some Erdős prizes have been awarded for "sufficient progress" even without full resolution
- This would require at minimum: a new bound h(N) ≤ √N + o(N^{1/4}), proven rigorously
- The dispersion paper draft CLAIMS this (√(log N) improvement) but the proof is not rigorous

## Verdict

**Prize eligibility: NO — not at current state.**

The formalization (Track 1) is honest, high-quality infrastructure but contains no new mathematics. The spectral note (Track 2) is a finite computation. The dispersion paper (Track 3) contains the only potentially new claim (Theorem 3), but:

- The proof is not rigorous as written
- It has not been scored, peer-reviewed, or Lean-formalized
- Even if correct, it falls far short of resolving the conjecture

### If the Dispersion Paper Were Made Rigorous

If Theorem 3 could be rigorously proved (the √(log N) improvement), the scoring would change from A0 to approximately **A2** ("compute, bound, extend — restricted"). This would be a publishable result and potentially the first improvement to the Lindström-type bound since 1969, which would be notable. But it would NOT qualify for the $1,000 prize, which requires resolution.

### Recommended Path Forward

1. **Immediate value (no prize):** Post the formalization to erdosproblems.com (already ready). This establishes presence, invites collaboration, and creates visibility. The erdosproblems_post_FINAL.md is correctly scoped and ready.

2. **Medium-term (possible publication):** Get the dispersion paper draft peer-reviewed by a combinatorialist. If Theorem 3 can be made rigorous, submit to a journal. Score it via erdos_scorer.py. This would move to A2/B2/C1.

3. **Long-term (prize-eligible):** The PPV reformulation, even if currently incomplete, identifies a new proof skeleton (information-theoretic rather than Cauchy-Schwarz). If this approach can be pushed further — perhaps by combining with the spectral rigidity connection — it could potentially close the gap. But this is speculative.

## Score Update Recommendation

| Track | Current Score | Proposed Score | Rationale |
|---|---|---|---|
| Classical formalization | A0 / B5 / C1 | A0 / B5 / C1 | No change — this is correctly assessed |
| Dispersion paper | UNSCORED | A2 / B0 / C1 (IF proven) | "Compute, bound" — new bound if rigorous; B0 because no Lean formalization |
| Spectral choreography | UNSCORED | A0 / B1 / C1 | Finite verification only; proof architecture but no compiled theorems |

---

*Assessment generated 2026-04-15. Bound by Formalization Integrity Protocol and Score-to-Verb Safety Table.*
