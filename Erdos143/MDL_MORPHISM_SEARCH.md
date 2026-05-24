# Orthogonal MDL-Morphism Search — Erdős #143 part (ii)

**Date:** 2026-05-16 · **Verdict:** FURTHER_WORK_NEEDED (no Stage-4 closure; two viable research threads identified)

## Why this search exists

KLL 2025 is the state of the art on Erdős #143: their GCD-graphs + Selberg-sieve framework proves $S_1(N) = \sum_{x \in A, x \le N} 1/x = o(\log N / \sqrt{\log \log N})$ for every well-separated set $A$. By Abel summation, this is exactly $(\log\log N)^{1/2 + \epsilon}$ short of forcing parts.ii ($\sum 1/(x \log x) < \infty$).

Two pathways forward: (a) strengthen KLL inside their framework (very hard — that's exactly what the literature is stuck on), or (b) take an **orthogonal** approach via a different morphism that bypasses the sieve machinery. This search ran pathway (b).

## The morphisms evaluated

Seven candidates were defined and triaged. M1 is the channel-occupancy baseline (validated in the prior session pass). M2–M7 are orthogonal candidates:

| # | Morphism | Gate 1 | Stage 3 signal | Verdict |
|---|---|---|---|---|
| M1 | Channel-occupancy (baseline) | GREEN but KLL-equivalent | (= KLL) | baseline only |
| M2 | Multiplicative Sidon analogue | YELLOW | STRONG: $\|A\|/\sqrt{N}$ decays to 0.15 at $N=10^8$ | promising but technically blocked |
| M3 | Log-position aperiodic point process | GREEN | MIXED: small log-shift lattice modes elevated, larger shifts suppressed | genuinely novel, messy |
| M4 | Anti-Beatty | RED | not tested | rephrasing, dropped |
| M5 | Shannon entropy on $\mathbb{1}_A$ | GREEN | WEAK: ~1 bit per element, no sub-counting structure | orthogonal philosophy, no signal |
| M6 | Koopman / transfer-operator spectral | GREEN | not tested (build cost too high for session) | most orthogonal, most ambitious |
| **M7** | **VN entropy on bipartite dilation channel** | **RED (added follow-up)** | **Sub-maximal entanglement: $S_{VN}/\log\|A\| \to 0.05$, effective rank → 1.4** | **NOTATIONAL — WS arithmetic washed out by universal quantum inequalities** |

Perplexity ranking from the M2–M6 comparative pass (1 = most promising): **M6 #1, M3 #2, M2 #3, M5 #4, M4 #5**. M7 was evaluated in a separate follow-up pass and rated NOTATIONAL.

## The M7 follow-up: an important negative result

After the initial six-candidate pass, M7 was added: define a bipartite Hilbert space $\mathcal{H}_A \otimes \mathcal{H}_K$, encode well-separation as a sparsity constraint on the bipartite amplitude matrix, and apply von Neumann entropy / strong subadditivity / Araki–Lieb. This connects directly to the H² portfolio's `PinchingIdentity.lean` (the §101 patent shield establishing $S_{VN}(\rho) = H_{Shannon}(\text{eigenvalues})$).

Both signals converged on NOTATIONAL:

- **Perplexity Q2 verdict:** "M7 is a notational reformulation that repackages the same eigenvalue constraints, not an orthogonal route... SSA and pinching do not source any new inequality specific to the WS arithmetic structure; they give generic entropy inequalities valid for all bipartite states with a given Gram matrix. The arithmetic content is entirely in bounding the Gram-matrix entries, which is precisely where KLL already invoke Selberg's sieve."
- **Numerical experiment (greedy WS, indicator α):** the bipartite state is approximately a product state. $S_{VN}/\log|A| \to 0.05$ at $|A| = 1000$; effective rank $\to 1.4$. Mutual information $I(A:K) \approx 0.66$ bits — nearly zero entanglement. But this near-product structure is *not* WS-specific — it's a generic property of the dilation indicator (for $x$ near $N$, only $k=1$ produces a valid dilate). Any sparse set gives the same near-product structure. WS arithmetic is washed out.

This is an honest negative result for the H² portfolio: the KvN / VN-entropy infrastructure is genuinely valuable for §101 patent shielding and for spectral-physics problems, but it does NOT translate into a new attack on parts.ii via universal quantum-information machinery. The arithmetic structure of WS sets lives in congruences and GCDs, which universal entropy inequalities don't see.

Where M7 might still work: a **non-universal** quantum-information inequality specific to multiplicative number theory — currently undeveloped. Building one would be a research project comparable in scope to M6's Koopman framework. Both routes require establishing arithmetic-specific operator theory rather than borrowing generic quantum-info tools.

## Why M2 looked best numerically but isn't the right thread

The $\|A \cap [1,N]\| / \sqrt{N}$ ratio for greedy near-WS at $N=10^8$ is **0.15** — well below the Sidon threshold of $O(1)$. Across decades, the ratio decays monotonically from 2.5 (at $N=10^3$) to 0.15 (at $N=10^8$). This is **sub-Sidon behavior** experimentally: greedy WS sets are sparser than the M2 conjecture would require.

But as Perplexity flagged: the analytical path for proving a Sidon-rate bound for WS sets converges back to KLL's GCD-graph counting unless a genuine sum–product or new incidence bound materializes. Current technology gives no such bound. M2 is the right experimental signal but the wrong analytical hook.

## Why M6 is the most promising research direction

M6 is the only candidate that is **structurally** different from KLL: instead of counting edges in a GCD graph, it asks for a spectral gap of the dilation transfer operator $T_k$ acting on a suitable function space. If such a gap exists, exponential decay of correlations follows, which would imply much stronger bounds than KLL achieves — convergence of $\sum_{x \in A} x^{-s}$ for $s$ near 1, far beyond parts.ii.

The reason no one has done this is operational: there's no canonical Markov partition or symbolic coding for WS sets the way there is for continued fractions / Gauss map. Building one is itself a research project (~6–12 months for an expert in thermodynamic formalism).

**Portfolio bias caveat:** The H² portfolio favors Koopman / von Neumann framings (RPQEC, Mendoza Conjecture, KvN bridge). Perplexity rank-1 for M6 may be inflated by this bias. The Gate-5 verdict explicitly acknowledges this and discounts accordingly — but even discounted, M6 remains the cleanest orthogonal architecture in the menu.

## What about M3 (log-process)?

M3 is genuinely orthogonal but the experimental signature was mixed. Lattice-mode power at small integer log-shifts $t = 1, 2, 3$ came in at 8.78, 2.64, 5.98 — well above the Poisson baseline of 1.0 — while $t \ge 4$ mostly suppressed. This suggests greedy WS construction has natural log-spacing preferences at ratios $e, e^2, e^3$ that aren't aperiodic.

This is a real signal but may be a **greedy-artifact**: the greedy construction samples log-uniformly, so log-positions inherit some uniform-spacing structure. The right test is M3 on an explicit number-theoretic WS set (e.g., $\{\lfloor \alpha^n \rfloor : \alpha > 1\text{ with strong Diophantine properties}\}$) where the construction doesn't bias the log-positions. Deferred to follow-up.

## What about M5 (entropy)?

M5's prediction was entropy-rate zero for $\mathbb{1}_A$, which would imply counting function $n^{o(1)}$ and parts.ii would follow. The numerical entropy per element comes out to ~1 bit (consistent with finite positive entropy rate), so this prediction isn't supported by the greedy construction.

The honest reading is that WS is a **global multiplicative** constraint, not a finite-radius local pattern, so standard subshift-of-finite-type entropy machinery doesn't apply directly. Adapting Ornstein-Weiss / Gurevich entropy to nonlocal multiplicative interactions is its own research project.

## What about M4 (anti-Beatty)?

Dropped at Gate 1. Perplexity made the case crisply: "the missing ingredient is a transference principle connecting large harmonic mass to Beatty correlation in this multiplicative context; there is no off-the-shelf machinery." Without that transference, M4 reduces to known density constraints — no novelty.

## Stage 4 (Lean) decision

**SKIP.** No candidate produces a clean conditional reduction lemma in the sense of `lemmas.parts_ii_of_strengthened_S1` (line 130 of 143.lean). Adding speculative conditional lemmas for M3/M5/M6 hypotheses would inflate the file with unproven, unmotivated abstractions. The existing conditional reduction lemma already crystallizes the KLL-quantitative gap; the M-candidates point in different architectural directions that aren't reducible to a single one-line Lean statement.

If the M6 thread later produces a precise "spectral-gap hypothesis ⟹ parts.ii" statement, that would warrant a new conditional reduction lemma. Not in this pass.

## Honest scope

This is **research-direction work, not a proof**. The orthogonal search:

- ✅ Surveyed five alternate morphism architectures against KLL
- ✅ Ran prior-art checks on each (no rediscovery of an existing failed attempt)
- ✅ Identified the experimental signature each predicts and tested 3 of them
- ✅ Produced a ranked menu suitable for follow-up
- ❌ Did not close parts.ii (which remains open)
- ❌ Did not produce a new conditional reduction lemma worth formalizing
- ❌ Did not bypass KLL — KLL stays the state of the art

## Files

- `MDL_MORPHISM_INVENTORY.json` — structured per-candidate record with Gate 1/2/3 verdicts
- `MDL_RESEARCH_DIRECTIONS.md` — hand-off-ready single-page descriptions of M2 and M6 threads
- `MDL_perplexity_query.md` + `MDL_perplexity_comparative.md` — Stage 1 prior-art pass query and full Perplexity response
- `MDL_experiments.py` + `MDL_experiments_results.json` — Stage 3 numerical experiments (M2/M3/M5)
- `.claude/skills/erdos-problem-solver/state/143_morphism_candidate_M1.json` — baseline channel-occupancy (preserved)
- `.claude/skills/erdos-problem-solver/state/143_pipeline_state.json` — pipeline state extended with `mdl_morphism_search` key

## Pipeline state update

`143_pipeline_state.json` gets a new top-level `mdl_morphism_search` section with sub-keys: `started_at`, `completed_at`, `candidates_evaluated`, `gate_5_verdict`, `most_promising_thread`, `artifacts[]`. Base pipeline state (`current_stage=7 DONE`) is unchanged — this is an extension, not a reset.

## Next actions (your call)

1. **External hand-off** — reach out to Lamzouri or Lichtman with `MDL_RESEARCH_DIRECTIONS.md` to ask: (a) is the M6 Koopman framing genuinely orthogonal to KLL's sieve, or already implicit? (b) is the M2 sub-Sidon experimental signature compatible with their density bounds? Publisher gate + Crackpot-Scrub before sending.
2. **Follow-up session on M3** — re-run the log-process autocorrelation experiment on explicit number-theoretic WS constructions (Beatty-like with strong Diophantine $\alpha$) to verify greedy-artifact hypothesis.
3. **Follow-up session on M6** — design a finite-sample dilation transfer operator (e.g., on a finite WS subset of $[1, 1000]$) and compute its spectrum to give M6 a falsifiable numerical signature.
4. **Drop the thread** — accept FURTHER_WORK_NEEDED and revisit only if KLL gets strengthened or a new sieve technique appears.
