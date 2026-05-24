# Research Directions — Erdős #143 part (ii)

Two orthogonal threads identified by the MDL-morphism search (2026-05-16). Both target the open Erdős conjecture: $\sum_{x \in A} 1/(x \log x) < \infty$ for every well-separated set $A \subseteq (1, \infty)$ (countably infinite with $|kx - y| \ge 1$ for all $x \ne y \in A$ and all integers $k \ge 1$).

Reference baseline: Koukoulopoulos–Lamzouri–Lichtman 2025 (arXiv:2502.09539), which proves $\sum_{x \in A, x \le N} 1/x = o(\log N / \sqrt{\log \log N})$ via GCD-graphs + Selberg sieve. By Abel summation, KLL falls a $(\log\log N)^{1/2 + \epsilon}$ factor short of forcing parts.ii.

---

## Thread 1 — Koopman / transfer-operator spectral architecture (theoretical, highest ceiling)

**Idea.** Replace the GCD-graph counting with a spectral analysis. Define a dilation transfer operator $T_k$ on a suitable function space (a quasi-compact operator on a Banach space carrying a Markov partition of well-separated sets, in the style of Baladi–Mayer's treatment of the Gauss map for continued fractions). If $T_k$ admits a spectral gap, the resulting exponential decay of correlations bounds the count of near-collisions $|kx - y| < 1$ exponentially in the gap parameter — much stronger than KLL's $o(\log N / \sqrt{\log \log N})$.

**What this would yield.** If the spectral gap is genuine, $\sum_{x \in A} x^{-s}$ would converge for some $s < 1$ — *far* stronger than parts.ii. parts.ii follows as a side corollary; the stronger Erdős bound $S_1(N) \ll \log N / \sqrt{\log \log N}$ follows too.

**The hard part.** No canonical Markov partition for well-separated sets exists in the literature. Building one is the research project. By analogy: continued fractions have the Gauss map as a natural symbolic coding; well-separated sets under integer dilation have no obvious analogue. Candidate constructions:

1. Encode $A$ as a sequence of "log-scale shifts" under multiplication by integer $k$, and define $T_k$ on the orbit space.
2. Set $T_k f(x) = f(kx) / k$ (transfer-operator dual to multiplication) and ask for the function space on which this is quasi-compact.
3. Use thermodynamic formalism (Dolgopyat-type estimates) on the underlying multiplicative shift.

**Why it might fail.** The dilation action $k \mapsto kx$ has no natural invariant measure on $(1, \infty)$ comparable to Gauss measure on $[0, 1)$. The spectral gap may not exist, or may be an artifact of any chosen "canonical" measure.

**Why it might succeed.** Multiplicative number theory has precedent in Tao–Teräväinen on Sarnak's conjecture and Frantzikinakis–Host's ergodic approach to multiplicative functions. The spectral framework for *multiplicative* dynamical systems is genuinely developing. The orthogonality to KLL is real: their proof never invokes spectral theory of operators.

**Experimental probe (not yet done).** Construct a finite-sample transfer operator on a WS subset of $[1, 1000]$. Compute its spectrum. If the gap is large and stable across sample sizes, the framework is on track. If spectrum is dense or gap shrinks with sample size, the construction is the wrong one.

**Key references.**
- KLL 2025 (arXiv:2502.09539) — what to beat.
- Baladi, *Positive Transfer Operators and Decay of Correlations* — the toolkit.
- Tao, Sarnak's conjecture work — multiplicative ergodic methods.
- Frantzikinakis–Host, "The logarithmic Sarnak conjecture" — closest existing multiplicative-ergodic analogue.

**Estimated effort.** 6–18 months for a researcher fluent in thermodynamic formalism + analytic number theory. Not session-scale.

---

## Thread 2 — Multiplicative Sidon analogue (experimental signal, technically blocked)

**Idea.** Treat well-separated sets as a multiplicative analogue of Sidon (B₂) sets. The condition $|kx - y| \ge 1$ for all $k \ge 1$ resembles a "no near-collision under integer dilation" constraint, the multiplicative cousin of the additive Sidon constraint $a + b \ne c + d$. If the analogue holds quantitatively, $|A \cap [1, N]| \ll \sqrt{N}$ (Sidon rate), which makes parts.ii trivial.

**Experimental support (this session).** Greedy near-WS construction at $|A| = 1500$, scaled to $N = 10^8$:

| $N$ | $\|A \cap [1,N]\|$ | $\|A\|/\sqrt{N}$ |
|---|---|---|
| $10^3$ | 80 | 2.53 |
| $10^4$ | 246 | 2.46 |
| $10^5$ | 544 | 1.72 |
| $10^6$ | 864 | 0.86 |
| $10^7$ | 1201 | 0.38 |
| $10^8$ | 1500 | 0.15 |

The ratio decays monotonically and stays well below the Sidon threshold. Greedy WS is sub-Sidon. Not a proof — could be a greedy artifact — but the signal is consistent across six decades.

**The hard part.** The analytical handle is missing. Classical additive Sidon machinery (Cilleruelo, Ruzsa, Green–Ruzsa Fourier-analytic) targets additive collisions $a + b = c + d$, not the multiplicative-dilation condition $|kx - y| < 1$. Adapting requires either:

1. A new incidence/energy bound for the relation $\{(x, y, k) : |kx - y| < 1\}$.
2. A multiplicative analogue of Fourier-analytic Sidon bounds — currently no such analogue is established with Sidon-quality sharpness.
3. A sum–product approach on $\mathbb{R}$ or $\mathbb{Q}$ to show that too many near-collisions force unintended additive structure.

**Why it might fail.** Once you reduce a Sidon-style estimate to counting double-collisions $|kx_1 - y_1| < 1$ and $|kx_2 - y_2| < 1$, you tend to recover almost exactly the same inequalities KLL use for bounding GCD-graph edges. The "Sidon framing" then becomes a relabeling rather than a new bound.

**Why it might succeed.** The numerical signal *is* sub-Sidon, not merely KLL-rate. If the gap between the greedy experimental rate ($\|A\|/\sqrt{N} \to 0$) and the KLL rate ($\|A\|/\log N \to 0$) is genuine for all WS sets, the analytical tools for capturing the Sidon-quality bound must exist somewhere — even if not in the current incidence-geometry literature.

**Experimental follow-up (not yet done).** Test the $\|A\|/\sqrt{N}$ ratio on explicit non-greedy WS constructions to confirm the signal isn't a sampling artifact. Possible test sets:
- Number-theoretic WS like $\{p_n^2 : p_n \text{ prime}\}$ if WS under integer dilation
- Algebraic-irrational generated sets like $\{\lfloor \alpha^n \rfloor\}$ for $\alpha$ with strong Diophantine properties
- WS sets constructed by inverse-image under integer-multiplication

**Key references.**
- Cilleruelo, "Sidon sets in $\mathbb{N}^d$" — additive B₂ machinery.
- Erdős–Bateman 1962, J. Number Theory — multiplicative Sidon precursors.
- Ruzsa, B₂-set surveys — Fourier-analytic approach.
- KLL 2025 — direct comparison: their bound vs the Sidon-rate conjecture.

**Estimated effort.** 3–12 months IF a Sidon-style sieve technique materializes. Unbounded otherwise.

---

## Two-thread strategy

These threads are complementary:

- **Thread 1 (M6 Koopman)** is the high-risk, high-reward research direction. Builds new infrastructure (transfer operator + function space) that could generalize beyond #143.
- **Thread 2 (M2 Sidon)** is the high-signal, low-machinery direction. Falsifiable in months via experimentation; either yields the bound or shows the limits of additive-style methods.

Recommended pursuit order: **Thread 2 experimentally** (cheap follow-up sessions to verify the sub-Sidon signal isn't artifactual), then **Thread 1 theoretically** (requires substantial commitment).

## Candidate external collaborators

- **Dimitris Koukoulopoulos** (KLL lead) — best-positioned to evaluate whether either thread is genuinely orthogonal to their sieve framework or already implicit in their proof structure. Email: dimitris.koukoulopoulos@umontreal.ca (per Univ. Montreal directory; verify before sending).
- **Youness Lamzouri** (KLL co-author) — multiplicative function specialist; specifically qualified to assess M6 Koopman framing.
- **Jared Lichtman** (KLL co-author) — primitive-sets specialist; specifically qualified to assess M2 Sidon analogue.
- **Terence Tao** (long shot, but Tao's blog has covered #143 and adjacent multiplicative dilation problems) — would have an opinion on whether the Sidon-rate conjecture is plausible.
- **John Baez** (#baez persona) — for category-theoretic adversarial read on whether the M6 Koopman framing is a genuine new attack or a categorial relabeling.

**Outreach prerequisites.** Each email through the H² Publisher pipeline + Crackpot-Scrub. The technical content is general (no patent surface, no portfolio-specific machinery), so the gate should pass cleanly. Disclose use of Perplexity for the prior-art pass per AI-acknowledgment standards.

## What stays open

Erdős #143 part (ii) remains an open conjecture. KLL is state of the art. This document is a research-direction hand-off, not a proof.
