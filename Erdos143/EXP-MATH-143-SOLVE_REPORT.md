# EXP-MATH-143-SOLVE — Erdős Problem 143 / DeepMind formal-conjectures PR

**Date:** 2026-05-16 · **Verdict:** PARTIAL PROGRESS — parts.i status upgraded; parts.ii **remains open** but the precise quantitative gap is now pinned down.

## Bottom line up front

This session did **not** solve parts.ii of Erdős Problem 143 — it remains an open problem. What we did do:

1. **parts.i** ($\liminf$ density $= 0$ for well-separated sets): upgraded from `research open` to `research solved`, citing Koukoulopoulos–Lamzouri–Lichtman (2025, arXiv:2502.09539). They prove the density-zero result alongside the $o(\log n)$ harmonic-sum bound, so the literature status is "solved on paper, Lean body still sorry."

2. **parts.ii** ($\sum 1/(x \log x) < \infty$): the precise quantitative gap blocking a proof is now identified. By Abel summation, KLL's bound $S_1(N) \ll \log N / \sqrt{\log\log N}$ yields only $S_2(N) \ll \sqrt{\log\log N}$ — unbounded in principle. Closure requires strengthening the loglog exponent in KLL from $1/2$ to $1+\epsilon$, i.e. $S_1(N) \ll \log N / (\log\log N)^{1+\epsilon}$. No such strengthening is known.

3. **New formal artifact:** added a conditional reduction lemma `erdos_143.lemmas.parts_ii_of_strengthened_S1` to the Lean file. It states the implication "strengthened KLL bound ⟹ parts.ii" formally. Body is `sorry`, but the statement crystallizes the gap and the docstring carries the Abel-summation proof sketch.

## The problem (recap)

A set $A \subseteq (1, \infty)$ is **well-separated** if it is countably infinite and satisfies $|kx - y| \ge 1$ for all $x \neq y \in A$ and all integers $k \ge 1$. Erdős asked whether well-separation forces:

| Part | Statement | Status (post-session) |
|------|-----------|-----------------------|
| i | $\liminf \lvert A \cap [1,x] \rvert / x = 0$ | **research solved** (KLL 2025) |
| ii | $\sum_{x \in A} 1/(x \log x) < \infty$ | research open (gap pinned) |
| KLL variant | $\sum_{x \in A,\, x < n} 1/x = o(\log n)$ | research solved (Lean body `sorry`) |

## The Abel-summation calculation (heart of the session)

Define $S_1(N) := \sum_{x \in A,\, x \le N} 1/x$ and $S_2(N) := \sum_{x \in A,\, x \le N} 1/(x \log x)$. Abel summation against the monotone weight $1/\log t$:

$$
S_2(N) \;=\; \frac{S_1(N)}{\log N} \;+\; \int_2^N \frac{S_1(t)}{t (\log t)^2} \, dt.
$$

**Plug in KLL's actual bound** $S_1(t) \ll \log t / \sqrt{\log\log t}$:

- First term: $S_1(N)/\log N \ll 1/\sqrt{\log\log N} \to 0$. ✓
- Integral: $\ll \int_2^N \frac{dt}{t \log t \sqrt{\log\log t}}$. Substitute $u = \log t$, then $v = \log u$: integral becomes $\int dv / \sqrt{v}$, which **diverges** like $2\sqrt{\log\log N}$.

So KLL's bound gives only $S_2(N) \ll \sqrt{\log\log N}$ — unbounded.

**Plug in hypothetical strengthening** $S_1(t) \ll \log t / (\log\log t)^{1+\epsilon}$:

- First term: $S_1(N)/\log N \ll 1/(\log\log N)^{1+\epsilon} \to 0$. ✓
- Integral: $\ll \int \frac{dt}{t \log t \, (\log\log t)^{1+\epsilon}}$. Same substitutions yield $\int dv / v^{1+\epsilon}$, which **converges** for $\epsilon > 0$. ✓

So **the precise gap** between KLL and parts.ii is a strengthening of the loglog exponent from $1/2$ to $1+\epsilon$ — i.e. by a factor of $(\log\log N)^{1/2+\epsilon}$.

## Numerical evidence

Greedy-random near-well-separated set construction (log-uniform sampling, online WS-check with $k_{\max} = 40$). Scaled up to $|A| = 2000$ elements in $(1.5, 10^{10})$:

| $N$ | count | density | $S_1$ | $S_1/\log N$ | $S_2$ | $S_2 / \sqrt{\log\log N}$ |
|-----|-------|---------|-------|---------------|-------|----------------------------|
| $10^3$ | 63 | $6.3 \times 10^{-2}$ | 0.669 | 0.097 | 0.208 | 0.150 |
| $10^6$ | 783 | $7.8 \times 10^{-4}$ | 0.742 | 0.054 | 0.218 | 0.134 |
| $10^{10}$ | 2000 | $2.0 \times 10^{-7}$ | 0.742 | 0.032 | 0.218 | 0.123 |

$S_2$ stabilizes around 0.22 across **six decades** of $N$. The ratio $S_2 / \sqrt{\log\log N}$ even decreases mildly. This is strong experimental evidence for parts.ii — but it is not a proof, and it cannot be one: a sufficiently exotic well-separated set could in principle violate the conjecture without being captured by greedy random sampling.

(Caveat: the final rigorous all-pairs WS-check up to $k = 50$ revealed 3 mild violations at $k \in \{44, 45, 47\}$ with sub-unit slack $0.16$–$0.35$. The construction is **near-well-separated** rather than exactly well-separated, since enforcing $k \to \infty$ in greedy construction is not feasible. The experiment is qualitative evidence, not certification.)

The cluster construction I tried first (consecutive integers in clusters at $N_k = N_0 \cdot \text{growth}^k$) is **invalid** for integer growth ratios: cross-cluster pairs have $y/x$ exactly integer at $k = \text{growth}^{j-i}$, violating WS. The spot-check passed only because it only sampled within-cluster pairs.

## New Lean artifact: conditional reduction

```lean
@[category API]
theorem erdos_143.lemmas.parts_ii_of_strengthened_S1
    (h : ∃ ε > (0 : ℝ), ∀ A : Set ℝ, WellSeparatedSet A →
      (fun n : ℕ => ∑' x : A, if (x : ℝ) < n then 1 / (x : ℝ) else 0)
        =O[atTop] (fun n : ℕ =>
          Real.log n / (Real.log (Real.log n)) ^ ((1 : ℝ) + ε))) :
    ∀ A : Set ℝ, WellSeparatedSet A →
      Summable fun (x : A) ↦ 1 / (x * Real.log x) := by
  sorry
```

The body is `sorry`. The hypothesis is currently unknown — strengthening KLL by $(\log\log N)^{1/2+\epsilon}$ is open. The docstring contains the Abel-summation proof sketch and clearly labels the hypothesis as "currently open."

Why include a sorry'd lemma? It crystallizes the gap as a formal statement. A future Lean prover working on parts.ii now knows: prove the strengthened $S_1$ bound (or any sufficient analogue) and the analytic Abel-summation argument formally, and parts.ii closes.

## Compile

`lake env lean FormalConjectures/ErdosProblems/143.lean` → exit 0. Exactly 4 `sorry` warnings (parts.i, parts.ii, KLL variant, new conditional lemma). No errors.

## What this is NOT

- **Not a proof of parts.ii.** parts.ii is an open Erdős problem. KLL is the state of the art and falls short by a $(\log\log N)^{1/2+\epsilon}$ factor.
- **Not a complete formalization of parts.i.** The Lean body is `sorry`; KLL's GCD-graph machinery is a 60-page paper, not session-scoped.
- **Not a guarantee of upstream acceptance.** The DeepMind formal-conjectures maintainers may or may not want the conditional reduction lemma in their canonical file — it goes beyond the original Erdős statement.

## What this IS

- **Honest status update for parts.i.** Now matches the literature.
- **Precise quantitative localization of parts.ii.** The gap is named: one loglog factor between KLL exponent $1/2$ and sufficient exponent $1 + \epsilon$.
- **Formal scaffolding for future work.** The conditional reduction lemma is a stable target.
- **Six decades of numerical evidence** consistent with parts.ii holding.

## Pipeline state

| Stage | Status | Notes |
|-------|--------|-------|
| 0 — Selection | done | user-specified |
| 1 — Morphism + numerical | done | WS3 greedy 2000 pts to $10^{10}$ |
| 2 — Strategy | done | Abel-summation analytic threshold |
| 3 — Compute | done | analytical |
| 4 — Lean edit + build | done | 4 sorry warnings, exit 0 |
| 5 — Gate 5 quorum | done | two rounds; both caught corrections / clarified gap |
| 6 — Artifacts | done | this report |

## Files

- Edited: `Math/formal-conjectures-143-pr-wt/FormalConjectures/ErdosProblems/143.lean` (branch `codex/erdos143-formal-conjectures-kll`)
- Backup: `…/143.lean.bak`
- Build log: `/tmp/erdos143_build3.log`
- Numerical scripts: `/tmp/erdos143_experiment.py`, `/tmp/erdos143_ws_experiment.py`, `/tmp/erdos143_ws3_big.py`
- Results JSON: `Math/erdos-experiments/Erdos143/EXP-MATH-143-SOLVE_RESULTS.json`
- This report: `Math/erdos-experiments/Erdos143/EXP-MATH-143-SOLVE_REPORT.md`
- SHA-256: `Math/erdos-experiments/Erdos143/EXP-MATH-143-SOLVE_RESULTS.sha256`

## Next actions (your call)

1. **Commit + upstream PR?** Open PR against DeepMind formal-conjectures with the parts.i status upgrade, parts.ii docstring, and conditional reduction lemma. Publisher gate + Crackpot-Scrub on docstrings first.
2. **Research thread for parts.ii (not session-scale):** can KLL's GCD-graph framework be sharpened by a $(\log\log N)^{1/2+\epsilon}$ factor? Or is there a structural argument (e.g., multiplicative-energy decomposition) that bypasses the Abel-summation route?
3. **Drop the conditional lemma?** If you'd rather keep the Lean file aligned strictly with Erdős's three original questions, I can remove the new lemma.
