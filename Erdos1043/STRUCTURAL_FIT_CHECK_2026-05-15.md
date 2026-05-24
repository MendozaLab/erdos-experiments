# Structural Fit Check: WS-01-CENTER-STRIP-CANCELLATION (EHP114) applied to Erdős #1043

Date: 2026-05-15
Plan track: Probe 3 (WS-01 cross-problem cold-test extension) of `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)"
Claim ceiling: **Internal structural-fit diagnostic only.** Not a proof, not a Lean statement, not a solution claim, not a public-facing artifact.

## Problem statement (canonical, verbatim)

From `Math/formal-conjectures/FormalConjectures/ErdosProblems/1043.lean` lines 30–67:

> Let $f \in \mathbb{C}[x]$ be a monic polynomial. Must there exist a straight line $\ell$ such that the projection of $\{ z : |f(z)| \leq 1 \}$ onto $\ell$ has measure at most $2$? Pommerenke [Po61] proved that the answer is no.

The conjecture is marked `research solved`, negatively, by Pommerenke 1961. The weak variant (measure ≤ 3.3 on some line) does hold. References: EHP58 (Erdős–Herzog–Piranian 1958), Po59/Po61 (Pommerenke 1959/1961). Already formalized in Lean by Alexeev using Aristotle.

## Parameter space

- Monic complex polynomials $f \in \mathbb{C}[X]$ of arbitrary degree ≥ 1. The constraint "roots in the unit disk" is *not* part of the canonical #1043 statement (in contrast to #114 and #1041, which both require roots in the closed unit ball). The level set $\{|f(z)| \leq 1\}$ depends on the polynomial's coefficients without any disk constraint on the roots.
- Zero curve: $\{|f(z)| = 1\}$, the lemniscate. Same curve family as #114/#1041, but the parameter space is wider (no root-disk constraint).
- Objective functional: $\inf_{\|u\|=1} \mu_1(P_u(\{|f| \leq 1\}))$, the minimum 1D Lebesgue measure of the *orthogonal projection* of the closed sublevel set onto some line $\ell = \mathbb{R} \cdot u$. The projection direction $u \in S^1$ adds an extra unit-vector parameter.

The structural observation: #1043 lives in the same EHP-paper lineage as #114 and #1041, with the same lemniscate boundary, but the objective is a *projected* measure of the *closed* sublevel set rather than a sublevel-set-internal quantity (component count, path length). The projection direction $u$ is an extra parameter that does not affect the zero curve itself but adds a 1D family of derived measurements.

## Feature 1 — Zero-curve / level-set boundary

**Verdict: YES_STRUCTURAL_FIT**

The boundary $\{F = 0\}$ with $F(z) := |f(z)|^2 - 1$ is again the lemniscate, the same curve as #114/#1041. The closed sublevel set $\{|f(z)| \leq 1\}$ is bounded by the same curve as the open sublevel set $\{|f(z)| < 1\}$ (the closure vs. open distinction matters for the projection's measure but not for the curve's identity). The zero curve is real-analytic away from critical points $\{p'(z) = 0\}$.

One caveat versus #114/#1041: the absence of a root-disk constraint means the lemniscate can be unbounded in $\mathbb{C}$ (for polynomials with large leading coefficient's behavior, $\{|f| \leq 1\}$ is bounded but for low-degree or specific coefficient choices the sublevel-set component count and projected measure can vary widely). The zero curve is well-defined locally for any monic $f$; WS-01's per-box wall test is local, so the global unboundedness does not block per-box gate machinery. The fit at the *local* level is identical to #114.

## Feature 2 — Wall-separation gate / coercivity condition

**Verdict: PARTIAL_STRUCTURAL_FIT**

The wall-separation gate operates locally on a box on $\{F = 0\}$ regardless of objective. So WS-01's wall LHS/RHS comparison ports unchanged. *However*, the connection between certifying the lemniscate's local geometry and bounding $\mu_1(P_u(\{|f| \leq 1\}))$ is more indirect than in #114:

- For #114 the lemniscate's connected components map directly to the sublevel-set components whose lengths are summed.
- For #1041 the lemniscate is the boundary along which a path is constructed; the path-length conjecture lives inside the sublevel set.
- For #1043 the projection $P_u$ collapses the 2D sublevel set onto a 1D line; the projected measure depends on the *width* of the sublevel set across the projection direction at each line-parameter $s \in \ell$. Two points on the lemniscate that share a common projection contribute to the projection's measure as a single interval.

What this means structurally: a wall-separation gate that certifies the local geometry of $\{F = 0\}$ does not, by itself, bound the projected measure — Pommerenke's 1961 negative result is *precisely* a construction of polynomials whose lemniscate-bounded sublevel set has projection > 2 in every direction. WS-01 can still certify the *wall* (the curve's local non-degeneracy), but the conjecture's failure mechanism (Pommerenke-style adversarial constructions) is upstream of where WS-01 helps. The gate machinery applies; its diagnostic value for the conjecture is downstream-limited.

Hence `PARTIAL`: the architecture ports, but the cold test on #1043 measures wall-separation portability across an objective change, not WS-01's ability to close the conjecture. This is the *correct* test for the cross-problem portability question.

## Feature 3 — Branch-point validation via interval IVT + bounded-away gradient

**Verdict: YES_STRUCTURAL_FIT**

The branch-point gate per failure box uses the same two interval-arithmetic facts as #114/#1041:

1. $F_{iv}$ straddles zero (interval IVT)
2. $|\nabla F|$ bounded away from zero at box center

Both are local properties of the polynomial-derived $F(z) = |f(z)|^2 - 1$, evaluable with the same `mpmath.iv` or `inari` machinery regardless of which sublevel-set-derived quantity is the objective. The projection direction $u$ is parameterized over $S^1$ but does not enter the gate's box-level computation; it only enters the downstream measure-computation step.

The gate certifies branch points of the lemniscate; the projection direction $u$ produces a 1-parameter family of *secondary* questions about how the certified lemniscate-local geometry projects onto $\ell$. The gate machinery is unchanged.

## Aggregate verdict

**`PARTIAL_FIT_PROCEED — WS01_APPLICABLE_POSSIBLE`**

Features 1 and 3 are full YES. Feature 2 is PARTIAL because the wall test's diagnostic value for the conjecture is downstream-limited by the projection-measure objective (Pommerenke's 1961 negative result lives outside WS-01's reach). The gate machinery itself ports cleanly; the cold-test exercise on #1043 tests *whether the architecture survives an objective change on the same boundary geometry*.

The expected cold-test outcome: WS-01 likely closes a similar fraction of boxes as on #1041, because the wall computation depends only on the polynomial and the local box geometry, not on the projection direction. The diagnostic value of that closure for the *conjecture* is lower than on #1041, but the diagnostic value for the *portability question* is exactly the same.

## Recommended next step

**Proceed to cold-test step B**, using the same kernel as #1041 plus a projection-direction parameter sweep over $u \in S^1$ at 6 discrete angles ($\theta \in \{0, \pi/6, \pi/3, \pi/2, 2\pi/3, 5\pi/6\}$). The projection direction is recorded for each box but does not enter the per-box gate.

## Sources cited

Local files read for this fit check:

- `Math/formal-conjectures/FormalConjectures/ErdosProblems/1043.lean` (lines 1-70) — canonical Lean problem statement
- `Math/erdos-experiments/cross_problem_ws01_applicability_probe.py` (lines 1-300) — probe schema and classifier
- `Math/erdos-experiments/Erdos1038/STRUCTURAL_FIT_CHECK_2026-05-14.md` (lines 1-139, used as comparative template)
- `Math/erdos-experiments/Erdos1041/STRUCTURAL_FIT_CHECK_2026-05-15.md` (lines 1-92, sibling reference)

External canonical reference (not refetched this session, locked via Lean file metadata):

- erdosproblems.com/1043 (Erdős's problem entry)
- Pommerenke, Ch., *On metric properties of complex polynomials*, Michigan Math. J. (1961), 97-115 — negative resolution
- Pommerenke, Ch., *On some problems by Erdős, Herzog and Piranian*, Michigan Math. J. (1959), 221-225
- Alexeev's Lean formalization at `https://github.com/plby/lean-proofs/blob/main/src/v4.24.0/ErdosProblems/Erdos1043.lean`

Plan source:

- `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)" Probe 3 (lines 318-326)
- `~/.claude/plans/yes-wondrous-blum-agent-a4ae5ad039883b7ba.md` (per-candidate execution method)
