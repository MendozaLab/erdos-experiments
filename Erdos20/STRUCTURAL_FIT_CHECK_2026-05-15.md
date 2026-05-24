# Structural Fit Check: WS-01-CENTER-STRIP-CANCELLATION (EHP114) applied to Erdős #20

Date: 2026-05-15
Plan track: Probe 3 (WS-01 cross-problem cold-test extension) of `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)"
Role: **Explicit negative control.** Confirms that the structural-fit checker is not permissive — i.e. it can correctly identify a problem where WS-01 does not apply.
Claim ceiling: **Internal structural-fit diagnostic only.** Not a proof, not a Lean statement, not a solution claim, not a public-facing artifact.

## Problem statement (canonical, verbatim)

From `Math/formal-conjectures/FormalConjectures/ErdosProblems/20.lean` lines 62–83:

> Let $f(n,k)$ be minimal such that every $F$ family of $n$-uniform sets with $|F| \geq f(n,k)$ contains a $k$-sunflower. Is it true that $f(n,k) < c_k^n$ for some constant $c_k > 0$ and for all $n > 0$?

A sunflower with kernel $S$ is a collection of sets in which all possible distinct pairs of sets share the same intersection $S$ (the "petals" outside the kernel are pairwise disjoint). The conjecture is marked `research open`, AMS category 5 (combinatorics, not analysis). References: erdosproblems.com/20, Wikipedia Sunflower (mathematics).

## Parameter space

- $n, k \in \mathbb{N}$ — discrete integer pair indexing the function $f$
- $F$ ranges over set families of $n$-uniform sets (set of sets, each set has cardinality $n$)
- $S$ is a fixed kernel set

The parameter space is **discrete and combinatorial**. There are no real or complex coordinates, no polynomial coefficients, no continuous deformation parameter. The objects of interest are integer-valued functions ($f(n,k)$, set cardinalities) and set-theoretic predicates (pairwise intersection equals $S$).

For comparison, the EHP114/EHP1041/EHP1043 parameter space is monic complex polynomials of degree $n$ — real-analytic in $2n$ real coordinates, with a well-defined zero curve of the polynomial's level set $\{|f(z)| = 1\}$ in $\mathbb{C}$.

## Feature 1 — Zero-curve / level-set boundary

**Verdict: NO_STRUCTURAL_FIT**

There is no $F(z) = 0$ zero-curve. The extremal object is a *combinatorial set family* witnessing the lower bound $f(n,k)$, not a level set of a real-analytic function. The conjecture asks about the asymptotic growth rate of an integer function indexed by two integers. No real-analytic or complex-analytic function is involved at any step of the problem definition. There is no manifold, no codimension, no normal direction, no tangent space.

The closest analog one could attempt is: define an indicator function on the discrete space of $n$-uniform set families and ask about its level set. But (a) such an indicator is not continuous (let alone analytic) in any natural topology on the discrete space, (b) the conjecture's quantity $f(n,k)$ is the *minimum* set-family cardinality witnessing a property, which is computed by an integer-valued $\sInf$ operation, not by tracking a continuous boundary curve, and (c) no interval-arithmetic-style enclosure machinery would apply to an integer-valued function on a finite/discrete space anyway.

This is the *intended* mismatch — the brief specifies #20 as the explicit negative control. The structural-fit checker correctly returns NO at the first feature.

## Feature 2 — Wall-separation gate / coercivity condition

**Verdict: NO_STRUCTURAL_FIT**

A wall-separation gate requires a continuous parameter space with a smooth (or at least real-analytic) boundary curve separating two regions, and a gradient bounded away from zero perpendicular to the boundary so the curve does not develop critical points in the gate's collar. The sunflower problem has none of these:

- No continuous parameter space (set families are discrete, $f(n,k)$ is integer-valued)
- No boundary curve to apply the wall to
- No gradient of any continuous quantity
- No collar geometry, no normal direction

The Razborov / Naslund-Sawin / Alweiss-Lovett-Wu-Zhang advances on this problem use combinatorial / probabilistic / pseudorandom techniques (e.g. spread approximations) that have no analog in WS-01's branch-point + interval-IVT toolkit.

## Feature 3 — Branch-point validation via interval IVT + bounded-away gradient

**Verdict: NO_STRUCTURAL_FIT**

Interval IVT requires a continuous real-valued function whose interval-enclosure can be checked for sign changes inside a box. The sunflower problem has no such function:

- $f(n,k)$ is integer-valued
- Set-family cardinalities are integer-valued
- The Sunflower predicate is a discrete logical condition, not a real-valued quantity

Bounded-away gradient requires a smooth function with a non-vanishing derivative. There is no derivative to bound.

## Aggregate verdict

**`NO_STRUCTURAL_FIT — COLD TEST SKIPPED`**

All three features are NO. The sunflower problem is a discrete combinatorial extremal-set problem; WS-01 is an interval-arithmetic wall-separation rewrite for real-analytic level-set boundaries. The two architectures address fundamentally different mathematical objects.

This is the **intended outcome for the negative control**. The structural-fit checker correctly identifies that WS-01 does not apply to #20 — confirming the checker is not permissive. Without this negative control we could not distinguish "WS-01 generalizes to the entire Erdős corpus" from "the structural-fit checker says YES to everything regardless of whether the architecture applies."

## Recommended next step

**Cold test skipped.** No box-data generation or probe invocation is performed. This document, the accompanying `_RESULTS.json`, and the `.sha256` are the entire artifact triple for #20 in this cross-problem extension.

The existing #20 sunflower work (`Math/erdos-experiments/Erdos20/EXP-MATH-ERDOS20-PER-CORE-CLOSURE-*`, `EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-*`, `SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md`) uses transfer-matrix and combinatorial-core-closure techniques that are appropriate for the discrete structure and orthogonal to WS-01. Those artifacts are not touched by this fit check.

## Sources cited

Local files read for this fit check:

- `Math/formal-conjectures/FormalConjectures/ErdosProblems/20.lean` (lines 1-86) — canonical Lean problem statement
- `Math/erdos-experiments/cross_problem_ws01_applicability_probe.py` (lines 1-300) — probe schema, confirming what inputs the architecture requires
- `Math/erdos-experiments/Erdos1041/STRUCTURAL_FIT_CHECK_2026-05-15.md` (sibling reference, used for YES-verdict comparison)
- `Math/erdos-experiments/Erdos20/SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md` (referenced; existing #20 work uses orthogonal techniques)

External canonical reference (not refetched this session, locked via Lean file metadata):

- erdosproblems.com/20 (Erdős's problem entry)
- Wikipedia: Sunflower (mathematics)

Plan source:

- `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)" Probe 3 (lines 318-326)
- `~/.claude/plans/yes-wondrous-blum-agent-a4ae5ad039883b7ba.md` (per-candidate execution method; #20 designated negative control)
