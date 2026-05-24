# Structural Fit Check: WS-01-CENTER-STRIP-CANCELLATION (EHP114) applied to Erdős #1041

Date: 2026-05-15
Plan track: Probe 3 (WS-01 cross-problem cold-test extension) of `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)"
Claim ceiling: **Internal structural-fit diagnostic only.** Not a proof, not a Lean statement, not a solution claim, not a public-facing artifact. This is an in-portfolio decision document about whether to invest cold-test effort on #1041 as part of the cross-problem WS-01 portability boundary test.

## Problem statement (canonical, verbatim)

From `Math/formal-conjectures/FormalConjectures/ErdosProblems/1041.lean` lines 19–69:

> Let $f(z) = \prod_{i=1}^{n} (z - z_i) \in \mathbb{C}[x]$ with $|z_i| < 1$ for all $i$. Conjecture: Must there always exist a path of length less than 2 in $\{ z \in \mathbb{C} \mid |f(z)| < 1 \}$ which connects two of the roots of $f$?

Companion statement (the Erdős–Herzog–Piranian Component Lemma) at lines 41–54:

> If $f$ is a monic degree $n$ polynomial with all roots in the unit disk, then some connected component of $\{z \mid |f(z)| < 1\}$ contains at least two roots with multiplicity.

The component lemma is marked `research solved`, the path-length conjecture `research open`. Reference: Erdős, P. and Herzog, F. and Piranian, G., *Metric properties of polynomials*. J. Analyse Math. (1958), 125-148.

## Parameter space

- Monic complex polynomials $f \in \mathbb{C}[X]$ of degree $n$ with `f.rootSet ℂ ⊆ Metric.ball 0 1`.
- Real dimension: $2n$ root-location coordinates, modulo overall permutation symmetry — **identical** to the parameter space of EHP58 #114 (where the conjecture asks about the length of the *longest* connected component of the same sublevel set).
- Zero curve: $F(z) := |f(z)|^2 - 1 = 0$, i.e. the lemniscate $\{|f(z)| = 1\}$. **Identical** to #114.
- Objective functional: $\inf_f \big[ \mathcal{H}^1(\gamma) \big]$ where $\gamma$ is a path in the sublevel set joining two roots. The component lemma guarantees such a $\gamma$ exists; the path-length conjecture asks whether $\inf_\gamma \mathcal{H}^1(\gamma) < 2$ is universal. **Differs** from #114's component-length-sum objective.

The key structural observation: #114 and #1041 are *the same paper*, *the same parameter space*, *the same zero curve*, and differ only in which functional of the sublevel-set geometry is bounded (sum of component lengths vs. minimum path length between roots). WS-01 operates on the wall test for the F=0 curve, which is shared.

## Feature 1 — Zero-curve / level-set boundary

**Verdict: YES_STRUCTURAL_FIT**

The lemniscate $\{|f(z)| = 1\} = \{F(z) = 0\}$ with $F(z) := |f(z)|^2 - 1$ is the boundary of the sublevel set $\{|f(z)| < 1\}$. This is *literally the same* zero curve as #114. The objective functional (path length vs. component length) does not change which curve forms the boundary; it only changes which derived measurement of the sublevel set is being bounded. Every WS-01-style interval-arithmetic argument about the lemniscate that holds for #114 holds identically here. The codimension is 1 (a real algebraic curve in $\mathbb{C} \cong \mathbb{R}^2$), the curve is real-analytic away from critical points $\{\nabla F = 0\}$, and the parameter-space variation under perturbation of root locations is identical to #114's. There is no structural daylight between #114 and #1041 at this feature.

## Feature 2 — Wall-separation gate / coercivity condition

**Verdict: YES_STRUCTURAL_FIT**

The EHP114 wall test requires the zero curve to split cleanly across the normal direction in a thin collar at the box center, with no critical points $\nabla F = 0$ inside the collar. The bounded-away $|\nabla F|$ condition guarantees that $F = 0$ defines a smooth 1-manifold with a well-defined tangent direction $t(u)$.

For #1041 the same gate applies *unchanged*. The lemniscate $\{F = 0\}$ is the same curve; the critical points $\{\nabla F = 0\} = \{p'(z) = 0\}$ are the critical points of the polynomial, which sit inside the unit disk by Gauss-Lucas. The wall-separation requirement — that a box on the lemniscate does not contain a critical point and has $|\nabla F|$ bounded away from zero — is a property of the curve and the polynomial, not of the objective functional. WS-01's center-strip-cancellation rewrite (replacing $\sup_r |F(0,r)|$ with the quadratic-plus-cubic Taylor remainder at a validated branch point) operates on the same wall LHS/RHS comparison.

Note: the *path-length objective* introduces an additional question of how to certify a path through the sublevel set of length $< 2$. That question is downstream of the wall test; WS-01 only certifies the wall, not the path. So WS-01 ports to the wall layer of #1041 cleanly; the conjecture-closing step (constructing the short path) is a separate layer that WS-01 does not address. This is consistent with #114, where WS-01 certifies wall separation per cell and leaves the global component-length aggregation to the atlas-integration layer.

## Feature 3 — Branch-point validation via interval IVT + bounded-away gradient

**Verdict: YES_STRUCTURAL_FIT**

The EHP114 branch-point gate per failure box requires:

1. $F_{iv}$ straddles zero (interval IVT delivers a zero of $F$ in the box)
2. $|\nabla F|$ bounded away from zero at the box center (zero curve is a smooth 1-manifold; tangent $t$ well-defined)

For #1041 the inputs to both checks are identical to #114: $F(z) = |f(z)|^2 - 1$, $\nabla F = 2 \operatorname{Re}(\bar f \cdot f')$ horizontal + $2 \operatorname{Im}(\bar f \cdot f')$ vertical, both evaluable interval-arithmetically with `mpmath.iv` (or Rust `inari` at the certified-Rust scale). The center-strip-cancellation Taylor remainder $\tfrac{1}{2} |F_{tt}| R^2 + \tfrac{T_3}{6} R^3 + R_n$ uses second and third directional derivatives of $F$, which are again polynomial-derived and interval-evaluable. There is no part of the gate machinery that depends on the objective functional. The same `cross_problem_ws01_applicability_probe.py` consumer schema applies without modification.

## Aggregate verdict

**`YES_STRUCTURAL_FIT — PROCEED`**

All three features map identically onto #1041 because #1041 shares the parameter space, zero curve, and gradient structure of #114. The only difference is which sublevel-set-derived quantity is the conjecture's objective; that difference lives downstream of the WS-01 wall test and does not affect the gate machinery. The cold test should reproduce #114-like closure on a small Chebyshev-like base polynomial with a 2D box parameter sweep, modulo finite-N quantitative variation in the constants $|F_{tt}|$, $T_3$, $R_n$.

The strongest possible structural-fit verdict: **WS-01 should port to #1041 by direct copy with only the box-data input changing.** If the cold test fails despite this, WS-01 is more brittle than even paper-sibling portability would predict, which would itself be informative.

## Recommended next step

**Proceed to cold-test step B** (build candidate failure-box generation, then run the probe). No pre-condition gates needed beyond the standard interval kernel that already works for #1038 — the polynomial-level-set kernel from EHP114's reference implementation reuses cleanly.

## Sources cited

Local files read for this fit check:

- `Math/formal-conjectures/FormalConjectures/ErdosProblems/1041.lean` (lines 1-71) — canonical Lean problem statement
- `Math/erdos-experiments/Erdos114/EHP114_REGULAR_SLICE_WALL_TAYLOR_TARGET_2026-05-06.md` (referenced via #1038 structural-fit check) — EHP114 wall-separation gate failure pattern, used for direct comparison
- `Math/erdos-experiments/cross_problem_ws01_applicability_probe.py` (lines 1-300) — probe schema and classifier algorithm
- `Math/erdos-experiments/Erdos1038/STRUCTURAL_FIT_CHECK_2026-05-14.md` (lines 1-139) — #1038 fit-check template used as comparative reference

External canonical reference (not refetched this session, locked via Lean file metadata):

- erdosproblems.com/1041 (Erdős's problem entry)
- Erdős, P., Herzog, F., and Piranian, G., *Metric properties of polynomials*, J. Analyse Math. (1958), 125-148 — same paper as #114

Plan source:

- `~/.claude/plans/yes-wondrous-blum.md` section "Intractability test program — 2026-05-15 (POC framing)" Probe 3 (lines 318-326)
- `~/.claude/plans/yes-wondrous-blum-agent-a4ae5ad039883b7ba.md` (per-candidate execution method)
