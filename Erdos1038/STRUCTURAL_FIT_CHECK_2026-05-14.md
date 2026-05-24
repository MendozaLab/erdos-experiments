# Structural Fit Check: WS-01-CENTER-STRIP-CANCELLATION (EHP114) applied to Erdős #1038

Date: 2026-05-14
Plan track: Track B1 of `~/.claude/plans/yes-wondrous-blum.md` (Tao-gap-closing plan)
Claim ceiling: **Internal structural-fit diagnostic only.** Not a proof, not a Lean statement, not a solution claim, not a public-facing artifact. This is an in-portfolio decision document about whether to invest B2 + B3 cold-test effort on #1038 or pivot to an alternative target.

## Problem statement (canonical, verbatim)

From `Math/formal-conjectures/FormalConjectures/ErdosProblems/1038.lean` lines 34-41 (DeepMind formal-conjectures, citation-locked to Tao25 blog post and erdosproblems.com/1038):

> What is the infimum of `|{x ∈ ℝ : |f x| < 1}|` over all nonconstant monic polynomials `f` such that all of its roots are real and contained in `[-1,1]`?

Companion statements in the same file:

- Part (ii): `sup = 2·2^(1/2)` (proved in Tao25)
- `inf < 1.835` (research solved upper bound)
- `inf ≥ 2^(4/3) - 1` (research solved lower bound, ≈ 1.5198, from Borwein-Erdélyi-Kós 1999)

Conjectured infimum from current local Model-Mayhem work: `11/6 ≈ 1.8333…`, with best computational upper-bound witness `M ≈ 1.83639` from a recovered N=200 atom-plus-cloud configuration (`ERDOS1038_CURRENT_STATUS_2026-05-08.md` lines 20-22). The decimal `1.83639` is reproducible but not interval-certified.

**Discrepancy with plan's rough description:** the plan paraphrased #1038 as "for monic polynomial `p(z)` of degree n with all roots in `[-1, 1]`, infimum of `m(p) := measure({x ∈ [-1, 1] : |p(x)| ≤ 1})`." The canonical statement differs on two points: (1) `x` ranges over all of ℝ, not `[-1,1]`; (2) the inequality is strict `< 1`, not `≤ 1`. The canonical statement is used throughout this document.

## Parameter space

The natural parameter space is the space of probability measures on `[-1,1]`, via the identification

```text
root measure  μ = (1/n) Σ δ_{x_i}  on [-1,1]
log potential V_μ(x) = ∫ ln|x - t| dμ(t)
sublevel set  Ω_μ = {x ∈ ℝ : V_μ(x) < 0}  (since |f(x)| < 1 ⇔ V_μ(x) < 0/n after monicness)
objective     M(μ) = |Ω_μ|  (Lebesgue measure)
```

(See `Math/Math-Problems/Erdos-Standard/erdos1038-fast/euler_lagrange_analysis.py` lines 1-40 and `MORPHISM_BRIDGE_WHITEPAPER_CORE.md` § 2.)

For a fixed atom-plus-cloud ansatz μ = (1-w)·δ_{+1} + w·μ_eq([-1, -1+δ]) the effective finite-dimensional parameter space is **`(w, δ) ∈ (0, 1) × (0, 2)`**, with the sublevel boundary points `{a_i, b_i}` and the auxiliary internal turning point `x_0` (between primary and secondary wells) being functions of `(w, δ)`. The full infinite-dimensional EL problem is on the space of probability measures on `[-1,1]` (Saff-Totik 1997 weighted-equilibrium framework).

For comparison, the EHP114 parameter space is `(complex polynomial coefficients) → real dimension 2n − 5 after symmetry`, with the wall-separation gate operating on individual hard cells in coefficient × root-location space (see `Math/erdos-experiments/Erdos114/EHP114_REGULAR_SLICE_WALL_TAYLOR_TARGET_2026-05-06.md`).

## Feature 1 — Zero-curve / level-set boundary

**Verdict: YES_STRUCTURAL_FIT**

The sublevel set boundary `∂Ω_μ = {x ∈ ℝ : V_μ(x) = 0}` is the direct analog of EHP114's lemniscate `|p(z)| = 1`. In coordinates, define

```text
F(x; μ) = V_μ(x) - 0 = ∫ ln|x - t| dμ(t)
```

so the zero-curve is `{F = 0}`. For a finite atom-plus-cloud configuration this is a finite set of real boundary points `{a_1 < b_1 < a_2 < b_2 < …}` (each component of Ω_μ is an open interval `(a_i, b_i)`). For the continuous-limit atom-plus-cloud ansatz `μ* = (1-w)·δ_{+1} + w·μ_eq` the boundary collapses to two intervals (a primary well anchored at +1 of width ≈ 1.836 and a tiny secondary well near -1 of width ≈ 0.035; `ERDOS1038_PHYSICS_MORPHISMS_PREPRINT.md` line 95).

Note one important degeneracy: in #1038 the zero-set is 0-dimensional in `x` (finitely many real roots of `V_μ`) for each fixed μ, whereas in #114 the zero-curve is 1-dimensional in `z` (a real algebraic curve in ℂ). The analog of the 1-dimensional zero-curve is the graph `{(μ, x) : F(x; μ) = 0}` in the product space (parameter × spatial). This is 1-dimensional along the spatial direction when μ varies in a 1-parameter family. The structural analog holds; the codimension story is the same as #114 if we treat μ as a deformation parameter rather than the spatial variable.

## Feature 2 — Wall-separation gate / coercivity condition

**Verdict: PARTIAL_STRUCTURAL_FIT**

The EHP114 wall-separation gate requires the zero-curve to split into two sheets across the normal direction with no critical points (`F_t = 0`) inside a thin collar at `s = 0`. The #1038 analog is the requirement that `V_μ'(x)` does not vanish in a collar around each boundary point `a_i, b_i`, so that the boundary moves smoothly under perturbation of μ and components do not merge or split.

The good news: at the boundary points of a well-conditioned atom-plus-cloud configuration, `V_μ'(a_i) < 0` (V increasing into the sublevel interval from outside on the left) and `V_μ'(b_i) > 0` (V decreasing into the interval from outside on the right) — these signs are required for the EL derivation and are observed numerically (`euler_lagrange_analysis.py` lines 21-23, 180-186). This gives bounded-away-from-zero `|V_μ'(x_b)|` at the boundary, which is the direct analog of the wall-separation condition.

The bad news, and why this is `PARTIAL` not `YES`: the Model-Mayhem preprint explicitly flags **non-smoothness of the sublevel-measure functional from boundary components merging/splitting under perturbation** (`ERDOS1038_PHYSICS_MORPHISMS_PREPRINT.md` line 71: "Direct L-BFGS-B re-optimization at these N converges to a boundary minimum … due to non-smoothness of the sublevel measure functional (boundary components merge/split under perturbation, making finite-difference gradient estimates unreliable).") This is the same failure mode as EHP114's wall-separation failure — `F` intervals straddling zero inside the collar — but it is parametric (along μ) rather than purely spatial (along x). A wall-separation gate in #1038 would need to certify that within a neighborhood of the extremal `(w*, δ*)` the secondary well does not collapse onto the primary well and that no new boundary points spawn. This is more delicate than the EHP114 gate because the failure mode is parametric merge/split, not just sign-flipping of `Fx`.

Additionally, the recent variational-boundary-kernel robustness audit (`ERDOS1038_CURRENT_STATUS_2026-05-08.md` lines 102-119) showed that the sharp boundary kernel separates hard controls (60/72) but is brittle on near-perturbations (39/108), while the smoothed kernel stabilizes near-perturbations (324/324) but loses hard-control rejection (9/216). This brittleness-vs-bluntness tradeoff is exactly the wall-separation gate problem that EHP114's center-strip cancellation was designed to handle, which is a positive sign — but the gate has not yet been ported. Hence `PARTIAL`, with the structural reason for upgrading to `YES` being that the failure mode is the same family as EHP114's, even if the parametric/spatial direction differs.

## Feature 3 — Branch-point validation via interval IVT + bounded-away gradient

**Verdict: PARTIAL_STRUCTURAL_FIT**

The EHP114 gate requires every failure box to have `F_iv` straddling zero (interval IVT delivers a zero) AND `|grad F|` bounded away from zero at the box center (zero-curve is a smooth 1-manifold; tangent direction `t(u)` is well-defined).

For #1038 the analog is, for each candidate boundary point `x_b` of `V_μ`:

1. **Interval IVT for the zero-curve in `x`.** `V_μ(x)` is real-analytic in `x` away from the support of μ. Standard interval-arithmetic evaluation of `V_μ([x_lo, x_hi])` is feasible — the log kernel `ln|x - t|` has a known monotone interval rule away from `x = t`, with controlled treatment of the near-pole region. Confirming `V_μ([x_lo, x_hi])` straddles zero in a small box around each `x_b` is a tractable interval computation. This piece is fully checkable.
2. **Bounded-away `|V_μ'(x_b)|`.** Numerically observed `|V_μ'(a_i)|, |V_μ'(b_i)|` are O(0.1 – 1) on the recovered configurations (`euler_lagrange_analysis.py` § Section 2). Interval-arithmetic certification of a positive lower bound on `|V_μ'|` in a collar is again tractable — the derivative `V_μ'(x) = ∫ (x - t)^{-1} dμ(t)` has a Cauchy-kernel interval-arithmetic recipe. This piece is also checkable.

So the spatial-direction (boundary-point smoothness in `x`) version of Feature 3 is **YES**.

The reason this is `PARTIAL` rather than `YES`: the EHP114 gate is really validating branch points of the **zero-curve as a function of the parameter** (where the curve becomes singular and tangent direction breaks down). For #1038 the corresponding objects are the **critical (w, δ) where components merge or split** — and these are exactly the loci where `V_μ'(x_b) → 0` as `x_b` slides into a neighboring boundary point. The gate is checkable for fixed μ, but the question of whether the **set of admissible μ** has bounded-away non-degenerate boundary structure across the full neighborhood of the extremal is the parametric question flagged in Feature 2. The same gate machinery can be applied (interval IVT in `(w, δ)` plus bounded-away `|V_μ'|`), but the certification is over a 2D parametric box rather than a 1D spatial box, which is a quantitative scale-up rather than a qualitative obstruction.

## Aggregate verdict

**`PARTIAL_FIT_PROCEED_WITH_CAVEAT`**

All three features map cleanly onto the #1038 setup: there is a well-defined `F = V_μ` zero-curve (Feature 1: YES), a wall-separation gate analog with a known but tractable failure mode (Feature 2: PARTIAL, identified failure family matches EHP114's), and interval-IVT-plus-gradient validation is feasible in both spatial and parametric directions (Feature 3: PARTIAL, gate machinery ports, scale increases). Crucially, no feature is `NO_STRUCTURAL_FIT`. The problem is structurally the right shape for WS-01-CENTER-STRIP-CANCELLATION to apply.

The caveats are real but mechanical, not architectural:

1. The wall-separation failure mode in #1038 is parametric (μ-direction merge/split) more than spatial (x-direction sign flip). The center-strip cancellation will need to be redesigned to operate over a `(w, δ)`-box rather than a spatial collar.
2. The current Model-Mayhem variational-boundary-kernel attempts already exhibit the brittleness-vs-bluntness tradeoff that the EHP114 gate addresses. The port should improve on the current sharp-kernel approach.
3. Local upgrade pre-conditions before B2/B3 should include: (a) interval-arithmetic implementation of `V_μ` and `V_μ'` (this does not exist locally yet; the current implementations use float64 + finite differences); (b) a clean statement of the parametric wall-separation gate in `(w, δ)` coordinates; (c) interval root-enclosure for boundary points `{a_i, b_i}` rather than `brentq` float root-finding.

## Recommended next step

**Proceed to B2 + B3 cold test on #1038, with three pre-condition gates.**

Sequence:

1. **B2 pre-condition (interval V_μ kernel).** Implement an interval-arithmetic version of `V_μ` and `V_μ'` for atom-plus-cloud configurations. Use `mpmath` or a Rust `inari`-backed kernel (the existing `erdos1038_fast` Rust evaluator already uses `inari` for the sublevel-measure scan — extend to the EL gradient kernel). Pass condition: interval enclosure of `V_μ(x_b)` straddles zero by less than `1e-8` width on the recovered N=120 / N=200 configurations.
2. **B2 main (parametric wall-separation gate in (w, δ)).** Build the interval-IVT + bounded-away-`|V_μ'|` gate over a small 2D box around `(w*, δ*) ≈ (0.174, 0.200)`. Pass condition: gate certifies smoothness of the boundary in the entire box, with explicit margins on `|V_μ'(x_b)|` and on the gap between the primary and secondary wells.
3. **B3 cold test.** Compare the parametric wall-separation gate output for #1038 against an analogous gate applied to a non-extremal control measure (e.g., uniform measure or arcsine on `[-1, 1]`). The control should fail or produce visibly weaker margins; the extremal should pass with a clear margin. This is the cold-test analog of EHP114's hard-cell `(6,4)` regular-slice diagnostic.

If any of the three gates fails, pivot to an alternative (see below).

## Alternative candidate (only if B2/B3 fails)

**Primary fallback: Erdős #20 (sunflower / Szemerédi).**

Rationale (one sentence): #20 has a discrete combinatorial extremal structure with a Razborov-style approximation lemma that admits an explicit "wall" between extremal and non-extremal configurations, providing a different (and possibly easier) shape of wall-separation gate than the continuum #1038 problem.

**Secondary fallback: Erdős #233 (KvN / dynamical-systems variant).**

Rationale: #233 is scoped to the H² Koopman-von Neumann bridge that the portfolio already has live infrastructure for (RPQEC / Model-Mayhem KvN runners), and the wall-separation analog there is operator-spectral rather than measure-theoretic, which avoids the parametric merge/split failure mode of #1038 entirely.

The pivot should be triggered only if the B2 pre-condition (interval V_μ kernel) cannot reach `1e-8` enclosure width on the recovered configurations, since that is the cheapest of the three gates and a hard failure there indicates the interval-arithmetic machinery cannot carry the rest of the gate.

## Sources cited

Local files read for this fit check:

- `Math/formal-conjectures/FormalConjectures/ErdosProblems/1038.lean` (lines 1-70) — canonical Lean problem statement
- `Math/erdosatlas-workbench/data_local/deepmind_erdos_statements.json` (entry `problem_id = 1038`) — DeepMind formal-conjectures embedding text
- `Math/Math-Problems/Erdos-Standard/erdos1038-fast/ERDOS1038_CURRENT_STATUS_2026-05-08.md` (lines 1-146) — current status, N=200 recovery, variational-kernel robustness audit
- `Math/Math-Problems/Erdos-Standard/erdos1038-fast/euler_lagrange_analysis.py` (lines 1-391) — explicit EL derivation, gradient function G(t), boundary-derivative signs
- `Math/Math-Problems/Erdos-Standard/erdos1038-fast/MORPHISM_BRIDGE_WHITEPAPER_CORE.md` (lines 1-100) — logarithmic potential reformulation, Thomson/Riesz morphism, three proof bridges
- `Math/Math-Problems/Erdos-Standard/erdos1038-fast/ERDOS1038_PHYSICS_MORPHISMS_PREPRINT.md` (lines 71, 95, 137, 163) — non-smoothness of sublevel functional, two-component structure, turning-point condition, Saff-Totik regime
- `Math/Math-Problems/Erdos-Standard/erdos1038-fast/analyze_structure.py` (lines 1-80) — N-point optimizer structure analysis
- `Math/erdos-experiments/Erdos114/EHP114_REGULAR_SLICE_WALL_TAYLOR_TARGET_2026-05-06.md` (lines 1-97) — EHP114 wall-separation gate failure pattern (used for comparison; not part of #1038 inputs)
- `Math/erdos-experiments/results/erdos-1038/tao_ansatz_eps0.0005_inari_RESULTS.json` (Tao ansatz inari-backed numerical results, listing only; not opened for content here)

External canonical reference (not refetched this session, locked via Lean file):

- erdosproblems.com/1038 (Erdős's original problem entry)
- Tao, T. *Sublevel Sets of Logarithmic Potentials*, Terry Tao's Blog, Dec. 2025 (PDF at `terrytao.wordpress.com/wp-content/uploads/2025/12/erdos-1038-1.pdf`)

Plan source:

- `~/.claude/plans/yes-wondrous-blum.md` (Track B1 of the Tao-gap-closing plan; read by reference only — not directly opened this session because the plan text is in the user prompt).
