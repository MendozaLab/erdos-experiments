# Cross-Problem WS-01 Generalization Outcome — Erdős #1038 Cold Test

Date: 2026-05-14
Plan track: Track B4 of `~/.claude/plans/yes-wondrous-blum.md` (Tao-gap-closing plan)
Claim ceiling: **Internal diagnostic only.** This document records the outcome of a Python-level cold-test of the cross-problem WS-01 applicability probe (Track B2) on a synthetically-constructed #1038 boundary slice (Track B3). It is NOT a proof of #1038, NOT a Lean statement, NOT a public artifact, and NOT a claim that "WS-01 generalizes." The diagnostic is limited to `mpmath.iv` at `dps=30` on N=8 perturbed Chebyshev-like roots; the B1 structural-fit verdict's 1e-8 interval-enclosure precondition is NOT met by this run.

## Background

The Tao-gap-closing plan asks whether the wall-separation rewrite WS-01-CENTER-STRIP-CANCELLATION (the Branch-Centered Moving-Frame Collar reduction developed in Erdős #114, CELL-02-03, owner family 4488:3) is portfolio-canonical (re-usable across multiple Erdős problems) or #114-specific (a one-trick reduction). Tracks B1 (structural fit), B2 (problem-agnostic probe), B3 (kernel + cold-test input), and B4 (this document) form a four-step cold-test loop on a target problem outside #114's home cell. The target chosen for this loop is Erdős #1038 — Tao's December 2025 sublevel-set problem on monic polynomials with real roots in [-1, 1].

Canonical #1038 statement (per B1, citation-locked to DeepMind formal-conjectures `Erdos1038.parts.i`): the infimum of `|{x ∈ ℝ : |f(x)| < 1}|` over nonconstant monic polynomials `f` such that all roots are real and lie in `[-1, 1]`. Known bounds: `2^(4/3) − 1 ≈ 1.5198 ≤ inf < 1.835`; conjectured `11/6 ≈ 1.8333…`; best local Model-Mayhem upper-bound witness `M ≈ 1.83639` at N=200 (not interval-certified).

B1's verdict was `PARTIAL_FIT_PROCEED_WITH_CAVEAT`: all three #114 features (zero-curve, wall-separation gate, branch-point IVT) port mechanically onto #1038, with the caveat that the failure mode in #1038 is parametric (boundary components merge/split in `(w, δ)` parameter space) rather than purely spatial (sign flips in `x`). B2 produced a problem-agnostic probe `cross_problem_ws01_applicability_probe.py` that consumes a boundary-slice JSON in #114's `remaining_unresolved[*]` schema and returns one of four verdicts: `WS01_APPLICABLE`, `WS01_PARTIAL`, `WS01_INAPPLICABLE_QUANTITATIVE`, `WS01_INAPPLICABLE_STRUCTURAL`. The probe's regression test on #114's 16-box CELL-02-03 owner-family-4488:3 reference run passes 16/16 at conservative 2R, matching the 2026-05-09 reference run. Track B3 is the construction of a minimal interval-arithmetic V_μ kernel for #1038 and the generation of a small set of candidate failure boxes through that kernel.

## What was done

A minimal `mpmath.iv` interval kernel for the #1038 boundary functional was built at `Erdos1038/build_1038_boundary_slice_cold_test.py`. The kernel evaluates `V_μ(x) = log|p(x)|`, `V_μ'(x) = Σ 1/(x − r_i)`, `V_μ''(x) = −Σ 1/(x − r_i)²`, and `V_μ'''(x) = 2·Σ 1/(x − r_i)³` over interval inputs at `dps=30`, with explicit pole detection (returns `None` if the interval contains a root, so callers can reject the box). The base polynomial is the N=8 Chebyshev-like configuration with roots `cos(kπ/9)` for `k = 1, …, 8`. Perturbations are parameterized by `(w, δ)` on the grid `w ∈ {0.10, 0.13, 0.16, 0.19}` × `δ ∈ {0.10, 0.15, 0.20}` (12 candidate parameter pairs around B1's recommended center `(0.174, 0.200)`). For each pair, the leftmost `k = ⌈w·N⌉` roots are shifted rightward by `δ`, the rightmost real zero of `log|p(x)|` outside the root support is located by float bisection at `mp.dps=50`, and a box of half-width `R = 1e-4` is constructed around that boundary point. Each box's failure-record schema fields are then filled from interval evaluations:

- `f_interval` = interval enclosure of `log|p|` over the box
- `inequality.center_fn_abs_lower` = lower bound on `|V_μ'(x_b)|` at the box center
- `inequality.tangent_radius_upper = R`
- `inequality.f_tt_abs_upper` = upper bound on `|V_μ''|` over the box
- `inequality.third_directional_upper` = upper bound on `|V_μ'''|` over the box
- `inequality.normal_remainder_bound = (T3/6) · R³` (interval Taylor remainder upper bound)
- `inequality.wall_lhs = sup_box |F| + R_n`
- `inequality.wall_rhs = R · lower(|V_μ'|)` with safety factor `S = R`
- `inequality.center_strip_f_abs_upper = sup_box |F|`
- `reason = "third_order_wall_separation_failed"` (matching #114's tag verbatim)
- `source_ownership_key = f"w={w:.3f}/d={delta:.3f}"` (the (w, δ) pair as the per-box stand-in for #114's owner-family key)

The 12 boxes were written to `Erdos1038/proof_path/EXP-MATH-ERDOS-1038-WS01-COLD-TEST-20260514-01/INPUT_boundary_slice_1038_2026-05-14.json` and fed to the B2 probe with `--problem-id 1038`. The probe wrote the standard three-file contract (`_RESULTS.json` + `_REPORT.md` + `_RESULTS.sha256`) into the same artifact folder.

## Result

The B2 probe returned `WS01_APPLICABLE` with the following aggregate numbers (per `EXP-MATH-ERDOS-1038-WS01-APPLICABILITY-PROBE-20260514-01_RESULTS.json`):

- Failure boxes processed: **12**
- Branch-point validation rate: **1.0** (12/12 boxes: `f_interval` straddles zero AND `center_fn_abs_lower > 0`)
- Closure rate at conservative 2R: **1.0** (12/12 close `new_LHS_2R < wall_rhs`)
- Closure rate at tight R (cumulative): **1.0**
- Structural-difference rate: **0.0** (no branch-point validation failures, no ill-defined gradients)
- Median required factor (RHS/LHS) before WS-01: **~0.9999** (the original wall test sits right at threshold)
- Median required factor (RHS/LHS) after WS-01 at 2R: **~2503.7** (WS-01 buys ~2500× margin over the original test)
- Still failing: 0
- Branch point not found: 0
- Branch point not well-defined: 0

The per-box outcomes table in the probe's REPORT.md shows `wall_LHS_before ≈ 1.0e-3 ≈ wall_rhs` (the original wall test is at threshold because in 1D the supremum of `|F|` over a box of half-width R scales as `R · |F'|`, exactly the RHS) while `wall_LHS_after_2R ≈ 4e-7` after the WS-01 rewrite eliminates the linear-in-R term. The post-WS-01 LHS is dominated by `0.5 · |F''|_upper · (2R)² ≈ 0.5 · 20 · 4e-8 ≈ 4e-7`, three orders of magnitude below the wall RHS.

## What generalized

Three pieces of the #114 architecture ran cold on #1038 without modification:

The probe's classification machinery (`classify_box` in `cross_problem_ws01_applicability_probe.py`) accepted the synthesized #1038 boxes as-is and produced sound outputs. The schema validation passed for all 12 boxes — every required field (`f_interval`, `inequality.center_fn_abs_lower`, `tangent_radius_upper`, `f_tt_abs_upper`, `third_directional_upper`, `normal_remainder_bound`, `wall_lhs`, `wall_rhs`) maps directly onto an interval-arithmetic computation against `V_μ` and its first three derivatives.

The branch-point validation logic — `f_interval` straddles zero AND `center_fn_abs_lower > 0` — fires cleanly. In #1038's 1D setting these correspond respectively to the interval IVT for `V_μ` having a zero in the box (the real boundary point `x_b`) and the Cauchy-kernel derivative `V_μ'` having a nonzero lower bound at `x_b` (which B1 numerically observed to be O(0.1–1) on recovered configurations and which we measured here at O(10) on the perturbed N=8 base).

The WS-01 quantitative gain — eliminating the linear-in-R term `|F_t(box-midpoint)| · R` in the wall LHS by re-anchoring at a validated zero of `F` — is exactly what the synthesized #1038 boxes show. The before-vs-after numbers (`~1e-3` → `~4e-7` for the LHS, against an unchanged `~1e-3` RHS) recover the same qualitative story as #114's CELL-02-03 4488:3 reference run.

## What didn't generalize (or required adaptation)

Three pieces required either honest re-interpretation or were not exercised by this cold test:

The 1D-vs-2D geometric frame is the first adaptation. In #114 the "tangent" and "normal" directions are genuine 2D unit vectors in coefficient × root-location space, and the wall test bounds `F` separately in those two directions. In #1038 the boundary is 0-dimensional in `x` for fixed `μ`, so the entire `(t, n)` frame collapses: there is no spatial tangent direction. The kernel reports sentinel values (`tangent = (1, 0)`, `normal = (0, 1)`) and the 1D Taylor argument becomes degenerate — `wall_LHS_before` is just `sup|F| + R_n` with no genuine `|F_t|`-vs-`|F_n|` distinction. The probe's algorithm still produces a valid closure verdict because the WS-01 trick (subtracting the linear-in-box-radius term that vanishes at a zero of `F`) is identical in 1D and 2D; only the geometric interpretation differs.

The parametric merge/split failure mode flagged by B1's Feature 2 was NOT exercised by this cold test. The (w, δ) grid is a small neighborhood that the kernel handles without any boundary components merging, so the probe never had to face the "parametric wall-separation gate over a 2D (w, δ) box" that B1 specified as the real B2 main gate. The 12 boxes here are 12 independent 1D boundary-slice tests at 12 different polynomial configurations, not a single 2D parametric gate. The cold test therefore validates that WS-01's *algebraic* form (the Taylor rewrite of the wall LHS) ports cleanly, but it does NOT validate that the *parametric* version of the gate ports — that remains a stretch target requiring an inari Rust kernel and the full 2D interval evaluation.

The 1e-8 enclosure precondition specified by B1 is NOT met. We ran at `mpmath.iv dps=30` with `R = 1e-4`, which gives enclosures on the order of `1e-3` for `F`, `1e-2` for `F'`, and `1e-1` for `F''`. These are far wider than the B1 stretch target. The diagnostic is honest about being a smaller-scale validity check, not a certified port. Additionally, three of the 12 (w, δ) pairs (rows 3, 6, 9 in the per-box outcomes table) produced identical boundary points because `⌈w·8⌉ = 2` for `w ∈ {0.13, 0.16, 0.19}` — the discrete shift-count function maps the lower end of our grid onto a single perturbed configuration. So the effective diversity of the 12 boxes is more like 6 (three (w, δ) sweep widths × two δ values that map to distinct k_shift counts). This is a thin-diversity artifact of the synthetic perturbation, not a bug, and it is recorded honestly in the per-box numbers.

## Implication for the Atlas / Collider claim

On this evidence, WS-01 looks **architecture-portable but not yet portfolio-canonical**. The probe's classification algorithm is genuinely problem-agnostic (it ran cold on #1038 with no #114-specific tweaks), and the branch-point validation + Taylor-rewrite logic produces sound, meaningful outputs on a different problem's boundary functional. That is a real positive datum for the Atlas / Collider thesis that wall-separation reductions are reusable across the Erdős corpus.

But "architecture-portable" is a much weaker claim than "portfolio-canonical." A portfolio-canonical reduction would survive the parametric merge/split failure mode that B1 flagged as #1038's distinguishing feature; this cold test did not exercise that failure mode. A portfolio-canonical reduction would also clear the B1 stretch precondition of 1e-8 enclosure on the recovered N=120 / N=200 configurations; this cold test ran at N=8 with `dps=30` and box half-width `1e-4`, three or four orders of magnitude shy of the precondition.

The honest read is: WS-01 is `#114 + plausible-on-#1038-cold-test`. To upgrade to portfolio-canonical we need (a) the inari Rust kernel for `V_μ` and derivatives at the recovered N=200 boundary points reaching 1e-8 enclosure, (b) the genuine 2D parametric gate over `(w, δ)` boxes, and (c) a non-extremal control comparison (uniform or arcsine measure on `[-1, 1]`) that visibly fails or shows weaker margins than the extremal. Without all three, the cold-test verdict here is signal, not conclusion.

## Honest scope

This is a Python-level diagnostic at modest precision. It is not a proof of #1038, not a Lean statement, and not a public-facing artifact. The claim ceiling explicitly bans calling this run "WS-01 generalizes" without the three stretch targets above. The kernel did run, the probe did produce a sound verdict, and the verdict was `WS01_APPLICABLE` on 12/12 boxes at conservative 2R — but the underlying boxes are synthesized from a perturbed N=8 base polynomial in a 1D-collapsed schema, not a full 2D parametric gate on a recovered Model-Mayhem N=200 witness. The right summary is "the architecture's algebraic core ports; the geometric and parametric scale-up does not yet."

## Artifacts

- B3 kernel + generator: `Math/erdos-experiments/Erdos1038/build_1038_boundary_slice_cold_test.py`
- Generated boundary-slice JSON: `Math/erdos-experiments/Erdos1038/proof_path/EXP-MATH-ERDOS-1038-WS01-COLD-TEST-20260514-01/INPUT_boundary_slice_1038_2026-05-14.json`
- Probe output (three-file contract): `Math/erdos-experiments/Erdos1038/proof_path/EXP-MATH-ERDOS-1038-WS01-COLD-TEST-20260514-01/EXP-MATH-ERDOS-1038-WS01-APPLICABILITY-PROBE-20260514-01_RESULTS.json` + `_REPORT.md` + `_RESULTS.sha256`
- B1 prior work: `Math/erdos-experiments/Erdos1038/STRUCTURAL_FIT_CHECK_2026-05-14.md`
- B2 probe (read-only here): `Math/erdos-experiments/cross_problem_ws01_applicability_probe.py`

## Sources cited

- B1 verdict and feature-by-feature structural fit check (read-only this session)
- B2 probe script, including its regression test against #114 CELL-02-03 owner-family-4488:3 reference run (16/16 at conservative 2R, verdict `WS01_APPLICABLE`; regression passed in this session before the #1038 cold test ran)
- DeepMind formal-conjectures `FormalConjectures/ErdosProblems/1038.lean` (canonical statement source, cited via B1)
- `mpmath` 1.3.0 interval arithmetic at `dps=30` (the kernel substrate)
