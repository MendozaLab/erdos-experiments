# Proof Recipes

## 2026-05-05 - EHP114 local mixed-remainder closure scaffold

Compiled source: `erdos-experiments/Erdos30/scratch/Ehp114LocalMixedRemainderScratch.lean`

### Recipe: isolate the analytic bottleneck behind component-budget algebra

Use this when a hard local-stability proof decomposes into several certified
or externally targeted estimates and one remaining mixed term.

Pattern:

1. Define the local objects as abstract interfaces first: radial deficit, shape
   quadratic, mixed remainder, admissibility, and local cone.
2. Package each estimate as a named proposition rather than hiding it in the
   final theorem statement.
3. Prove the final reserve inequality from those named propositions by pure
   linear arithmetic.
4. Leave the real mathematical bottleneck as one named target proposition.

In the current #114 lane this appears as `RadialCertificate14`,
`ShapeConeCertificate14`, `MixedRemainderAbsorption14`, and
`local_deficit_reserve_from_components`. The key lesson is that the theorem
frontier is not the algebraic splice; it is the uniform mixed-remainder
absorption theorem on the local admissible cone.

### Mathlib / repo pitfall

Lean 4.27 in this project accepted explicit `axiom` declarations for abstract
analytic objects, while `constant` declarations in the same position failed to
parse. For scratch theorem-target scaffolds in this repo, use `axiom` for
unimplemented analytic interfaces and document the claim ceiling in the file
header.

## 2026-04-20 - Erdos30_BFR support lemmas

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_BFR.lean`

### Recipe: range bound via injective image plus `Finset.card_le_card`

Use this when a Sidon map is already known to be injective and the remaining job is only to count how many target values are possible.

Pattern:

1. Build the image finset explicitly.
2. Prove it sits inside `Finset.range bound`.
3. Finish with `Finset.card_le_card` on the subset proof.

In the current BFR batch this appears as `distinctSums_subset_range` and `erdos_turan_counting_bound`. The useful split is to keep the subset lemma separate from the final cardinality inequality so the same subset proof can feed more than one counting argument.

### Recipe: shifted-family intersection at most one

For additive-combinatorics arguments, define a translated copy
`shifted A t := A.image (fun a => a + t)` and count collisions between two translates through a Sidon difference lemma.

Pattern:

1. Prove `card_shifted` by `Finset.card_image_of_injective`.
2. Prove each shift stays inside a concrete range (`shifted_subset_range`).
3. For `card (shifted A t ∩ shifted A s) <= 1`, convert membership in the intersection into
   `a1 + t = a2 + s`, rearrange to `a1 - a2 = s - t`, and invoke the imported Sidon distinct-differences lemma.

In the current BFR batch this appears as `shifted`, `card_shifted`, `shifted_subset_range`, and `shifted_inter_card_le_one`.

### Mathlib / repo pitfall

Do not import `Erdos30_Complete` just to reach a Sidon difference lemma. In this repo that import drags a stale dependency chain. Import `Erdos30_Lindstrom` directly and reuse `sidon_distinct_differences`.

## 2026-04-20 - Dense finite rigidity scratch targets

Compiled source: `erdos-experiments/Erdos30/scratch/Erdos30_IntervalOccupancyTarget.lean`

### Recipe: convert an initial-segment counting statement into a sorted-prefix identity

Use this when the mathematical claim is "the elements of `A` below the `i`-th order statistic are exactly the first `i+1` sorted elements," and Lean needs that turned into a cardinality theorem on a filtered finset.

Pattern:

1. Set `l := A.sort (· ≤ ·)` and `a := orderedElement A i`.
2. Express the initial segment as a boolean filter on the list: `p x := decide (x ≤ a)`.
3. Prove `(l.filter p).toFinset = intervalSlice A 0 a` by extensional `simp`.
4. Compute the filter length by splitting `l` as `take (i+1) ++ drop (i+1)`.
5. Show the prefix contributes all its length using sortedness and `rel_get_of_le`.
6. Show the suffix contributes zero using strict sortedness plus `rel_of_mem_take_of_mem_drop`.
7. Convert `toFinset.card` back to list length via nodup.

In the current scratch target this is the proof of `ordered_prefix_card_target`. The reusable lesson is that `List.countP`, `take_append_drop`, and `sort_toFinset` are the clean bridge between order-statistic statements and finset cardinality statements.

### Recipe: isolate an honest algebraic obstruction with a correction-term split

Use this when a heuristic target is "almost" obtained from a proved counting theorem, but one residual term still needs separate control.

Pattern:

1. Package the proved error term as `q`.
2. Package the residual correction term as `r`.
3. Rewrite the target expression as `-q + r`.
4. Apply `Int.natAbs_add_le`.
5. Keep the correction term explicit instead of hiding it in an over-optimistic target.
6. Add a promotion lemma saying that any standalone bound on `r` upgrades the corrected statement to the desired final shape.

In the current scratch target this appears as `dense_sidon_ordered_element_with_correction` and `dense_sidon_ordered_element_of_correction_bound`. This is the right pattern whenever the proof architecture is sound but one geometric estimate is still missing.

### Recipe: transport an external ordered-element theorem to a prefix-count statement at cutpoints

Use this when the literature gives an asymptotic theorem for the `i`-th ordered element, while the local formalization wants a statement about prefix counts evaluated at `t = a_i`.

Pattern:

1. Keep the external theorem as an explicit axiom/interface with citation and honest scale.
2. Prove the exact local identity `|(A ∩ [0, a_i])| = i+1`.
3. Rewrite the literature theorem by replacing `i+1` with that exact prefix count.
4. Treat the result as a cutpoint consequence, not as a uniform discrepancy theorem on all prefixes.

In the current scratch target this appears as `dense_sidon_ordered_element_external` together with `ordered_prefix_card_target`, yielding `dense_sidon_prefix_cutpoint_external`. The lesson is that exact local combinatorics can still extract useful structure from an imported asymptotic theorem even when the stronger local discrepancy theorem is out of reach.

### Recipe: package nearby-prefix control into an index-free nonterminal theorem

Use this when the local combinatorics naturally produce a theorem "between consecutive ordered elements," but the downstream statement wants to quantify over `t` directly rather than over a bracketing index.

Pattern:

1. Prove an internal bracketing lemma: if `0 < s(t) < |A|`, then `t` lies between `a_{s(t)-1}` and `a_{s(t)}`.
2. Define the local index `i := s(t) - 1` as a `Fin A.card`.
3. State the successor-index equality explicitly as an equality of `Fin` terms. Do not expect `omega` to rewrite `⟨i.1 + 1, _⟩` to `⟨s(t), _⟩` inside a larger application on its own.
4. Handle the empty-prefix boundary `s(t) = 0` separately by showing `t < a_0`; the contradiction is that otherwise `a_0 ∈ A ∩ [0,t]`.
5. Combine the internal and empty-prefix constants with an additive envelope to get a theorem for all nonterminal prefixes `s(t) < |A|`.

In the current scratch target this appears as `dense_sidon_internal_prefix_external`, `dense_sidon_empty_prefix_external`, and `dense_sidon_nonterminal_prefix_external`. The key lesson is that one-step geometric control plus an honest boundary split is enough to recover a usable theorem for all nonterminal prefixes without overclaiming the terminal regime.

### Recipe: terminal-prefix closure via last-element control plus floor drift

Use this when a prefix statement is already proved for all nonterminal prefixes and the only missing case is when the prefix has swallowed all of `A`.

Pattern:

1. Let `iLast := ⟨|A|-1, _⟩` and use the external ordered-element theorem at `a_last`.
2. Prove combinatorially that if `|A ∩ [0,t]| = |A|`, then `a_last ≤ t`; otherwise `intervalSlice A 0 t` is a strict subset of `A`.
3. Bound the upper side by `t ≤ n`.
4. Convert the floor deficiency relation `|A| + L = Nat.sqrt n` into the real drift bound
   `n - |A|√n ≤ (L+1)√n`.
5. Package the result as a terminal-prefix theorem, then dominate the nonterminal theorem by the same drift term to get a single all-prefix statement.

In the current scratch target this appears as `dense_sidon_terminal_prefix_external` and `dense_sidon_positive_card_prefix_external`. The lesson is that the endpoint does not need a sharper external theorem; it needs the correct center and a separate floor-error accounting.

### Recipe: sum an external pointwise theorem before chasing closed forms

Use this when the literature gives a pointwise bound for each ordered element and the next natural structural consequence is a global mass-balance statement.

Pattern:

1. Rewrite the target mass difference as a sum of pointwise deviations using `Finset.sum_sub_distrib`.
2. Apply `Finset.abs_sum_le_sum_abs`.
3. Bound each summand by the same external error term with `Finset.sum_le_sum`.
4. Collapse the constant sum to `|A| * E` by `simp`.
5. Only afterwards decide whether it is worth simplifying the profile sum to a closed form.

In the current scratch target this appears as `dense_sidon_ordered_mass_external`. The key lesson is that the mathematically meaningful theorem is the summed deviation bound; the closed form for the profile sum is secondary bookkeeping and should not block the structural result.

### Recipe: bridge ordered-index sums back to the finset sum through `List.ofFn`

Use this when a theorem is naturally proved by summing over `i : Fin |A|`, but the user-facing statement should talk about `A.sum id`.

Pattern:

1. Prove locally that `List.ofFn (fun i => orderedElement A i) = orderedElements A` via `List.ofFn_getElem`.
2. Take `List.sum` of both sides.
3. Rewrite the left side with `List.ofFn_eq_map` and `Fin.sum_univ_eq_sum_range`.
4. Rewrite the right side from `orderedElements A = A.sort (· ≤ ·)` to `A.sum id` using `List.sum_toFinset` and `Finset.sort_nodup`.
5. Cast this local `ℕ`-valued theorem up to `ℝ` only when the external theorem needs it.

In the current #30 lane this appears as `sum_orderedElement_eq_sum` in the local file and `dense_sidon_finset_mass_external` in the scratch file. The lesson is that the `Fin`-indexed proof architecture and the finset-facing theorem can stay cleanly separated, with one small bridge lemma handling the translation.

### Recipe: make the affine profile sum explicit only after the structural theorem is done

Use this when the core theorem already controls
`|A.sum id - ∑ i, ((i+1) : ℝ) * s|` and you want the readable arithmetic center
without reopening the main proof.

Pattern:

1. Keep the structural theorem in the `Fin`-indexed profile form.
2. Prove a separate closed-form lemma with `Fin.sum_univ_eq_sum_range`.
3. Rewrite the natural-number sum by `Finset.sum_add_distrib` and `Finset.sum_range_id`.
4. Cast only once with `Nat.cast_sum`, then `simpa [Nat.cast_add]`.
5. Rewrite the main theorem with the closed-form lemma as a final corollary.

In the current #30 lane this appears as `dense_sidon_profile_sum_explicit` and `dense_sidon_finset_mass_explicit_external`. The lesson is that closed-form arithmetic is bookkeeping, not the backbone of the rigidity argument.

### Recipe: switch from floor-deficiency to symmetric `| |A| - sqrt(n) |` drift at the terminal boundary

Use this when an external ordered-element theorem is stated against
`max(0, sqrt(n) - |A|)`, but the prefix theorem you want must still make sense
for sets with `|A| > floor(sqrt(n))`.

Pattern:

1. Keep the ordered-element error term in the literature deficiency variable.
2. For the terminal prefix case, separate the global drift `n - |A| * sqrt(n)`
   from the local ordered-element error.
3. Rewrite `n` as `(sqrt(n))^2`.
4. Bound `sqrt(n) - |A|` by `|(A.card : ℝ) - sqrt(n)|`.
5. Package the terminal theorem with symmetric drift, then dominate the
   nonterminal theorem by adding one extra `sqrt(n)` step.

In the current #30 lane this appears as
`sidon_in_range_terminal_prefix_external` and the tightened
`sidon_in_range_positive_card_prefix_external`. The lesson is that the ordered
element theorem can stay at the published deficiency scale while the prefix
package uses a stronger, symmetric bookkeeping term at the endpoint, then
compresses the all-prefix wrapper to `max(| |A| - sqrt(n) |, 1) * sqrt(n)`.

### Recipe: absorb a set-dependent endpoint drift into an `n`-only fourth-root term on the super-floor corridor

Use this when a prefix theorem still carries a term like
`max(| |A| - sqrt(n) |, 1) * sqrt(n)`, but you are allowed to assume
`floor(sqrt(n)) ≤ |A|`.

Pattern:

1. Split the density gap into lower-side and upper-side control.
2. Use `floor(sqrt(n)) ≤ |A|` plus `Real.real_sqrt_le_nat_sqrt_succ` to show
   the lower-side deficit is at most `1`.
3. Use `lindstrom_bound` to cap the upper-side excess by `Nat.sqrt (Nat.sqrt n) + 1`.
4. Reassemble the symmetric gap as an absolute value bounded by
   `Nat.sqrt (Nat.sqrt n) + 1`.
5. Use the same super-floor hypothesis to bound the truncated literature
   deficiency `max(0, sqrt(n) - |A|)` by `1`.
6. Feed both estimates into the existing prefix theorem to replace the
   set-dependent endpoint term by a pure ambient correction
   `((Nat.sqrt (Nat.sqrt n) : ℝ) + 1) * sqrt(n)`.

In the current #30 lane this appears as
`realGapFromSqrt_le_fourthRoot_of_superfloor`,
`realDeficiencyFromSqrt_le_one_of_superfloor`, and
`sidon_in_range_superfloor_prefix_external`. The lesson is that once the set is
already known to live in the super-floor corridor, Lindström is strong enough
to turn the endpoint drift from geometric bookkeeping about `A` into a pure
fourth-root ambient term, which is the first honest absorption step toward a
no-drift prefix theorem.

### Recipe: collapse the super-floor fourth-root correction into the single `n^(7/8)` scale

Use this after the previous recipe, when the theorem already has the form

`prefix error ≤ ambient(n) + C * n^(7/8) + C * n^(3/4)`

with `ambient(n)` independent of the set.

Pattern:

1. Prove the ambient correction is itself `O(n^(7/8))` for `n ≥ 1`.
2. Do this with exponent monotonicity on base `n ≥ 1`, not with ad hoc case
   splits: `n^(1/2) ≤ n^(7/8)` and `n^(3/4) ≤ n^(7/8)`.
3. Rewrite the mixed square-root product using `Real.sqrt_eq_rpow` and
   `Real.rpow_add_of_nonneg`.
4. Absorb the remaining `C * n^(3/4)` term the same way.
5. Package the result with a new coarse constant, keeping the theorem honest
   about positivity assumptions on `n`.

In the current #30 lane this appears as
`superfloor_ambient_le_two_sevenEighths` and
`sidon_in_range_superfloor_prefix_coarse_external`. The lesson is that once the
endpoint drift has been made ambient-only, the final collapse to a single
`n^(7/8)` scale is mostly exponent bookkeeping. That is the first compiled
prefix theorem in the current lane with no explicit set-dependent drift and no
separate fourth-root correction term.

### Recipe: port ordered-element difference control from the dense corridor to the maximizer-friendly lane

Use this when a theorem for `DenseSidonAtScale` already controls
`a_j - a_i - (j-i) sqrt(n)` and you want the same statement for `SidonInRange`.

Pattern:

1. Reuse the external ordered-element theorem twice, once at `i` and once at `j`.
2. Package the common error term as a local `E`.
3. Rewrite `(j-i)` as `((j+1) - (i+1))` by `Nat.cast_sub`.
4. Subtract the two ordered-element approximations and use `abs_sub`.
5. Only after the index-difference theorem is in place, rewrite it at cutpoints
   using `ordered_prefix_card_target`.
6. On the super-floor corridor, absorb the remaining deficiency term exactly as
   in the prefix theorem: `sqrt(def) ≤ 1` and `n^(3/4) ≤ n^(7/8)`.

In the current #30 lane this appears as
`sidon_in_range_index_difference_external`,
`sidon_in_range_cutpoint_interval_external`, and
`sidon_in_range_superfloor_cutpoint_interval_coarse_external`. The lesson is
that the interval-rigidity layer should be built from cutpoints outward: first
get the two-index ordered displacement theorem, then convert it into a cutpoint
interval theorem, and only then try to say something about arbitrary `t`
between cutpoints.

### Recipe: collapse internal nearby-prefix control to the single `n^(7/8)` scale on the super-floor corridor

Use this when the file already has:

1. an index-free internal-prefix theorem of the form
   `sqrt(n) + C * n^(7/8) + C * sqrt(def) * n^(3/4)`, and
2. the super-floor bound `realDeficiencyFromSqrt A n ≤ 1`.

Pattern:

1. Stay in the genuinely internal regime `0 < s(t) < |A|`, so no terminal drift
   is needed.
2. Absorb the leading `sqrt(n)` step into `n^(7/8)` using `n ≥ 1`.
3. Bound `sqrt(def) ≤ 1`, then absorb the remaining `n^(3/4)` term into
   `n^(7/8)` by exponent monotonicity.
4. Package the result as an index-free theorem, not a theorem parameterized by
   a bracketing index.

In the current #30 lane this appears as
`sidon_in_range_superfloor_internal_prefix_coarse_external`. The lesson is that
once the theorem is restricted to internal prefixes, the endpoint bookkeeping
disappears and the super-floor corridor naturally supports a clean no-extra-drift
statement at the literature `n^(7/8)` scale.

### Recipe: isolate the super-floor boundary branches instead of hiding them in the global wrapper

Use this once the file already has:

1. a coarse super-floor all-prefix theorem,
2. a coarse internal-prefix theorem, and
3. the endpoint absorption ingredients `realGapFromSqrt_le_fourthRoot_of_superfloor`,
   `realDeficiencyFromSqrt_le_one_of_superfloor`, and
   `superfloor_ambient_le_two_sevenEighths`.

Pattern:

1. Prove the empty-prefix coarse theorem directly from the existing
   `sidon_in_range_empty_prefix_external` plus the same `sqrt(n)` and
   `n^(3/4)` absorption used in the internal theorem.
2. Package the empty and internal branches together into a clean nonterminal
   theorem by a `by_cases` split on `0 < |A ∩ [0,t]|`.
3. Prove the terminal coarse theorem separately by absorbing
   `realGapFromSqrt A n * sqrt(n)` into the ambient fourth-root term, then into
   the same single `n^(7/8)` scale.
4. Keep these theorems even if a stronger global wrapper already exists; the
   point is architectural honesty about where the bulk and endpoint control live.

In the current #30 lane this appears as
`sidon_in_range_superfloor_empty_prefix_coarse_external`,
`sidon_in_range_superfloor_nonterminal_prefix_coarse_external`, and
`sidon_in_range_superfloor_terminal_prefix_coarse_external`. The lesson is that
the super-floor prefix package is cleaner when the boundary cases are explicit
theorems rather than hidden inside one broad all-prefix statement.

### Recipe: collapse the super-floor mass theorem to the honest `n^(11/8)` envelope

Use this when the file already has:

1. the maximizer-friendly explicit mass theorem with error
   `|A| * (n^(7/8) + sqrt(def) * n^(3/4))`, and
2. the super-floor control `floor(sqrt(n)) ≤ |A|`.

Pattern:

1. Start from the explicit mass theorem, not from a fresh summation-by-parts
   argument. The structural work is already done there.
2. Absorb `sqrt(def) ≤ 1` exactly as in the coarse prefix theorems.
3. Use Lindström plus the super-floor hypothesis to bound `|A|` by a constant
   multiple of `sqrt(n)`.
4. Multiply the `n^(7/8)` term by `sqrt(n)` to get the natural mass scale
   `n^(11/8)`.
5. Bound the lower-order `sqrt(n) * n^(3/4)` term by the same `n^(11/8)`
   envelope.

In the current #30 lane this appears as
`sidon_in_range_superfloor_finset_mass_coarse_external`. The lesson is that the
mass theorem should be treated as an envelope theorem at `n^(11/8)` scale, not
as a prefix-style rigidity statement at `n^(7/8)` scale.

### Recipe: change the mass center after the envelope theorem, not before it

Use this when exact data says a density-adjusted center is better for mass, but
the proof architecture already naturally lands on the older `sqrt(n)`-centered
template.

Pattern:

1. Prove the coarse mass envelope first around the center you get for free from
   summing the ordered-element theorem.
2. Introduce the new target center only afterwards as a separate rewrite
   theorem.
3. Bound the shift between the two centers by factoring it into
   `(|A| + 1) * |sqrt(n) - |A|| * sqrt(n)`.
4. On the super-floor corridor, absorb `|sqrt(n) - |A||` into the same ambient
   fourth-root term already used in the prefix theorems.
5. Use the existing super-floor cardinality bound `|A| = O(sqrt(n))` to keep
   the whole shift at the same `n^(11/8)` scale.
6. Finish with one triangle inequality, so the center change is visibly a
   quantitative refinement of the same envelope theorem rather than a new proof
   architecture.

In the current #30 lane this appears as
`sidon_in_range_superfloor_finset_mass_density_adjusted_external`. The lesson
is that the honest place to introduce the density-adjusted center is after the
mass envelope has already been proved. The data says density adjustment helps
mass, but not prefix, so the theorem architecture should reflect that split
directly.

### Recipe: add extremal-surface vocabulary before proving compatibility

Use this when exact finite diagnostics suggest a phenomenon lives on the
maximal surface rather than across a generic dense layer.

Pattern:

1. Add a cardinality-maximal predicate as a `Prop`, not as a computable `h(n)`
   function. For #30 this is `IsMaximalSidonInRange A n`.
2. Add a slackened cardinality predicate separately. For #30 this is
   `NearExtremalSidonInRange A n δ`.
3. Prove only projection lemmas back to the existing base package first, such
   as `.sidonInRange`, and a zero-slack inclusion from maximal to near-extremal.
4. Name the observables as definitions before stating any optimizer theorem.
   For #30 these are `prefixResidualAfterGeneralDrift A n t` and
   `densityAdjustedMassDeviation A n`.
5. Keep optimizer-selection, Pareto-frontier, and compatibility claims out of
   this first patch. The point is to make the next theorem stateable without
   smuggling finite experimental conclusions into Lean.
6. Add only wrapper theorems that reuse already-compiled envelopes under the
   new surface predicate. For #30 this is
   `maximal_sidon_in_range_superfloor_prefix_mass_joint_envelope_external`.

In the current #30 lane this appears in
`Erdos30_IntervalOccupancyTarget.lean`. The lesson is that an exact-data
surface effect needs vocabulary before proof ambition: first define the surface
and observables, then decide whether any compatibility theorem is honest.

## 2026-04-28 — Transdimensional Painter axiom-discharge salvo (29 → 22)

Compiled source: `Math/Lean4/transdimensional-painter/TransdimensionalPainter/{Layer3a_SpectralUnitary,Layer4_TakensEmbedding,Layer5_LemniscateGeometry}.lean`

### Recipe: norm-preservation-to-eigenvalue-on-circle for unitary operators

Use this when proving `‖λ‖ = 1` for an eigenvalue of a Hilbert-space unitary operator from the inner-product preservation hypothesis.

Pattern:

1. Get `‖U v‖ = ‖v‖` from `⟪U v, U v⟫_ℂ = ⟪v, v⟫_ℂ` via `inner_self_eq_norm_sq_to_K` on both sides, equating, casting to ℝ via `exact_mod_cast`, and closing with `nlinarith [norm_nonneg ..]`.
2. Substitute the eigenvector equation `U v = λ • v` and use `norm_smul` to get `‖λ‖ * ‖v‖ = ‖v‖`.
3. Cancel `‖v‖ ≠ 0` (from `v ≠ 0` via `norm_ne_zero_iff`) using `mul_right_cancel₀`.

Mathlib API used: `inner_self_eq_norm_sq_to_K`, `norm_smul`, `norm_ne_zero_iff`, `mul_right_cancel₀`. The `inner_self_eq_norm_sq_to_K` rewrite is version-stable across Mathlib v4.20+. Lives at `eigenvalue_on_circle` in `Layer3a_SpectralUnitary.lean`.

### Recipe: distinct-eigenvalue orthogonality via complex-conjugate symmetry

Use this for `⟪v₁, v₂⟫ = 0` between unitary-operator eigenvectors with distinct eigenvalues, when the surrounding context already has `eigenvalue_on_circle` available.

Pattern:

1. Case-split on `v₁ = 0` and `v₂ = 0` separately (`simp` closes both).
2. From inner-product preservation `⟪U v₁, U v₂⟫ = ⟪v₁, v₂⟫`, rewrite via `hv₁`, `hv₂`, `inner_smul_left`, `inner_smul_right` to get `starRingEnd ℂ λ₁ * (λ₂ * ⟪v₁, v₂⟫) = ⟪v₁, v₂⟫`.
3. Factor as `(starRingEnd ℂ λ₁ * λ₂ - 1) * ⟪v₁, v₂⟫ = 0` using `linear_combination`.
4. To show `starRingEnd ℂ λ₁ * λ₂ ≠ 1`: assume equality, derive `starRingEnd ℂ λ₁ * λ₁ = 1` from `‖λ₁‖ = 1` (via `Complex.normSq_eq_norm_sq` + `Complex.inv_def` + `inv_mul_cancel₀`), multiply, conclude `λ₁ = λ₂` contradiction.
5. Close with `(mul_eq_zero.mp hfact).resolve_left (sub_ne_zero.mpr hcoeff)`.

Mathlib API used: `inner_smul_left`, `inner_smul_right`, `linear_combination`, `Complex.inv_def`, `Complex.normSq_eq_norm_sq`, `inv_mul_cancel₀`, `sub_ne_zero`, `mul_eq_zero`. The `Complex.inv_def` rewrite chain produces `starRingEnd ℂ λ = λ⁻¹` when `‖λ‖ = 1` — package this as a private helper if you need it more than once.

### Recipe: HasSum-lifted operator power via per-eigenspace induction

Use this for `U^N = Id` from `λ_i^N = 1` on a HasSum-style EigenDecomposition (where each `x ∈ H` decomposes as `HasSum (fun i => proj i x) x`).

Pattern:

1. Private helper `pow_apply_proj : (U^n) (proj i x) = (eigenval i)^n • proj i x`, by induction on `n`. Inductive step: `pow_succ'` to peel one `U`, then `ContinuousLinearMap.mul_apply`, `ContinuousLinearMap.map_smul`, `d.eigen_eq`, `smul_smul`, `← pow_succ`.
2. Apply `HasSum.mapL (U^N)` to `d.span_complete x : HasSum (proj · x) x`, getting `HasSum (fun i ↦ (U^N)(proj i x)) ((U^N) x)`.
3. Each summand collapses via `pow_apply_proj` + hypothesis `(eigenval i)^N = 1` + `one_smul`.
4. Two HasSums of the same family in a T2 space close via `HasSum.unique`.

Mathlib API used: `pow_succ'`, `ContinuousLinearMap.mul_apply`, `ContinuousLinearMap.map_smul`, `smul_smul`, `one_smul`, `HasSum.mapL` (`= ContinuousLinearMap.hasSum`), `HasSum.unique`, `funext`. Lives at `power_eq_id_of_roots_of_unity` in `Layer3a_SpectralUnitary.lean`.

### Recipe: rotational-symmetry preservation for level-set lemniscates

Use this for `ω · z ∈ {z : |P(z)| = 1}` when `ω` is an N-th root of unity and `P` is invariant under multiplication by N-th roots.

Pattern:

1. `simp only [ErdosLemniscate, Set.mem_setOf_eq, PN, Polynomial.eval_sub, Polynomial.eval_pow, Polynomial.eval_X, Polynomial.eval_C] at *` to unfold the level-set definition and the polynomial evaluator into bare `‖z^N - 1‖ = 1`.
2. Prove `ω^N = 1` via `← Complex.exp_nat_mul`, then a `show` that pins the argument to `2 * ↑Real.pi * I` (use `field_simp` alone — NOT `field_simp; ring`, which over-runs and leaves "No goals to be solved").
3. Close `ω^N = 1` with `Complex.exp_two_pi_mul_I`.
4. Apply `mul_pow`, `hωN`, `one_mul` to reduce `(ω · z)^N - 1 = z^N - 1`, close with hypothesis.

Mathlib API used: `Complex.exp_nat_mul`, `Complex.exp_two_pi_mul_I`, `mul_pow`. The `field_simp; ring` redundancy bug bit Wave 1 — `field_simp` already discharges polynomial identities of this shape; do not chain `ring` after it.

### Recipe: scope-honest doc-comment when statement is weaker than prose

Use this when a previously-axiomatized statement turns out to be the trivial corollary of a stronger axiom in the same file (e.g. axiom A implies the weak version directly), but the doc-comment was selling the strong version.

Pattern:

1. Verify the statement is weaker than the prose (e.g. `≤ 1` vs `< 1`).
2. Convert axiom to theorem with the trivial proof (`exact le_of_eq (strong_axiom args)` or similar).
3. Rewrite the doc-comment to be explicit about (a) what was actually proved (the weak version), (b) why the strong version cannot be discharged in-tree (Mathlib gap, axiom dependency chain), (c) what would need to land in Mathlib for the strong version, with full citations (Ransford 1995, Saff–Totik 1997).
4. Do NOT silently leave the prose claiming the strong version. The integrity gain is the doc-comment, not the one-line proof.

Lives at `capacity_decreases_under_perturbation` in `Layer5_LemniscateGeometry.lean`. The audit-honesty pattern: if the statement and the prose disagree, fix the prose, not the statement.

### Pitfall: `Int.cast_natAbs` collapse breaks `((q.natAbs : ℕ) : ℂ) = ±(q : ℂ)`

Wave 3 attempted to discharge `periodic_of_finite_roots_of_unity` (rational-eigenvalue spectral periodicity) and ran aground on this one cast bridge. In Mathlib v4.27.0, `mod_cast` and `push_cast` both rewrite `((q.natAbs : ℕ) : ℂ)` to `(|q| : ℂ)` via `Int.cast_natAbs`, leaving residual goals of the form `↑q.natAbs = ↑|q|` that aren't definitionally equal in a `rfl`-able way. Multiple attempts (`exact_mod_cast h1`, `conv_rhs; rw; push_cast; rfl`, `push_cast; rfl` against an `omega`-derived ℤ-fact) all hit the same wall.

**Workaround pattern (untested, for next wave):** reformulate the rationality hypothesis using `q : ℚ` directly, and use `q.den` (always positive, no `natAbs` casework needed) instead of `(q : ℤ).natAbs`. This sidesteps the simp-set fight entirely. The reduction to `power_eq_id_of_roots_of_unity` (which IS a clean theorem) is sound; only the helper bridge needs reformulating.

**Lesson:** when Mathlib's simp set has a normalization rule that collapses your goal in an unwanted direction, the answer is usually to pick a different representative of the data, not to fight the simp set. Reformulate, don't out-tactic.

### Recipe: foreground build-verification after parallel salvo

Use this after dispatching parallel subagents to discharge axioms in independent files.

Pattern:

1. Run `lake build` in foreground (NOT background piped through `tee | tail`, which masks the lake exit code with the pipeline's tee/tail exit code = 0 even when lake fails).
2. If build fails, read the full log (NOT just `tail -25` — errors might be earlier than the tail window).
3. Surgical fixes for cast/tactic-redundancy errors are usually < 5 lines per error site.
4. If a fix re-fails twice, revert to axiom with a refined doc-comment. Do not block the entire build for a single non-discharged axiom.

The pipe-tail-masking-exit-code bug is real: `lake build 2>&1 | tee /tmp/log | tail -25` returns the exit code of `tail`, which is 0 even when lake errored. Always read the actual log content for build status, not the wrapping pipeline's exit code.

---

## Recipe: finite packet certificate via Boolean enumeration (Erdos #30 Singer57)

Use this when an experiment packet has exported a small exact finite face and the goal is a compiled Lean certificate of the finite facts, not a theorem-level asymptotic result.

Pattern:

1. Pin the Lean lane before formalizing. For Google Formal Conjectures compatibility, use `leanprover/lean4:v4.27.0` and mathlib input rev `v4.27.0`.
2. Verify packet sidecars first (`shasum -a 256 -c *_RESULTS.sha256`) and copy only the finite witness literals into Lean.
3. Define small executable predicates over `List Nat`: membership, counting, no-duplicates, positive differences, ordered modular differences, set equality, and bounded affine-image search.
4. State each packet claim as `predicate = true` and close with `native_decide`; do not introduce `axiom`, `sorry`, or theorem language beyond the finite certificate.
5. Compare unordered difference skeletons by finite set equality, not raw list order. Raw positive-difference lists are order-sensitive and can falsely reject translated or reordered skeletons.
6. For modular coverage, enumerate all ordered pairs against the full source list; tail-only enumeration misses reverse differences and undercounts the `(57,8,1)` PDS check.

Landed at `Erdos30_Singer57_Certificate.lean`: C1-C6 finite checks for the SHA-verified Singer57 V2 and n58 branch packets. Build target: `lake build Erdos30_Singer57_Certificate` under Lean 4.27.0.

---

## 2026-05-01 — Face/field abstract finite lemma

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField.lean`

Build target: `lake build Erdos30_FaceField` under Lean 4.27.0 / Mathlib v4.27.0.

### Recipe: exposed-face split by two finite observables

Use this when an exact finite packet shows that two observables select different
zero-temperature witnesses on the same extremal family, and the first formal
goal is only to justify the "face, not unique witness" language.

Pattern:

1. Keep the lemma abstract over a finite family `F : Finset α`; do not import
   Sidon, B_2[g], or experimental vocabulary into the proof.
2. Define minimization as a transparent predicate:
   `IsFieldMinOn F φ x := x ∈ F ∧ ∀ y ∈ F, φ x ≤ φ y`.
3. Define a split as two chosen minimizers `x`, `y` for two ordered observables
   with `x ≠ y`.
4. Prove `{x, y} ⊆ F` by unfolding pair membership with
   `simp only [Finset.mem_insert, Finset.mem_singleton]`.
5. Convert the unequal pair to cardinality two with `Finset.card_pair`.
6. Finish the face-size certificate with `Finset.card_le_card`.

The useful theorem names are:

- `Erdos.Collider.fieldSplit_card_two_le`
- `Erdos.Collider.fieldSplit_not_card_le_one`

### Pitfall: `Mathlib.Data.Finset.Basic` is too small for card lemmas

In the current Mathlib v4.27.0 layout, importing only
`Mathlib.Data.Finset.Basic` left `.card`, `Finset.card_pair`, and
`Finset.card_le_card` unavailable. Import `Mathlib` or a sufficiently broad
Finset card module before using this recipe.

---

## 2026-05-01 — Erdos #30 n=30 face/field certificate

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_N30_Certificate.lean`

Build target: `lake build Erdos30_FaceField_N30_Certificate` under Lean 4.27.0
/ Mathlib v4.27.0.

### Recipe: packet-row certificate from exposed witnesses

Use this when an exact packet exports a small number of field-selected witnesses
and the goal is to connect them to the abstract exposed-face lemma without
claiming a theorem about the full extremal problem.

Pattern:

1. Copy the packet witnesses as literal `List Nat` values and convert them to
   `Finset Nat` with `.toFinset`.
2. Certify each witness separately: cardinality, range inclusion, and the local
   predicate (`Erdos.Sidon.IsSidonSet` here) using `native_decide`.
3. Define a small `exposedWitnessFamily : Finset (Finset Nat)` containing only
   the packet-exported witnesses being certified. Name this honestly; it is not
   the full exact face unless all face states are included.
4. Define observables as finite functions over that witness family, matching the
   packet's zero-temperature selections.
5. Do not ask `native_decide` to solve `FieldSplit` directly: the existential is
   over the infinite type `Finset Nat`, so Lean cannot synthesize the needed
   global decision procedure.
6. Instead, explicitly provide the two witnesses, prove membership by finite
   `simp` over the witness family, and prove the minimum inequalities by cases
   on membership.
7. Apply `Erdos.Collider.fieldSplit_card_two_le` and
   `Erdos.Collider.fieldSplit_not_card_le_one` to get the face-size consequence.

The useful theorem names are:

- `Erdos30FaceFieldN30Certificate.prefixWitness_sidon`
- `Erdos30FaceFieldN30Certificate.massWitness_sidon`
- `Erdos30FaceFieldN30Certificate.jointWitness_sidon`
- `Erdos30FaceFieldN30Certificate.n30_prefix_mass_field_split`
- `Erdos30FaceFieldN30Certificate.n30_exposed_family_has_at_least_two_witnesses`
- `Erdos30FaceFieldN30Certificate.n30_exposed_family_not_singleton`

### Pitfall: stale `.olean` headers after toolchain/cache drift

The first certificate build hit an incompatible-header error for
`Erdos30_Sidon_Defs.olean`. Rebuilding the dependency with
`lake build Erdos30_Sidon_Defs` fixed the cache without changing source.

---

## 2026-05-01 — Erdos #30 Sidon n=20..30 face/field window certificate

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_Window_Certificate.lean`

Build target: `lake build Erdos30_FaceField_Window_Certificate` under Lean
4.27.0 / Mathlib v4.27.0.

### Recipe: generated window certificate from exact packet rows

Use this when a packet exports one or two field-selected witnesses for every
row in a small finite window and the proof goal is an audited finite
certificate, not a theorem about the whole extremal face.

Pattern:

1. Generate the Lean literals from the exact JSON packet so every row has
   explicit `PrefixWitness` and `MassWitness` definitions.
2. Certify each row independently: witness cardinality, range inclusion, Sidon
   predicate, and prefix/mass distinction. `native_decide` is appropriate for
   these finite Boolean checks.
3. Reuse a generic two-point family helper:
   `pairFamily x y := ({x, y} : Finset (Finset Nat))`.
4. Define two small observables that select the left and right endpoint of that
   certified pair. This keeps the field-response proof finite and avoids
   pretending the pair is the full exact face.
5. Prove each row's `FieldSplit` by explicit witness pair, then apply
   `Erdos.Collider.fieldSplit_card_two_le` for the two-witness consequence.
6. Aggregate the window theorem only after every row theorem is available.

The useful theorem names are:

- `Erdos30FaceFieldWindowCertificate.sidon_window_20_30_all_prefix_mass_split`
- `Erdos30FaceFieldWindowCertificate.sidon_window_20_30_all_exposed_pairs_have_two_witnesses`

### Pitfall: right-associated conjunctions

Lean parses chained conjunction statements as right-associated. A generated
proof for an 11-row window must therefore emit nested pairs in the same shape:
`⟨p20, ⟨p21, ... ⟨p29, p30⟩...⟩⟩`. A left-nested generated proof fails even
though all row lemmas exist.

### Pitfall: existential field splits over infinite carrier types

As in the n=30 certificate, do not ask `native_decide` to solve
`Erdos.Collider.FieldSplit` directly when the carrier is `Finset Nat`. Give the
two finite witnesses explicitly, prove pair-family membership by `simp`, and
let the abstract lemma carry only the finite face-size consequence.

---

## 2026-05-01 — Erdos #30 n=57/58 branch face/field certificate

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_57_58_Certificate.lean`

Build target: `lake build Erdos30_FaceField_57_58_Certificate` under Lean
4.27.0 / Mathlib v4.27.0.

### Recipe: branch certificate over exported witness skeletons

Use this when the exact packet has already exported a small branch/family table
and an earlier certificate has named the witness lists.

Pattern:

1. Import the finite witness certificate instead of copying all witness lists
   again. For this lane, reuse `Erdos30_Singer57_Certificate`.
2. Define only the branch witnesses that participate in the new field-response
   statement, such as `N57_3`, `N57_5`, `N58_2`, `N58_7`, and `N58_9`.
3. Certify local facts with `native_decide`: Sidon/range/card checks, distinct
   field witnesses, same/different positive-difference skeletons, and
   translation-chain identities.
4. Keep observable minimization over a two-point exported branch family unless
   the whole face has been imported into Lean. Do not call this a full face
   theorem.
5. State the packet-index checks explicitly as finite Boolean claims so a later
   reviewer can trace the Lean certificate back to the handoff tables.
6. Package the final theorem as a Boolean finite certificate, not as an
   asymptotic or extremal theorem.

The useful theorem names are:

- `Erdos30FaceField5758Certificate.n57_mass_joint_field_split`
- `Erdos30FaceField5758Certificate.n57_mass_joint_share_positive_difference_skeleton`
- `Erdos30FaceField5758Certificate.n58_prefix_mass_field_split`
- `Erdos30FaceField5758Certificate.n58_prefix_branch_has_new_positive_difference_skeleton`
- `Erdos30FaceField5758Certificate.n58_pareto_chain_persists_from_n57`
- `Erdos30FaceField5758Certificate.finite_57_58_branch_field_response_certificate`

### Pitfall: branch data is not automatically a full-face theorem

The `n=57/58` rows are valuable because the exported witness set is small and
structurally rich, but the Lean module still certifies the named exported
branch relationships only. The next stronger climb requires importing every
face member and every observable score needed to prove the minimizers directly
from the complete finite table.

---

## 2026-05-01 — Erdos #30 n=56..58 full exported-face certificate

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_56_58_FullFace_Certificate.lean`

Build target: `lake build Erdos30_FaceField_56_58_FullFace_Certificate` under
Lean 4.27.0 / Mathlib v4.27.0.

### Recipe: lift a complete exported face into Lean

Use this after a finite exact packet exports every member of a small ground
face and records the selected winners for several observables.

Pattern:

1. Generate one `Finset Nat` definition for every exported witness. Keep the
   witness names index-aligned with the source packet.
2. Define the exported face as the `List.toFinset` of those witnesses and prove
   its cardinality. This catches duplicate witness literals.
3. Separately check every witness has the packet cardinality, lies in `[0,n]`,
   and satisfies `Erdos.Sidon.IsSidonSet`.
4. Work over a finite index face, `Finset.range k`, for observable-selection
   theorems. This avoids carrying large `Finset Nat` equality terms through
   every minimizer proof while the indexed witness map preserves the link back
   to the exported face.
5. Convert packet score columns into rank-coded observables. These ranks
   preserve the packet ordering for prefix, density-adjusted mass, and joint
   scores on the exported face. They are not raw floating-point theorem
   statements.
6. Recompute prefix/mass/joint winners and Pareto minima from the ranks, then
   compare those lists to the packet winner lists with `native_decide`.
7. Prove `IsFieldMinOn` by explicit bounded index cases:
   `simp [n57FaceIndices] at hy; interval_cases y <;> native_decide`.
8. Use `Erdos.Collider.FieldSplit` on the full indexed face, not just a
   two-point pair, once the minimizer lemmas are available.

The useful theorem names are:

- `Erdos30FaceField5658FullFaceCertificate.n56_exported_face_card`
- `Erdos30FaceField5658FullFaceCertificate.n57_exported_face_card`
- `Erdos30FaceField5658FullFaceCertificate.n58_exported_face_card`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_winner_table_matches_packet`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_split_pattern_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_mass_winner_table_matches_packet`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_mass_split_pattern_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_winner_table_matches_packet`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_strict_nonwinner_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_zero_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_packet_probe_separation_is_actual_residual_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_full_prefix_zero_and_probe_separation_certificate`
- `Erdos30FaceField5658FullFaceCertificate.face_handoff_56_58_exact_prefix_probe_exact_mass_split_certificate`
- `Erdos30FaceField5658FullFaceCertificate.n56W3_prefix_count_table_matches_witness`
- `Erdos30FaceField5658FullFaceCertificate.n56W3_full_prefix_segment_bound_certificate`
- `Erdos30FaceField5658FullFaceCertificate.n56W3_full_prefix_residual_zero`

### Recipe extension: replace one rank-coded observable with an exact formula

Compiled helper source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_ExactObservables.lean`

Build target: `lake build Erdos30_FaceField_ExactObservables` under Lean
4.27.0 / Mathlib v4.27.0.

The first exact-observable lift is density-adjusted mass. Use the integer key
`|2 * sum(A) - n * (|A| + 1)|`, implemented as
`Erdos30FaceFieldExactObservables.densityAdjustedMassTwice`. This is exactly
twice the packet's mass deviation, so it preserves the minimizer order without
importing floats or rationals.

Pattern:

1. Keep the witness import and finite index face from the full-face recipe.
2. Define `nXXExactMassTwice i` by applying the exact helper to
   `nXXWitnessOfIndex i`; do not hardcode the theorem result as the definition.
3. Generate the expected exact integer table only in the theorem statement.
4. Prove the packet mass winner list and the rank/mass agreement with
   `native_decide`.
5. Rebuild field-split theorems with `ExactMassTwice` on one axis. This proves
   that the split survives when the mass field is the mathematical formula, not
   a packet score rank.

### Recipe extension: exact zero prefix probes

The first exact-prefix lift is a probe certificate, not yet the full residual
theorem. Define
`Erdos30FaceFieldExactObservables.prefixResidualProbeCard10 n t p` for the
packet-selected prefix location `t` and prefix count `p`.

For the `n=56..58` handoff face, every packet prefix winner has terminal probe
`t=n`, `p=10`. Prove these probes are exactly zero using:

1. A square-root lower bound such as `(n/10 : ℝ) ≤ Real.sqrt n`, discharged by
   `Real.le_sqrt` and `norm_num`.
2. The nonpositive terminal difference
   `(n : ℝ) - 10 * Real.sqrt n ≤ 0`, discharged by `nlinarith`.
3. `abs_of_nonpos` plus `linarith` to close the zero residual.

Then define `nXXExactPrefixProbe` and rebuild the split theorem with
`ExactPrefixProbe` on one axis and `ExactMassTwice` on the other. This proves
the split survives when both visible axes are exact finite formulas, while
leaving the full global prefix-residual theorem as the next climb.

### Recipe extension: full-prefix segment certificates for n=56..58

The full-prefix lift now covers every packet-selected exact-prefix winner in
the exported handoff face: `n56W3`, `n57W4`, `n57W5`, `n58W2`, `n58W8`, and
`n58W9`. The direct proof that unfolds `prefixCount` at every cutoff is too
expensive; Lean times out while repeatedly normalizing finite filters. Use the
certificate shape instead:

1. Generate the complete prefix-count table
   `(List.range 57).map (prefixCount n56W3)`.
2. Generate the same table as a compact model function
   `n56W3PrefixCountTable`.
3. Split the table into constant-prefix segments, for example `33..45` with
   prefix count `6`.
4. Prove each segment bound
   `|t - p * sqrt(56)| <= 10 * sqrt(56) - 56`
   using only the integer interval bounds on `t` and a rational lower bound
   for the square root (`7 <= sqrt(56) < 8` for n=56, and
   `15 / 2 <= sqrt(n) < 8` for n=57/n=58).
5. Package the segment lemmas as
   `nXXWY_full_prefix_segment_bound_certificate`.

This is the compiled bridge from packet-selected prefix probe to full prefix
profile bounds for the 56..58 handoff face. The current certificate includes
the table-match theorems and segment-bound theorems for all six selected
prefix winners.

### Recipe extension: turn segment bounds into full-prefix residual-zero facts

Once the segment bounds compile, the next move is not another scan. It is a
local algebra bridge:

1. Prove the selected witness has the exact terminal drift
   `prefixDrift n A = 10 * sqrt(n) - n`. The proof uses the exact card-10
   fact, `sqrt(n) < 8`, `abs_of_nonneg`, `max_eq_left`, and
   `(sqrt(n))^2 = n`.
2. For each constant prefix-count segment, prove the exact `prefixCount`
   equality by `interval_cases t <;> native_decide`.
3. Unfold `prefixResidualAt`, rewrite the drift and prefix count, normalize
   casts with `norm_num at hdev ⊢`, and close with the already compiled
   segment inequality.
4. Quantify over all cutoffs with `interval_cases t`, dispatching each literal
   cutoff to its segment-zero theorem.
5. Link each packet probe back to the real residual at the packet cutoff by
   unfolding `prefixResidualProbeCard10` and `prefixResidualAt`, rewriting
   drift, and proving the local prefix count with `native_decide`.

This compiled the full-prefix residual-zero certificates for `n56W3`,
`n57W4`, `n57W5`, `n58W2`, `n58W8`, and `n58W9`, plus actual positive residual
at packet cutoffs for every exact-prefix non-winner.

### Recipe extension: package full-prefix residuals as a scalar finite field

After pointwise residual-zero and positive-cutoff facts compile, package the
prefix field as a real finite maximum instead of leaving it as a predicate:

```lean
noncomputable def fullPrefixResidualMax (n : Nat) (A : Finset Nat) : ℝ :=
  (Finset.range (n + 1)).sup' (by exact ⟨0, by simp⟩)
    (fun t => prefixResidualAt n A t)
```

The useful helper facts are:

1. `fullPrefixResidualMax_eq_zero_of_fullPrefixResidualZero`: use
   `Finset.sup'_le` for the upper bound and `Finset.le_sup'` at `t = 0` for
   the lower bound.
2. `fullPrefixResidualMax_pos_of_pos_at`: put the positive cutoff into
   `Finset.range (n + 1)` and apply `Finset.le_sup'`.
3. `fullPrefixResidualMax_nonneg`: use nonnegativity at `t = 0` and
   `Finset.le_sup'`.
4. `fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at`: use `Finset.sup'_le`
   for a generated pointwise upper bound, then use `Finset.le_sup'` at the
   attaining cutoff to prove equality with the packet probe value.
5. Define `nXXFullPrefixResidualMax i` by applying the helper to
   `nXXWitnessOfIndex i`, then prove `IsFieldMinOn` from winner-zero plus
   global nonnegativity.

This compiled the scalar prefix-field certificates:

- `face_handoff_56_58_full_prefix_scalar_zero_certificate`
- `face_handoff_56_58_full_prefix_scalar_positive_certificate`
- `face_handoff_56_58_full_prefix_scalar_min_certificate`
- `face_handoff_56_58_full_prefix_scalar_exact_mass_split_certificate`

### Recipe extension: lift packet-probe joint to scalar full-prefix joint

Once the scalar maximum exists, the joint lift becomes a finite-face equality
problem rather than a new search.

Pattern:

1. For each non-prefix witness, generate segment upper bounds proving every
   cutoff residual is at most the packet exact-prefix probe value.
2. Use `fullPrefixResidualMax_eq_of_pointwise_le_and_eq_at` with the packet
   cutoff equality to prove the scalar full-prefix residual maximum equals the
   exact packet prefix probe.
3. Define `nXXScalarFullPrefixJointKey i` as
   `nXXFullPrefixResidualMax i * Real.sqrt n + nXXExactMassTwice i / 2`.
4. Prove `nXXScalarFullPrefixJointKey i = nXXExactJointKey i` for every
   exported witness. This transfers the old exact-probe joint comparisons to
   the scalar full-prefix observable.
5. Reuse the strict exact-joint comparisons to prove scalar joint winners.
6. If `nlinarith` stalls in the segment upper bounds, add tighter rational
   square-root brackets. In this row, `sqrt(56) < 15 / 2` and
   `sqrt(57), sqrt(58) < 23 / 3` were enough.

This compiled the scalar full-prefix joint certificates:

- `face_handoff_56_58_full_prefix_scalar_matches_packet_probe_certificate`
- `face_handoff_56_58_scalar_joint_matches_probe_joint_certificate`
- `face_handoff_56_58_scalar_full_prefix_joint_winner_certificate`

### Recipe extension: abstract scalar-joint minimizer transfer

The scalar-joint lift should not duplicate strict comparisons after the exact
probe joint minimizer has already compiled. Add the reusable congruence lemma
near `IsFieldMinOn`:

```lean
theorem isFieldMinOn_of_eq_on [LinearOrder β]
    {F : Finset α} {φ ψ : α → β} {x : α}
    (hmin : IsFieldMinOn F φ x)
    (heq : ∀ y ∈ F, ψ y = φ y) :
    IsFieldMinOn F ψ x
```

Then use the weighted-joint wrapper:

```lean
theorem isFieldMinOn_weightedJoint_of_prefix_eq_on
    {ι : Type*} {F : Finset ι} {x : ι}
    {prefixScalar prefixProbe mass : ι → ℝ} {weight : ℝ}
    (hmin :
      IsFieldMinOn F
        (fun i => prefixProbe i * weight + mass i) x)
    (hprefix : ∀ i ∈ F, prefixScalar i = prefixProbe i) :
    IsFieldMinOn F
      (fun i => prefixScalar i * weight + mass i) x
```

For each generated row, prove one face-local equality theorem:

```lean
theorem n58_full_prefix_residual_max_matches_exact_prefix_probe_on_face :
    ∀ i ∈ n58FaceIndices, n58FullPrefixResidualMax i = n58ExactPrefixProbe i
```

The scalar joint minimizer theorem then becomes one transfer from the exact
joint minimizer, instead of a second copy of all strict comparisons. This was
then used in the compiled `n=59` microcertificate:

- `Erdos30FaceField59FullFaceCertificate.n59_micro_scalar_full_prefix_joint_winner_certificate`
- `Erdos30FaceField59FullFaceCertificate.n59_micro_scalar_joint_matches_probe_joint_certificate`
- `Erdos30FaceField59FullFaceCertificate.n59_full_prefix_residual_max_matches_exact_prefix_probe_on_face`

The `n=59` row also exposed the first bracket-tightening lesson beyond the
original handoff: the old `sqrt(n) >= 15 / 2` lower bound is too weak for
some `p = 4` prefix segments; `sqrt(59) >= 23 / 3` and
`sqrt(59) < 39 / 5` close the generated segment bounds.

The same generator pattern then compiled for `n=60`:

- `Erdos30FaceField60FullFaceCertificate.n60_micro_scalar_full_prefix_joint_winner_certificate`
- `Erdos30FaceField60FullFaceCertificate.n60_micro_scalar_joint_matches_probe_joint_certificate`
- `Erdos30FaceField60FullFaceCertificate.n60_full_prefix_residual_max_matches_exact_prefix_probe_on_face`

The `n=60` row is the first useful warning that prefix selection is now a
surface phenomenon, not a point phenomenon: the packet has 21 prefix winners,
while exact mass and scalar/probe joint both select W43. The point-valued
minimizer transfer still works because the joint winner is unique at this row,
but the next rows must be stated as selected minimizer surfaces.

### Pitfall: rank-coded scores are an audit bridge, not raw analysis

The Lean certificate now proves the mass winners from the exact integer formula,
proves zero exact prefix probes for the prefix winners, proves strict positive
exact-prefix probe gaps for every non-winner, proves exact joint-key winners
from the exact prefix probe plus exact integer mass, proves that every
packet-selected prefix winner in n=56..58 has full prefix profile bounded by
the terminal drift, proves full-prefix residual zero for those winners, and
proves the probe non-winners have positive actual residual at their packet
cutoffs. The scalar prefix field now exists and compiles, and the scalar
full-prefix joint lift now compiles for the complete exported n=56..58 face
and the complete exported `n=59` and `n=60` faces. The next honest data climb
is not another unique-winner import: `n=61+` needs tie-aware minimizer language
or a narrower certificate surface.

### Recipe extension: switch from selected points to minimizer surfaces

Use this when a field has a tied selected face. Do not force a unique witness
if the packet reports multiple exact winners.

Minimal abstract shape:

```lean
def IsFieldMinimizerSet [LinearOrder β]
    (F S : Finset α) (φ : α → β) : Prop :=
  S ⊆ F ∧ ∀ x : α, x ∈ S ↔ IsFieldMinOn F φ x
```

The generator should then prove all selected indices minimize, and every
off-surface index is strictly beaten by at least one selected index. For
`n=61+`, keep the certificate index/scalar-based first; do not emit thousands
of full witness literals just to certify the selected surface over packet
observables.

This now compiles for `n=61` as the first tie-aware surface certificate:

- `Erdos30FaceField61JointSurfaceCertificate.n61_joint_surface_minimizer_certificate`
- `Erdos30FaceField61JointSurfaceCertificate.n61_joint_selected_surface_certificate`

The important mathematical shift is that the certified object is a selected
surface, not a selected point. For `n=61`, the packet-backed joint/Pareto
surface is `[85, 96, 113, 139]` over `Finset.range 152`; mass has one extra
mass-only point `79`. The Lean certificate proves the rank-table minimizer
surface over exported indices only. It deliberately does not import witness
literals or reprove the enumerator.

### Recipe extension: scale selected-surface certificates by index ranks

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_FaceField_61_64_JointSurface_Certificate.lean`

Build target: `lake build Erdos30_FaceField_61_64_JointSurface_Certificate`
under Lean 4.27.0 / Mathlib v4.27.0.

The same selected-surface pattern now compiles for `n=61..64` over exported
faces of size `152,398,1022,2360`. The certified joint surfaces have sizes
`4,6,10,11` respectively. This is the first useful scale step past the
single-row `n=61` certificate: it says the observable-selected object remains
a tiny proper surface of the exact face, not a single optimizer and not the
whole face.

Useful theorem names:

- `Erdos30FaceField6164JointSurfaceCertificate.n61_joint_surface_minimizer_certificate`
- `Erdos30FaceField6164JointSurfaceCertificate.n62_joint_surface_minimizer_certificate`
- `Erdos30FaceField6164JointSurfaceCertificate.n63_joint_surface_minimizer_certificate`
- `Erdos30FaceField6164JointSurfaceCertificate.n64_joint_surface_minimizer_certificate`
- `Erdos30FaceField6164JointSurfaceCertificate.n64_joint_selected_surface_certificate`

Pitfall: do not ask Lean to unfold a 2,360-branch rank function inside `simpa`.
The successful pattern is to prove the anchor selected index has rank zero by
`native_decide`, rewrite the anchor rank to `0`, and then apply the off-surface
positive-rank certificate. This keeps the certificate compact enough to compile
while still proving the minimizer-surface predicate.

### Recipe extension: certify selected-surface transitions by partition

Compiled sources:

- `erdos-experiments/Erdos30/lean/Erdos30_FaceField_61_64_JointSurface_Transition_Certificate.lean`
- `erdos-experiments/Erdos30/lean/Erdos30_FaceField_65_71_JointSurface_Transition_Certificate.lean`

Build targets:

- `lake build Erdos30_FaceField_61_64_JointSurface_Transition_Certificate`
- `lake build Erdos30_FaceField_65_71_JointSurface_Transition_Certificate`

Both compile under Lean 4.27.0 / Mathlib v4.27.0.

Use this after the selected-surface certificate compiles and the next question
is not "which indices minimize?" but "how did the selected surface change from
the previous row?" Keep the expensive witness lookup in the generator: classify
selected indices as direct previous-row, `+1` previous-row, or new, then emit
only compact index partitions and count tables to Lean.

Useful theorem names:

- `Erdos30FaceField6164JointSurfaceTransitionCertificate.n61_joint_surface_transition_by_indices`
- `Erdos30FaceField6164JointSurfaceTransitionCertificate.n64_joint_surface_transition_by_indices`
- `Erdos30FaceField6164JointSurfaceTransitionCertificate.all_surface_transition_counts_balance`
- `Erdos30FaceField6164JointSurfaceTransitionCertificate.all_surface_transitions_are_mixed`
- `Erdos30FaceField6164JointSurfaceTransitionCertificate.all_direct_previous_selected_surfaces_empty`
- `Erdos30FaceField6164JointSurfaceTransitionCertificate.joint_surface_transition_certificate_passes`
- `Erdos30FaceField6571JointSurfaceTransitionCertificate.n71_joint_surface_transition_by_indices`
- `Erdos30FaceField6571JointSurfaceTransitionCertificate.plus_one_previous_surface_count_table`
- `Erdos30FaceField6571JointSurfaceTransitionCertificate.new_surface_count_table`
- `Erdos30FaceField6571JointSurfaceTransitionCertificate.joint_surface_transition_certificate_passes`

The `n=65..71` extension keeps the same finite pattern over much larger exact
faces: direct previous-row selected-surface persistence remains empty,
`+1`-inherited counts are `32,24,85,74,179,125,324`, new selected counts are
`15,16,35,35,88,49,128`, and the joint-selected surface stays below one
twentieth of the exported face through `n=71`.

The key pattern is to make Lean prove a finite partition and its count
consequences, while the generator remains responsible for the artifact-backed
witness-label lookup. That preserves the honest boundary: this certifies
selected-surface transition bookkeeping, not the enumerator and not an
asymptotic Sidon theorem.

---

## 2026-05-01 — Erdos #30 n=58..71 ground-face branch summary certificate

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_GroundFaceBranch_58_71_Certificate.lean`

Build target: `lake build Erdos30_GroundFaceBranch_58_71_Certificate` under
Lean 4.27.0 / Mathlib v4.27.0.

### Recipe: certify the branch table before importing a huge face

Use this when the exact packet exports complete finite faces but the full
witness lift is too large for a useful first Lean certificate. The 58..71
packet reaches 203,840 exact faces at n=71, so the right foothold is a compact
table certificate, not hundreds of thousands of witness literals.

Pattern:

1. Generate finite summary functions for face counts, exported counts, skeleton
   counts, persistence counts, and observable-winner counts.
2. Prove exported counts match exact face counts for every row. This certifies
   that the packet is complete at the table level.
3. Prove branch persistence by adjacent edge: every post-58 row contains all
   previous faces and all previous faces shifted by `+1`, by count equality
   against the previous exact face count.
4. Prove skeleton counts strictly increase across the branch. This is the
   current formal handle on the branch becoming structurally richer.
5. Keep the selected surfaces as count-level facts: Pareto and joint-selected
   sets are proper subsets of the full face, and the joint count is no larger
   than the Pareto count. Do not assert equality; n=58 and n=59 have smaller
   joint counts than Pareto counts.
6. Close the certificate as a Boolean finite table theorem with `native_decide`.

The useful theorem names are:

- `Erdos30GroundFaceBranch5871Certificate.branch_table_face_counts`
- `Erdos30GroundFaceBranch5871Certificate.all_exports_are_complete_by_count`
- `Erdos30GroundFaceBranch5871Certificate.every_post58_face_contains_previous_face`
- `Erdos30GroundFaceBranch5871Certificate.every_post58_face_contains_previous_face_plus_one`
- `Erdos30GroundFaceBranch5871Certificate.skeleton_count_strictly_increases_58_71`
- `Erdos30GroundFaceBranch5871Certificate.pareto_surface_is_proper_subset_by_count_58_71`
- `Erdos30GroundFaceBranch5871Certificate.joint_count_no_larger_than_pareto_count_58_71`
- `Erdos30GroundFaceBranch5871Certificate.branch_summary_certificate_passes`

### Pitfall: a compact certificate should not overstate set equality

The branch packet supports complete exported counts, persistence counts, and
winner-count summaries. It does not by itself prove literal set containment of
joint winners inside Pareto winners in Lean. Until the actual selected witness
sets are imported, state the supported count fact: joint count no larger than
Pareto count.

### Pitfall: do not import massive faces just to prove the next climb exists

The full-face method worked for n=56..58 because the exported faces were small.
For n=58..71, the table grows from 10 to 203,840 faces. The table certificate
keeps the proof lane auditable while preserving the harder next target: exact
symbolic observables or a scalable compressed witness-set representation.

---

## 2026-05-02 - B+ attack-vector compile audit extraction

Builds rerun for the B+ audit packet:

- `lake build` in `Lean4/transdimensional-painter`: PASS under Lean 4.27.0 / Mathlib v4.27.0.
- `lake build` in `export-packets/erdos30`: PASS under Lean 4.24.0 / Mathlib commit `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`.
- `lake build` in `erdos-experiments/private-attacks/lean4`: PASS under Lean 4.24.0 / same Mathlib commit.
- `lake build` in `erdos-experiments/Erdos30`: PASS under Lean 4.27.0 / Mathlib v4.27.0, but the workspace includes `scratch/Erdos30_SpectralSidon.lean` with 4 `sorry`s. Public COMPILED status attaches only to zero-sorry targets.

Recipe extraction:

1. Treat finite packet certificates and theorem-level formalizations as different recipe families. The Erdos30 face-field certificates are excellent `native_decide` finite-table certificates, but they do not prove #30.
2. Keep scratch spectral files out of public proof counts until all `sorry`s are removed. A workspace build PASS is not enough when a target compiles with `sorry`.
3. For B+ positioning, the strongest formalization claim is not "solved" but "map finite packet -> Lean certificate -> proof target." This keeps the A/B/C distinction intact.
4. When a packet contains pass, control, and demotion examples, record all three together. The demotion is part of the recipe because it proves the pipeline can reject attractive morphisms.

Reusable pattern:

```
finite experiment packet
  -> SHA sidecar verification
  -> compact generated Lean predicates
  -> native_decide table certificate
  -> public claim only at finite-certificate scope
```

Pitfall: do not count API-reported or cached registry rows as public COMPILED numbers if the direct D1 verification query was unavailable in the session. Use the API observation for local audit context; rerun the direct D1 query before external copy.

## 2026-05-02 - Erdos30_Complete repair + sidon_elem_bound closure (Mathlib v4.27.0)

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_Complete.lean`, `erdos-experiments/Erdos30/lean/Erdos30_Lindstrom.lean`

Outcome: Erdos30_Complete failed to build under Mathlib v4.27.0 due to API drift (5 distinct kinds of breakage). The `sidon_elem_bound` axiom in Erdos30_Lindstrom — flagged in the prior recipe (2026-04-20) as "do not import Erdos30_Complete, use sidon_distinct_differences instead" — was closed by repairing Complete and adding a one-line bridge theorem. The earlier pitfall is now stale; importing Erdos30_Complete is safe again.

### Recipe: Sidon difference-counting via diff-injectivity → pigeonhole

Use this when the goal is `A.card * (A.card - 1) ≤ 2 * N` for Sidon A ⊆ Finset.range (N + 1).

Pattern:

1. Define the pair finset: `set pairs := (A ×ˢ A).filter (fun p : ℕ × ℕ => p.2 < p.1)`.
2. Prove diff-injectivity (`Set.InjOn`) on `↑pairs` for `diff_map p := p.1 - p.2`. The proof of difference-injectivity from sum-Sidon needs case analysis on orderings (`b₂ ≤ a₁` × `b₁ ≤ a₂`), with two of the four cases impossible (use `absurd h_sidon.right.symm (Nat.ne_of_lt hlt₁)` to dispatch).
3. Show `pairs.image diff_map ⊆ Finset.Icc 1 N` by destructuring membership; **must** insert `show 1 ≤ a - b ∧ a - b ≤ N` after `rw [Finset.mem_Icc]` to beta-reduce `diff_map (a, b)` so omega can see `a - b`. Without the `show`, omega sees an opaque function call and fails.
4. `Finset.card_image_of_injOn h_inj` collapses the image card to `pairs.card`.
5. `pairs.card = A.card * (A.card - 1) / 2` via bijection to `Finset.powersetCard 2 A`. The bijection function is `fun p _ => ({p.1, p.2} : Finset ℕ)`. The hard step is the injectivity branch of `Finset.card_bij` — see the next recipe for the v4.27 fix.
6. Combine: `pairs.card ≤ N` from steps 3-4; multiply by 2 and use `Even (A.card * (A.card - 1))` to drop the `/2`.

This is `sidon_difference_count` in `Erdos30_Complete.lean`. Reusable for any "Sidon-like" extremal bound where injective image-counting on differences gives the cardinality bound (#166 sum-free, #755 B_h[g], #1 distinct-subset-sums).

Pitfall: `Even` in Mathlib v4 unfolds to `∃ r, a = r + r` (not `∃ r, a = 2 * r`). The witness must satisfy the `r + r` form, not `2 * r`.

### Recipe: Mathlib v4.27.0 API drift patches for Erdos30_Complete

Five distinct breakages, each with a concrete fix. Apply this list when porting any Lean 4 file from a 2026-March pin (around `f897ebcf72`) to v4.27.0.

| Symptom | Where | Fix |
|---------|-------|-----|
| `Eq.symm h.right` has unexpected type in case-split that "is impossible" | `sidon_diff_injective` cases 1.2 and 2.1 | Replace `exact ⟨h.right.symm, h.left.symm⟩` with `exact absurd h.<side>.symm (Nat.ne_of_lt hlt)` — the original code was wrong even at proposition level; v4.27 just made the unsoundness visible. |
| `simp [Finset.mem_filter, Finset.mem_product] at hp` makes no progress and downstream `hp.1`, `hp.2.1`, `hp.2.2` fail | Anywhere mem of a filtered product is destructured | Replace `simp [...]` with `rw [Finset.mem_filter, Finset.mem_product]` (no auto-flatten). Then access nested form: `hp.1.1` (a ∈ A), `hp.1.2` (b ∈ A), `hp.2` (b < a). |
| `simp [Finset.mem_image, Finset.mem_filter, Finset.mem_product] at hd` over-simplifies, leaving `hd : (0, 0) ∈ pairs ∧ diff_map (0, 0) = d` (instantiated to (0,0)!) | `h_range` membership proof | Replace whole simp+obtain with `rw [Finset.mem_image] at hd; obtain ⟨⟨a, b⟩, hp, rfl⟩ := hd; rw [Finset.mem_filter, Finset.mem_product] at hp; obtain ⟨⟨ha, hb⟩, hlt⟩ := hp`. Step-by-step rather than one big simp. |
| `heq ▸ Finset.mem_insert_self a₁ _` fails with "rewrite did not find pattern" — heq is the unreduced `(fun p x => {p.1, p.2}) (a₁, b₁) h₁ = ...` from `Finset.card_bij` | Injectivity branch of `Finset.card_bij` | Beta-reduce heq: `change ({a₁, b₁} : Finset ℕ) = ({a₂, b₂} : Finset ℕ) at heq`. After this, `rw [← heq]` works in tactic-mode bridges. Alternative: avoid `▸` entirely and use rcases on `a₁ ∈ {a₂, b₂}` explicitly. |
| `Finset.pair_comm` typeclass instance stuck on `DecidableEq ?m` | Surjectivity branch of `Finset.card_bij`, when proving `{y, x} = {x, y}` | Provide explicit args: `exact Finset.pair_comm y x` (not bare `Finset.pair_comm`). The function takes `(a b : α)` explicitly; the typeclass `[DecidableEq α]` can't infer α without them. Also `change` to beta-reduce the goal first: `change ({y, x} : Finset ℕ) = ({x, y} : Finset ℕ)`. |
| `omega` fails to prove ring identities like `(2*m+1) * (2*m+1 - 1) = (2*m+1)*m + (2*m+1)*m` | `Even (A.card * (A.card - 1))` odd-card branch | omega doesn't handle multiplication of variables. Pull the Nat-sub out as a hypothesis, then ring: `have : A.card - 1 = 2 * m := by omega; rw [this, hm]; ring`. |
| `Nat.lt_succ_of_le (Nat.div_le_div_right this)` fails with "expected `a / t ≤ N / t`" | `Erdos30_Lindstrom.scaled_range` line 137 | After `simp [Finset.mem_range]`, the `< n + 1` form reduces to `≤ n` automatically. Drop the wrapper: `exact Nat.div_le_div_right this`. |

### Recipe: Bridging an axiom by importing a sibling-file theorem (after sibling repair)

Use this when a current file declares `axiom X` with comment "proved in Sibling.lean as Y, but Sibling has Mathlib drift", and Sibling has now been repaired.

Pattern:

1. Add `import Sibling` to the importer file.
2. Replace `axiom X (args) : T` with `theorem X (args) : T := Sibling.namespace.Y arg_remap` — the body is a one-liner.
3. Match arg names exactly to the axiom's signature (Lean's Eq doesn't care about display names but a good replacement matches for searchability).
4. Update the docstring with a closure date and reference to the proof location.
5. Remove the closed axiom from `AXIOM_INVENTORY.md` (or mark CLOSED), and from any per-file axiom-table in headers.

Concrete example: `sidon_elem_bound` in `Erdos30_Lindstrom.lean:122` was closed by:

```lean
-- before:
axiom sidon_elem_bound (A : Finset ℕ) (M : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (M + 1)) : A.card * (A.card - 1) ≤ 2 * M

-- after:
theorem sidon_elem_bound (A : Finset ℕ) (M : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (M + 1)) : A.card * (A.card - 1) ≤ 2 * M :=
  sidon_difference_count A M hS hA
```

Plus `import Erdos30_Complete` in the file header.

Pitfall: the prior recipe (2026-04-20) said "do not import Erdos30_Complete." That guidance was correct at the time because Complete was broken. After repair (2026-05-02), the import is safe and is the canonical bridge. **Recipes have shelf-lives; check the date and re-validate against the current Mathlib pin.**

## 2026-05-02 (afternoon) - order_diff_counting closure (Lindström 1969 §2)

Compiled source: `erdos-experiments/Erdos30/lean/Erdos30_Lindstrom.lean` (post-closure, lines 415-870)

Outcome: the `axiom order_diff_counting` was replaced by a fully-proven theorem. Diff +454/−19 over the file. The proof structure recurs in any "Sidon-like extremal counting under multiple-order differences" setting, so worth extracting cleanly.

### Recipe: Sigma-bijection cardinality for indexed pair sets

When you need to count a Finset defined as `{f(r,i) : 1 ≤ r ≤ ℓ, i + r < k}` and the index range varies with the outer index, do not flatten into a single Finset.range computation. Instead, model it as the underlying sigma type and bijection to it.

Pattern:
1. Define `orderPairs A ℓ : Finset (ℕ × ℕ)` as `(Finset.range A.card ×ˢ Finset.range A.card).filter (fun p => 0 < p.2 - p.1 ∧ p.2 - p.1 ≤ ℓ)` — flat, easy to count via `Finset.card_filter`.
2. Define `orderDiffs A ℓ : Finset ℕ` as `orderPairs A ℓ |>.image (fun p => oGet A p.2 - oGet A p.1)` where `oGet` is the sorted enumeration.
3. Card collapses through `Finset.card_image_of_injOn` once you have Sidon-distinct-difference injectivity.
4. The injection on differences pulls back to the sigma-type witness: `(r, i) ↦ {a_i, a_{i+r}}` with both pair members in A, and the Sidon hypothesis prevents collisions across distinct (r, i).
5. The arithmetic identity `(ℓ * (2 * A.card - ℓ - 1)) / 2 = ∑_{r=1}^ℓ (A.card - r)` falls out via `zify` + `Finset.sum_range_id_mul_two` once the count is in ℤ.

Pitfall: do NOT try to define `orderDiffs` as a finmap on a sigma type directly — Lean's elaboration in v4.27 stalls on the Σ-DepElim with `Finset.image`. The image-of-filter pattern keeps everything on `Finset (ℕ × ℕ)` and is much faster.

### Recipe: Row-decomposition telescoping for sum bounds

When you need to bound `∑_{r=1}^ℓ ∑_{i=0}^{k-r-1} (a_{i+r} - a_i)`, do not try to telescope the inner sum directly — it doesn't telescope (different shifts). Instead, partition the OUTER sum and reindex.

Pattern:
1. For fixed `r`, the inner sum equals `(a_{k-1} + a_{k-2} + ⋯ + a_{k-r}) − (a_0 + a_1 + ⋯ + a_{r-1})` by carrying the index across.
2. This is a sum of `r` paired terms `a_{k-1-j} - a_j` for `j = 0,…,r-1`. Each pair is between 0 and N (since `0 ≤ a_j ≤ a_{k-1-j} ≤ N`), so the inner sum is `≤ r·N`.
3. Outer sum over r: `∑_{r=1}^ℓ r·N = N · ℓ(ℓ+1)/2`.
4. In Lean, the cleanest implementation is `Finset.sum_range_sub` (the sub-form with the shifted index) — Mathlib v4.27 has the right shape via `(Finset.range (k - r)).sum (fun i => oGet A (i + r) - oGet A i)`.

Pitfall: Stay in ℤ for the row-decomposition step (use `zify`); the ℕ-saturating subtraction in `oGet A (i + r) - oGet A i` will silently bite you if you let it. Convert at the boundary, do the row identity in ℤ, convert back.

### Recipe: Inline sorted-enumeration helpers when imports cause v4.27 drift

The `Erdos30_OrderedElements.lean` file has a sorted-enumeration API (`orderedElements`, `orderedElement`, `orderedElements_sorted_lt`) but it has not been patched for v4.27 simp-lemma drift. When closing a downstream lemma that needs `oGet`-style access without wanting to widen the patch surface, prefer to inline a minimal local copy of the helpers (≤30 lines) into the consuming file rather than fix the upstream module.

Pattern: define `private noncomputable def oList (A : Finset ℕ) : List ℕ := A.sort (· ≤ ·)` and `private noncomputable def oGet (A : Finset ℕ) (i : ℕ) : ℕ := (oList A).getD i 0`. Prove the four facts you actually need (length, monotonicity, range-bound, membership) directly via `Finset.sort_sorted_lt` / `Finset.sort_perm` / `List.getD_eq_get`. Skip the high-API surface (`orderedElement_mem_take`, `intervalSlice` machinery) unless you actually consume it.

This costs ~30 lines per consumer; patching `OrderedElements.lean` for v4.27 would cost ~150 lines in a file you don't otherwise touch. The local-inline path also keeps the closure self-contained for review.

### Cross-cutting v4.27 lessons reinforced

- Tactics used: 22 `omega` + 6 `linarith`, all on legitimate Presburger / linear ℤ. No tactic-smuggling. `ring` for commutative-ring identities. `Finset.sum_range_id_mul_two` for the canonical Σr formula.
- `zify` is the right tool when ℕ-trunc subtraction and a pair-difference identity collide; do the algebra in ℤ, then push back.
- The `Finset.sum_range_sub` (in ℤ) lemma family is the workhorse for shifted-pair telescoping. Don't try to use the ℕ form — the saturation breaks the identity.

### Pitfall reaffirmed (from earlier 2026-05-02 entry)

The recipe earlier this session said "do not import `Erdos30_Complete`" — that was correct ONLY when Complete was broken. After today's repair, importing Complete is safe and is the canonical bridge. The lesson generalizes: recipes that say "do not import X" should be timestamped and re-validated whenever X is touched. Mark stale guidance with the date of the next-known repair.

## 2026-05-05 - PMF Tier-1 salvo recipes (#1, #30, #166, #755)

Compiled sources:

- `erdos-experiments/Erdos30/lean/Erdos1_DistinctSubsetSums.lean`
- `erdos-experiments/Erdos30/lean/Erdos166_SumFree.lean`
- `erdos-experiments/Erdos30/lean/Erdos30_SharpDiff.lean`
- `erdos-experiments/Erdos30/lean/Erdos755_BhG.lean`
- `erdos-experiments/Erdos30/lean/Erdos755_B3G.lean`
- `erdos-experiments/Erdos30/lean/Erdos755_BhG_General.lean`
- `erdos-experiments/Erdos30/lean/Erdos755_DifferenceCount.lean`

Outcome: four reusable elementary counting patterns were added to the arsenal. These are not SOTA breakthroughs, but they are high-yield proof infrastructure for the Atlas because they cover powerset injection, shift-injection, strict-upper Sidon difference injection, and arbitrary-length tuple fiber counting.

### Recipe: powerset-injection bound for distinct subset sums

Use this when a property says the map `s : Finset ℕ ↦ s.sum id` is injective on `A.powerset`, and every element of A lies in `Finset.range (N + 1)`.

Pattern:

1. Define the property as `Set.InjOn (fun s : Finset ℕ => s.sum id) ↑A.powerset`.
2. Use `Finset.card_powerset A` to rewrite `A.powerset.card = 2 ^ A.card`.
3. Prove every subset sum lies in `Finset.range (A.card * N + 1)` by bounding each summand with `N` and using `Finset.sum_le_card_nsmul`.
4. Use `Finset.card_image_of_injOn hD` to transfer cardinality from powersets to image.
5. Finish with `Finset.card_le_card` on the image subset.

This gives the standard bound `2 ^ A.card ≤ A.card * N + 1`. The valuable Lean lesson is to keep the image target as `Finset.range (A.card * N + 1)`; do not introduce an interval unless you need lower bounds too.

### Recipe: shift-injection for sum-free density

Use this when A is sum-free in `[0, N]` and you want the elementary bound `2 * A.card ≤ N + 1`.

Pattern:

1. Split on `A = ∅`.
2. Let `a := A.max' hne`.
3. Define the reflected image `A.image (fun x => a - x)`.
4. Prove the reflection is injective on A; `omega` sees injectivity once it has `x ≤ a` and `y ≤ a` from `A.le_max'`.
5. Prove A and the reflected image are disjoint. If `y = a - x` lies in A and `x ∈ A`, then `x + y = a ∈ A`, contradicting sum-free.
6. Bound both A and its reflected image inside `Finset.range (a + 1)`, use `Finset.card_union_of_disjoint`, then `a ≤ N`.

Pitfall: the empty case should be `subst hemp; simp`, not `simp [hemp]; omega`. The latter can leave no goals and trip the tactic-state checker.

### Recipe: strict-upper Sidon difference injection

Use this when converting the Sidon sum uniqueness predicate into the sharp difference bound `A.card * (A.card - 1) ≤ 2 * N`.

Pattern:

1. Define `strictUpper A := (A ×ˢ A).filter (fun p => p.2 < p.1)`.
2. Prove `Set.InjOn (fun p : ℕ × ℕ => p.1 - p.2) ↑(strictUpper A)`.
3. The injection proof requires four order cases. In two "crossed" cases, arithmetic alone is not enough; apply the Sidon predicate with permuted variables, then close the contradiction with `omega`.
4. Use `Prod.ext` explicitly. In the two non-crossed cases the equality directions differ, so prefer `Prod.ext hb.symm ha` or `Prod.ext ha.symm hb` over `ext`.
5. Count the strict-upper set via `Finset.offDiag_card` plus image-by-swap symmetry.
6. Avoid fragile rewrites like `Nat.mul_sub_one` inside hypotheses; move the target to the form Lean actually has and let `omega` handle the final Nat arithmetic.

This is the PMF-aligned finite close-packing shadow for Sidon: it gives the factor-2 sharpening over the raw ordered-sum bound without changing the mathematical statement.

### Recipe: arbitrary h-fold fiber counting with `Fintype.piFinset`

Use this when a B_h[g]-type condition bounds the number of ordered h-tuples with a fixed sum.

Pattern:

1. Model h-tuples as `Fintype.piFinset (fun _ : Fin h => A)`.
2. Define:

```lean
abbrev IsBhGSet (A : Finset ℕ) (h g : ℕ) : Prop :=
  ∀ s : ℕ,
    ((Fintype.piFinset (fun _ : Fin h => A)).filter
      (fun f : Fin h → ℕ => ∑ i, f i = s)).card ≤ h.factorial * g
```

3. Rewrite the tuple count with `Fintype.card_piFinset_const`.
4. Unpack membership with `Fintype.mem_piFinset`.
5. Bound the sum range by `∑ i, f i ≤ ∑ _i : Fin h, N = h * N`, then target `Finset.range (h * N + 1)`.
6. Apply `Finset.card_eq_sum_card_fiberwise`, then sum the uniform fiber bound.
7. Finish the constant rearrangement with `ring`.

This gives `A.card ^ h ≤ h.factorial * g * (h * N + 1)` for all h. The key Lean lesson is that `Fintype.piFinset` avoids hand-written nested products (`((A ×ˢ A) ×ˢ A) ×ˢ ...`) once h is symbolic.

### Recipe: ordered off-diagonal B_2[g] count

Use this when the full ordered B_2[g] predicate is already available and an unordered triangle proof would be over-engineered.

Pattern:

1. Use `A.offDiag` rather than a bespoke unordered pair set.
2. `Finset.offDiag_card` gives `A.card * A.card - A.card`; convert to `A.card * (A.card - 1)` with a case split on `A.card`.
3. Partition the off-diagonal pairs by sum with `Finset.card_eq_sum_card_fiberwise`.
4. Each off-diagonal fiber is a subset of the full ordered fiber, so the B_2[g] hypothesis bounds it by `2 * g`.
5. Finish with the same finite range `Finset.range (2 * N + 1)`.

This avoids the brittle upper-triangle bijection proof. The mathematical price is that the result is the honest ordered bound, not an illicit Sidon-style difference-injectivity generalization for g > 1.

### Recipe: symmetry quotient via a swapped shadow fiber

Use this when the native combinatorial statement counts unordered B_2[g]
representations, but the existing sum-counting proof expects ordered pairs.

Pattern:

1. Define three fibers over the same sum `s`:
   - `orderedFiber A s`: all `(a,b) ∈ A × A` with `a+b=s`
   - `unorderedFiber A s`: the upper triangle `a ≤ b`
   - `reversedStrictFiber A s`: the strict lower triangle `b < a`
2. Prove the reversed strict fiber injects into the unordered fiber by swapping
   coordinates: `(a,b) ↦ (b,a)`.
3. Prove the ordered fiber is covered by the union of the unordered fiber and
   the reversed strict shadow fiber using `by_cases hle : p.1 ≤ p.2`.
4. Combine `Finset.card_le_card_of_injOn`, `Finset.card_union_le`, and the
   native unordered hypothesis to get the ordered bound `≤ 2 * g`.
5. Feed that theorem into the existing ordered `b2g_sum_count`.

This lands at `Erdos755_SymmetryQuotient.unordered_to_ordered` and
`b2g_sum_count_unordered`. The useful lesson is that the missing term is not a
new extremal estimate; it is an orbit/shadow accounting term. For h=2, the
quotient is just upper triangle plus swapped strict lower triangle. General h
should be attacked later through multisets/orbit stabilizers, not by
hand-writing all permutations.

### Recipe: theorem-target algebra scaffold with explicit analytic axioms

Use this when a hard analytic or interval-certified theorem is not yet closed,
but the downstream closure algebra should be type-checked in Lean.

Pattern:

1. Define the local quotient type and all constants as transparent `def`s.
2. Define component proposition packages, e.g. `RadialCertificate14`,
   `ShapeConeCertificate14`, and `MixedRemainderAbsorption14`.
3. Prove only the algebraic splice theorem with `linarith`/`nlinarith`.
4. Keep the real analytic theorem as an explicit `axiom`, not a hidden `sorry`.
5. Add a separate weaker theorem-target proposition when experiments demote the
   strong route, e.g. `EpsScaledDeficit14`, so future agents attack the current
   bottleneck rather than a falsified stronger formulation.

This landed in `Ehp114LocalMixedRemainderScratch`. Build PASS means the theorem
target and algebraic dependencies are well typed. It does not promote the
analytic claim to COMPILED theorem status.

### Recipe: operator-norm residual perturbation

Use this when a numerical or finite-dimensional certificate has an unquantized
state `a`, a perturbed state `q`, a target `y`, and a bounded linear operator
`A`.

Pattern:

1. Rewrite the perturbed residual as
   `A q - y = (A a - y) + A (q - a)` using `rw [map_sub]` and `abel`.
2. Apply `norm_add_le` to split the old residual from the perturbation term.
3. Apply `ContinuousLinearMap.le_opNorm A (q - a)` to bound the perturbation by
   `||A|| * ||q - a||`.
4. If an external coefficient-error estimate is available, feed it through
   `mul_le_mul_of_nonneg_left` and `norm_nonneg A`.

This landed in `ErdosLean4V427.RH.BeurlingNymanMDL`. The compiled theorem is
finite Hilbert-space stability infrastructure for the Beurling-Nyman MDL
packet. It is not RH evidence and does not formalize zeta, fractional-part
basis functions, or componentwise rounding.

### Recipe: finite Euclidean coordinate-to-norm bound

Use this when a finite certificate has `k` active real coordinates and each
coordinate error is bounded by the same scalar.

Pattern:

1. State the vector in `EuclideanSpace ℝ (Fin k)`, not plain `Fin k → ℝ`; the
   former carries the Hilbert/L2 norm.
2. Rewrite the norm with `EuclideanSpace.norm_eq`.
3. Bound each squared coordinate by the scalar square using the coordinate
   hypothesis, `norm_nonneg`, and `nlinarith`.
4. Sum the coordinate bounds with `Finset.sum_le_sum`, then simplify the
   constant sum with `Fintype.card_fin` and `nsmul_eq_mul`.
5. Apply `Real.sqrt_le_sqrt`, then rewrite
   `sqrt((k : ℝ) * (Delta / 2)^2)` using `Real.sqrt_mul`, `Nat.cast_nonneg`,
   and `Real.sqrt_sq`.

This landed as `Erdos.RH_MDL.finite_coordinate_error_bound` and composes with
the operator-norm perturbation recipe as
`Erdos.RH_MDL.residual_bound_of_coordinate_error`.

### Recipe: certified predicate wrapper around a raw hypothesis

Use this when a theorem is already proved from a raw pointwise hypothesis, but
the paper or downstream API needs a named predicate as the stable interface.

Pattern:

1. Define a transparent predicate whose body is exactly the raw hypothesis,
   e.g. `CoordinatewiseQuantized k q a Delta := ∀ i, ...`.
2. Prove the paper-facing theorem by directly applying the raw-hypothesis
   theorem to the predicate assumption.
3. Do not put algorithmic content into the predicate unless it is already
   formalized; keep "certified output" separate from "how the certificate is
   produced."

This landed as `Erdos.RH_MDL.CoordinatewiseQuantized` and
`Erdos.RH_MDL.residual_bound_of_coordinatewise_quantized`.

### Recipe: orthogonal projection residual accounting

Use this when a finite Hilbert-space diagnostic splits a target vector into a
dictionary span plus an orthogonal residual.

Pattern:

1. State the dictionary block as a `Submodule ℝ E` with
   `[K.HasOrthogonalProjection]`.
2. Bound the residual norm by rewriting
   `y - K.starProjection y` as the orthogonal-complement projection via
   `Submodule.starProjection_orthogonal_val`.
3. Apply `Kᗮ.norm_starProjection_apply_le y` to get the
   norm-nonincreasing projection certificate.
4. For accounting identities, use
   `(K.starProjection_add_starProjection_orthogonal y).symm` to write
   `y = K.starProjection y + Kᗮ.starProjection y`.

This landed as
`Erdos.RH_MDL.orthogonal_projection_residual_norm_le` and
`Erdos.RH_MDL.orthogonal_projection_decomposition`. It is finite Hilbert-space
projection infrastructure only; it does not define primes, zeta, or RH.

## 2026-05-16 — Erdos1038 slack-augmented exchange (synthetic-infeasibility lesson)

**Source packet:** `EXP-MATH-ERDOS1038-SLACK-AUGMENTED-EXCHANGE-OBJECTIVE-20260516-01`
(satellite-resident at `Math/Math-Problems/Erdos-Standard/erdos-1038/`)

**Source skeleton code:** `Math/Math-Problems/Erdos-Standard/erdos1038-fast/{slack_lp_runner.py, lp_matrix_regen.py, run_slack_augmented_exchange.py}`

**Source Lean carriers (both COMPILED, Mathlib v4.24.0):**
- `Math/Lean4/private-mdl-workshop/Lean/Erdos1038_InfimumMeasure.lean`
- `Math/Lean4/private-mdl-workshop/Lean/Erdos1038_StrictComplementarity.lean`

### Recipe: do not promote a synthetic LP-input run to evidence on the certificate route

Use this when a Route Confound diagnosis names "regenerate the LP input
generator" as the live blocker, and you are tempted to validate an LP solver
end-to-end with a synthetic panel-constraint builder before the
interval-certified verifier panels are wired in.

Pattern:

1. The synthetic panel-constraint builder (`build_panel_constraints_synthetic`)
   linearizes `log|p(x)|` constraints on a uniform parameter grid without
   interval arithmetic. It is acceptable for *data-flow validation* of the LP
   solver wrapper, NOT for any certificate-bearing margin.
2. If you hard-floor the inactive-slack target (`s_I ≥ 1e-6`) against the
   synthetic LP, **expect HiGHS to report INFEASIBLE** — the synthetic LP
   simply does not have the structure that yields strict inactive slack at the
   real KKT point. This is a *negative control on data flow*, not a negative
   result on the math.
3. Tag the artifact with `synthetic: True` in JSON, prepend the REPORT.md with
   a ⚠️ SYNTHETIC banner, and write a separate SHA-256 for the amended
   REPORT.md before any downstream consumer touches it.
4. The real run is gated on integrating the Rust verifier's panel-extraction
   API into the LP-matrix regen step. Do NOT promote the synthetic 0/21 result
   to a Route Confound update; the Route Confound stays the same.

### Mathlib / repo pitfall

The cleanest sibling-file import pattern in `private-mdl-workshop/Lean/` is
the bare module name (e.g., `import Erdos1038_SublevelMeasure`), not
`import Lean.Erdos1038_SublevelMeasure`. The package's `srcDir := "."` means
the Lake module namespace has no `Lean.` prefix. Each new file must also be
registered in `lakefile.lean` as its own `lean_lib` with explicit `roots`,
or `lake build <Name>` will succeed-by-default on every sibling and silently
skip the new file.

### Reusable hypothesis

The slack-augmented LP pattern (objective penalty on inactive-slack indices
plus equality-form slack variables) is a valid generic primitive for any LP
where the underlying KKT point has tolerance-floor strict slack failures.
The pattern transfers as a *method prior*, not a structural law.

### Non-transfer warning

The pattern does NOT apply when:
- The blocker is not a degenerate KKT point but a true unbounded inactive
  alternative (then slack penalty cannot help — the LP is genuinely unbounded
  on the inactive face).
- The constraint matrix is rank-deficient at the optimum (then strict
  complementarity is structurally absent; no penalty α fixes this).
- The "slack" is measured in a non-standard norm or with rounding that
  obscures the floor (interval LP is required, not f64).

### Next falsifiable gate

`EXP-MATH-ERDOS1038-SLACK-AUGMENTED-EXCHANGE-OBJECTIVE-20260516-02` (next
session): replace `build_panel_constraints_synthetic` with the real Rust
verifier's panel-extraction API. Acceptance: 21/21 stress cases pass strict
inactive slack ≥ 1e-6 against interval-certified constraints. If 0/21 still,
escalate to chart-split + diagnostic per route-map kill criteria.

### Morphism classification

- Class: `route hypothesis` (the slack-augmented exchange route, not yet
  triangulated).
- Evidence tag: `[SKETCH]` (skeleton code + Lean axiom carriers, no closed
  certificate).
- This morphism has NOT transferred to another problem yet; do not promote
  to `FINITE_CERT_RECIPE` until 21/21 strict slack closes against real
  constraints AND a sibling Erdős problem with similar KKT-degeneracy
  structure inherits the slack-augmented exchange objective successfully.

### Compiled landings this session

- `Erdos1038_StrictComplementarity.lean` (CLEAN ✓ + COMPILED ✓, 0 sorries,
  2 axioms intentional, numerical-anchor lemma `threshold_exceeds_current_blocker`
  proven by `norm_num`).
- `Erdos1038_InfimumMeasure.lean` (CLEAN ✓ + COMPILED ✓, 0 sorries after
  closing `sublevel_measure_n1_eq_two` directly via `abs_lt` + `Real.volume_Ioo`).
  Seven axiomatic gates record the 8400m–8848m route obligations.
- v4.27.0 mirror at `MendozaLab.Erdos1038.{SublevelMeasure, InfimumMeasure, StrictComplementarity}`
  (DeepMind formal-conjectures alignment). Zero API drift between Mathlib
  v4.24.0 and v4.27.0 for these files.

### Recipe: SIGN-CONVENTION SANITY CHECK (the kind of bug that produces fake PASS results)

Use this before trusting any LP "21/21 PASS" result whose constraint matrix
comes from a different convention than scipy.linprog expects.

scipy.linprog uses `A_ub · x ≤ b_ub`. If you're encoding "panel positivity"
`Σ w · log|x - r_i| ≥ 0`, the matrix passed to scipy must be `-A_rows`
(negated), not `A_rows` itself. The synthetic builder gets this right by
construction (`A_ineq[p, i] = -np.log(v)`). A Rust binding that returns
A_rows in panel-positivity convention must be NEGATED by the Python
wrapper before scipy sees it.

The classic failure mode (caught 2026-05-16 in erdos1038 -02 and -03 runs):
the wrapper passes A_rows un-negated; scipy interprets the LP backwards;
the LP optimum chooses weights that VIOLATE the constraint by exactly
slack_floor (or by the natural-slack value). The "panel_pos=True" check in
`slack_lp_runner` is `s_opt ≥ -1e-12`, which is a tautology and doesn't
catch this. Required sanity check before trusting any LP run:

```
margin_at_baseline = A_rows @ w_baseline
margin_at_lp_opt   = A_rows @ w_lp_opt
```

If `margin_at_baseline >= 0` (panel positivity holds at baseline weights)
but `margin_at_lp_opt` has the opposite sign (negative values), the sign
convention is wrong somewhere. -02 and -03 retracted on this basis.

### Recipe: max-min LP for natural strict-inactive-slack diagnostics

When you want to test "is there a feasible weight vector with strict
positive panel-positivity margin on every panel," DON'T use the equality-
form slack-variable approach with sign-flipped objective penalty — it
hits scipy infeasibility at the boundary panel (machine-zero baseline
slack) and the natural-slack mode runs unbounded.

DO use the max-min LP via auxiliary variable t:

```
Variables: w ∈ ℝ^n, t ∈ ℝ
Constraints:
    A_rows · w ≥ t   (each panel)        ⇔ -A_rows · w + t ≤ 0
    Σ w_i = 1
    w_i ≥ 0
    t free
Objective: maximize t  (i.e. minimize -t)
```

`t*` (the LP optimum) is the BEST achievable worst-case panel-positivity
margin. Compare to threshold (e.g. 1e-6) for gate verdict.

This landed as `EXP-MATH-ERDOS1038-SLACK-AUGMENTED-EXCHANGE-OBJECTIVE-20260516-04`:
21/21 PASS_NATURAL_DOMINATION on the N=100 atom-cloud baseline, min t* =
9.468e-02 across stress cases. Constraint matrix interval-certified
(Rust + inari); LP solve f64 (scipy HiGHS).

### Recipe: feasibility-at-floor is not natural domination (LP audit pitfall)

Use this when a slack-augmented LP returns `min_inactive_slack = slack_floor`
uniformly across stress cases. The naïve reading is "gate closed: 21/21 PASS";
the correct reading is "the LP is sitting on the hard `s_inactive ≥ slack_floor`
bound, which is feasibility but NOT natural domination."

Pattern:

1. The slack-augmented LP has TWO mechanisms that push inactive slack up:
   (a) the `s_inactive ≥ slack_floor` bound (a hard constraint); and
   (b) the objective penalty `−α · sum(s_inactive)` (a soft preference).
   With (a) PRESENT, the LP optimizer prefers `s_inactive = slack_floor` exactly
   (because the penalty term wants `s_inactive` minimized subject to the
   bound). Returns `min_inactive_slack = slack_floor` exactly.
2. This DOES prove the LP is FEASIBLE at the floor — i.e., there exists at
   least one weight vector with inactive slack ≥ floor. That's a real
   existence claim. But it does NOT prove the OPTIMUM has inactive slack
   exceeding the floor naturally.
3. To test natural domination: REMOVE the hard bound (set `s_inactive ≥ 0`,
   keep the objective penalty), re-solve. If the LP optimum has
   `s_inactive ≥ slack_floor` without the bound being binding, that IS
   strict-complementarity evidence. If the optimum has `s_inactive ≈ 0`,
   the original blocker is confirmed and ordinary cuts + objective penalty
   are insufficient.
4. Per `/erdos-solve` v1.3.5 ("accept only if inactive alternatives are
   dominated by an interval-positive margin"), only the natural-domination
   reading counts as a 8400m gate close.

### Mathlib / repo pitfall (also)

`inari = "2.0"` is current stable on crates.io despite older docs citing
0.11. Empirical build success on stable Rust (no nightly needed for this
crate's ln/abs/arithmetic in default-features mode) trumps stale research.
Always run `cargo build --release` to verify before believing a research
finding about Rust crate stability.

### Compiled crate this session

- `erdos1038-fast` (Rust + PyO3 + rayon + inari 2.0): new
  `extract_panel_constraints` PyO3 function exposes interval-certified
  panel-positivity constraint matrices via IEEE 1788 `(panel − r).abs().ln()`.
  cargo test 6/6 ✔, maturin develop ✔. Smoke-test on N=100 atom-cloud:
  32 panels on `(z_left − 1, z_left) ∪ (z_right, z_right + 1)`, all panels on
  the legal-support side, min baseline slack at machine-zero (the boundary
  panel is at the active edge as expected).
