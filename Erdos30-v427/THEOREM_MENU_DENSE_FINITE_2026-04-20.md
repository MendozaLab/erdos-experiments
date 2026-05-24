# Theorem Menu for Erdős #30 — Dense Finite Rigidity

**Date:** 2026-04-20
**Scope:** Finite Sidon sets near the upper envelope
**Purpose:** Convert the current `#30` attack surface into a small set of lemma-level branches that match both the latest Sidon literature and the local formalization package.

**Status correction (later on 2026-04-20):** after checking Balasubramanian–Dutta, the current honest imported theorem interface is an ordered-element estimate at scale `O(n^{7/8}) + O(L^{1/2} n^{3/4})`. Candidate Lemma 1 below remains an aspirational local theorem target, not a theorem we are justified in importing at `sqrt(n)`-scale discrepancy.

**Computation correction (2026-04-22):** the exact scan
`EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22` shows that for every scanned
`10 ≤ n ≤ 50`, the true extremal size `h(n)` is strictly larger than
`floor(sqrt(n))`. So the current Lean package `DenseSidonAtScale A n L` is
studying a genuine dense-but-subextremal corridor, not the full maximizer
regime. That does not invalidate the rigidity lane, but it does mean future
local theorem statements must say explicitly which regime they live in.

**Maximizer check (later on 2026-04-22):** the exact scan
`EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22` then tested the widened
prefix and mass templates directly on all exact maximizers `A ⊆ [0,n]` for
`10 ≤ n ≤ 50`. Across 76,368 maximizers, the worst prefix residual after
subtracting the combined drift `max(||A| - sqrt(n)|, 1) · sqrt(n)` stayed
below `0.5304 · n^(7/8)`. So the new maximizer-friendly prefix wrapper is
pointed at the right regime. The mass center also stays bounded at the
literature-scale normalization, but with a larger observed constant
(`0.9505 · n^(11/8)`), so that side currently looks structurally coarser than
the prefix package.

**Calibration refinement (2026-04-22 and 2026-04-23):** the exact maximizer
packets now split the centering question cleanly. On the same `10 ≤ n ≤ 50`
window, even the raw prefix discrepancy
`|t - |A ∩ [0,t]| · sqrt(n)|`, with no drift subtracted at all, stayed below
`0.8928 · n^(7/8)`. But the first simple affine endpoint recentering candidate,
tested in `EXP-MM-030-RECENTERED-MAXIMIZER-PREFIX-DIAGNOSTICS-2026-04-23`, did
not improve the worst-case maximizer constant: it stayed at
`0.8856 · n^(7/8)`, much worse than the current theorem-aligned drift
(`0.5303 · n^(7/8)`). And the stronger density-adjusted affine prefix test,
recorded in `EXP-MM-030-DENSITY-ADJUSTED-MAXIMIZER-DIAGNOSTICS-2026-04-23`,
was worse still at `1.5981 · n^(7/8)`. So the empirical picture is now: prefix
slack is not removed by naive affine reanchoring, either by endpoint bridge or
by slope `n / |A|`.

## Why this branch

The old framing for `#30` was "search harder for better constructions" or "look for a new heuristic." That is no longer the highest-value reading of the problem.

The finite upper-bound frontier has moved in two directions:

1. The finite `n^{1/4}` coefficient has improved from the BFR `0.998` regime to `0.99703` (O'Bryant 2024) and then to `0.98183` (Carter–Hunter–O'Bryant 2025), with the strongest gain coming from computer-assisted combinatorics rather than a short new paper proof.
2. Dense finite Sidon sets appear to have much more internal rigidity than a naive "sparse random-looking object" model suggests. In particular, the `m`-th element and the sum of elements are both forced close to deterministic density-corrected templates.

That is exactly the geometry-constrained reading: the interesting objects are not arbitrary Sidon sets, but near-extremizers whose geometry is squeezed into a narrow corridor. This matches the physics intuition better than another unconstrained search.

The new exact scan sharpens that sentence. Right now the local formal corridor is
best read as "dense below the `sqrt(n)` floor" rather than "all near-maximizers."
Inside that corridor, the affine profile and mass center fit the exact data
well; across the full extremal problem, there is still a regime mismatch that
future statements need to absorb.

## Local package alignment

The current local package in `Erdos30/` already gives a clean base for this program:

- `Erdos30_Sidon_Defs.lean` supplies the canonical `IsSidonSet` definition.
- `Erdos30_Lindstrom.lean` and `Erdos30_BFR.lean` formalize the classical finite upper-bound chain.
- `PROOF_RECIPES.md` now records the reusable counting and shifted-family patterns extracted from the latest BFR-support batch.

What is missing is not more elementary counting. What is missing is a stability layer: lemmas that say a Sidon set close to the upper envelope must look almost evenly distributed, almost linearly spaced, and almost mass-balanced.

## Candidate Lemma 1 — Prefix Discrepancy Stability

**Status:** Aspirational local theorem target. Not currently justified by an imported literature theorem at `sqrt(n)`-scale discrepancy.

**Working statement**

Let `A ⊆ [n]` be Sidon with `|A| = floor(sqrt(n)) - L`, where `L` is small on the natural `n^{1/4}` or slightly larger scale. Then every initial segment `[1,t]` of macroscopic length contains

`|A ∩ [1,t]| = (t/n) |A| + error(t, n, L)`

with an error term substantially smaller than the main term.

**What it would mean**

This says dense Sidon sets cannot front-load or starve long prefixes. If a set is close to extremal, then its cumulative distribution function is already rigid on coarse scales. In physics language: the geometry is constrained at the coarse-grained level before any finer statement about individual points can even be true.

**Why this matters**

This is the gateway lemma we actually need. Once prefix discrepancy is controlled, the `m`-th element and sum-of-elements statements become natural corollaries rather than separate miracles.

**Execution mode:** `Lean-first`

This is the most natural formal next step because it fits the current package and it is strictly weaker than all-interval control:

- finite counting,
- prefix counting,
- Cauchy-Schwarz / discrepancy style arguments,
- no brute-force search required to state the result.

**Likely proof ingredients**

- BFR-style sum counting,
- discrepancy over prefixes,
- reuse of the `distinctSums_subset_range` / `card_le_card` pattern,
- eventually a formal interface for a literature theorem if the sharp constant is imported rather than reproved.

## Candidate Lemma 2 — Ordered Element Rigidity

**Status:** External interface, Balasubramanian–Dutta 2025/2026. Perplexity gate passed on 2026-04-20. Current honest imported scale is `O(n^{7/8}) + O(L^{1/2} n^{3/4})`, not a sharper local discrepancy corollary.

**Current formal state:** the external interface now splits into two honest layers. The older layer is the floor-`sqrt(n)` corridor package, which already yields compiled Lean consequences at seven levels: cutpoints `t = a_i`, nearby prefixes bracketed between consecutive ordered elements, an index-free nonterminal-prefix theorem for all `t` with `|A ∩ [0,t]| < |A|`, a positive-card all-prefix theorem for `t ≤ n` with the honest terminal drift `((L+1) : ℝ) · √n`, an ordered mass-balance theorem obtained by summing the external ordered-element estimate, the corresponding finset-mass theorem in the natural `∑ a ∈ A` language, and an explicit arithmetic-center rewrite of that mass theorem. The newer layer is the maximizer-friendly interface parameterized by `max(0, sqrt n - |A|)`, which now gives a general cutpoint theorem, a general positive-card all-prefix theorem with combined drift `max(|(A.card : ℝ) - sqrt(n)|, 1) · sqrt(n)`, a compiled super-floor refinement `sidon_in_range_superfloor_prefix_external` that first replaces this set-dependent endpoint term by the pure ambient correction `((Nat.sqrt (Nat.sqrt n) : ℝ) + 1) · sqrt(n)` whenever `floor(sqrt(n)) ≤ |A|`, then a coarse theorem `sidon_in_range_superfloor_prefix_coarse_external'` at the single `n^(7/8)` scale, and now also a maximizer-friendly interval-rigidity layer through `sidon_in_range_index_difference_external`, `sidon_in_range_cutpoint_interval_external`, and the coarse super-floor theorem `sidon_in_range_superfloor_cutpoint_interval_coarse_external`. That cutpoint layer has now been pushed all the way through the boundary split: `sidon_in_range_superfloor_internal_prefix_coarse_external`, `sidon_in_range_superfloor_empty_prefix_coarse_external`, `sidon_in_range_superfloor_nonterminal_prefix_coarse_external`, and `sidon_in_range_superfloor_terminal_prefix_coarse_external` together isolate the bulk, left edge, nonterminal package, and right edge, all at the same single `n^(7/8)` scale. The mass side is now packaged honestly as an envelope theorem too: `sidon_in_range_superfloor_finset_mass_coarse_external` gives a single `n^(11/8)` bound on the super-floor corridor, and `sidon_in_range_superfloor_finset_mass_density_adjusted_external` rewrites that same envelope around the density-adjusted center `n (|A| + 1) / 2`. That is still the right coarse scale for mass, but it now matches the exact April 23 split: density adjustment helps on mass while the prefix side still resists naive affine recentering. What remains open locally is no longer a missing case split. It is whether the density-adjusted mass constant can be tightened and whether prefix needs a genuinely different structural center rather than more endpoint bookkeeping.

**Small-`n` computation check:** in the exact floor-`sqrt(n)` scan
`EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22`, the largest observed prefix
residual beyond the deterministic `sqrt(n)` step stayed below
`0.7799 · n^(7/8)` and the largest observed mass deviation stayed below
`0.6409 · n^(11/8)`. That is not a proof of the literature-scale theorem, but
it is evidence that the current affine and mass templates are organizing the
data on the regime they are actually meant to describe.

**Maximizer computation check:** in the exact maximizer scan
`EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22`, the widened general prefix
wrapper also survives contact with the true extremal regime. The largest
observed prefix residual after subtracting the theorem-aligned drift
`max(||A| - sqrt(n)|, 1) · sqrt(n)` stayed below `0.5304 · n^(7/8)`, while the
largest observed mass deviation stayed below `0.9505 · n^(11/8)`. The right
reading is asymmetrical: the general prefix theorem now looks empirically well
centered on actual maximizers, while the explicit mass center is honest but not
obviously close to sharp.

The more surprising calibration fact is now asymmetric. For prefixes, the raw
discrepancy already sits at `0.8928 · n^(7/8)`, the endpoint-bridge recentering
stays at `0.8856 · n^(7/8)`, and the density-adjusted affine prefix test jumps
to `1.5981 · n^(7/8)`. So the next prefix target is no longer naturally "find a
bigger endpoint correction," but it is also not "recenter to slope `n / |A|`."
For mass, though, the density-adjusted center is a real improvement: the old
`sqrt(n)`-centered mass template peaks at `0.9505 · n^(11/8)`, while the
density-adjusted mass center drops to `0.5904 · n^(11/8)`. So the live theorem
split is now clean: density adjustment looks promising for mass, not for
prefix.

The newest exact packet sharpens that again at the observable level. In
`EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23`, the
prefix-best maximizer and the density-adjusted-mass-best maximizer coincide in
only `8` of the `41` scanned values of `n`. So the finite data is no longer
suggesting just "one bad affine center." It is now plausibly suggesting that
prefix and mass are genuinely different observables on the same extremal
objects, and that a future theorem may need different coordinates for those two
readouts rather than one global recentering.

The refined finite read is asymmetric but not monotone. The split stays common
throughout the window rather than obviously strengthening with `n`: it appears
in `16/21` cases for `10 ≤ n ≤ 30` and `17/20` cases for `31 ≤ n ≤ 50`. The
more robust signal is one-sided Pareto behavior. In the exact packet, the
mass-best witness has numerically zero prefix residual in `31/41` cases, while the
prefix-best witness has zero density-adjusted mass deviation in only `1/41`
case. So the safest current summary is "persistent mixed split with strong
one-sided flavor," not "cleanly widening large-n regime separation."

The KvN/Holevo analogy is useful only as a search lens here. The Sidon set is
the state, prefix and density-adjusted mass are readouts, and the finite packet
is warning that one readout does not expose the whole structure. The disciplined
translation is not "Sidon sets obey Holevo." It is: look for a compatibility
law between observables, and test it first against exact maximizers before
trying to formalize it.

The first finite compatibility candidate is one-sided. In the April 24
compatibility-summary rerun, the mass-best witness has prefix cost at most
`0.0884 · n^(7/8)` across the whole window, while the prefix-best witness has
density-adjusted mass cost as high as `0.2025 · n^(11/8)`. In normalized units,
the mass-best prefix penalty is no larger than the prefix-best density-adjusted
mass penalty in `40/41` values of `n`, with the only exception a numerical tie
at `n = 14`.

The first bounded beyond-50 probe keeps the same shape. In
`EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`, exact
maximizers for `51 ≤ n ≤ 55` satisfy the direct normalized comparison in `5/5`
values; the joint witness is mass-best in `4/5`. This is still structural
evidence, not a SOTA bound improvement.

The Rust exact-maximizer extension pushes that same test to `56 ≤ n ≤ 60`.
There the direct comparison holds in `4/5` values, with a real exception at
`n = 57`; the joint witness is mass-best in `3/5`, prefix-best in `2/5`, and a
third witness at `n = 59`. So the theorem menu should keep the exact-surface
compatibility seed alive, but not upgrade it to a clean monotone law.

The next Rust packet, `61 ≤ n ≤ 65`, returns to `5/5` on the direct comparison.
That made `n = 57` look isolated at first, but not irrelevant: it was the
first local pinch the seed had to explain. The joint witness is still mixed,
landing on a third witness at `n = 63` and `n = 65`.

The combined `51 ≤ n ≤ 65` exact-surface roll-up is `14/15` for the direct
comparison across `49,708` exact maximizers. That is the strongest current
reason to keep this theorem seed on the menu, and also the reason not to call it
monotone.

The `66 ≤ n ≤ 68` Rust boundary probe keeps the direct comparison in `3/3`,
making the `51 ≤ n ≤ 68` roll-up `17/18` across `115,354` exact maximizers.
That strengthens the seed, while the lone `n = 57` exception still blocks
monotone language.

The instrumented `n = 69` packet keeps the direct comparison again, making the
`51 ≤ n ≤ 69` roll-up `18/19` across `181,766` exact maximizers. The joint
witness is still a third witness at `n = 69`, so the menu item remains
one-sided compatibility, not optimizer selection.

The first-hit `n = 70` packet adds a second direct-comparison failure. The
current first-hit `51 ≤ n ≤ 70` roll-up is `18/20` across `298,968` exact
maximizers, with failures at `n = 57` and `n = 70`. The follow-up top-k
frontier packet corrects the interpretation: `n = 70` has joint witnesses with
zero prefix residual and zero density-adjusted mass deviation, so it is a
tie-selection artifact, not a face-level incompatibility. That keeps the menu
item alive, but changes the intended theorem shape from "almost monotone" to
"face-aware exact-surface compatibility."

`PINCH_ANALYSIS_57_70_2026-04-29.md` is now the active diagnostic for this
menu item. It reads `n = 57` as the remaining small witness-handoff pinch and
`n = 70` as a first-hit tie-selection artifact once top-k witnesses are
examined.

The near-maximizer pilot sharpens the boundary. In
`EXP-MM-030-NEAR-MAXIMIZER-COMPATIBILITY-PILOT-2026-04-24`, the same
compatibility pattern remains strong on the `h(n)` layer for `10 ≤ n ≤ 30`
(`20/21` values), but weakens sharply on the `h(n)-1` layer (`9/21` values).
So the current best reading is that the one-sided signal is attached to the
exact extremal surface, not yet to generic dense Sidon stability.

That boundary is now recorded as an explicit theorem seed in
`EXTREMAL_SURFACE_COMPATIBILITY_SEED_2026-04-24.md`: any future compatibility
lemma should condition on exact maximizers, or on a near-extremal predicate
strong enough to behave like the exact surface. The generic dense-Sidon version
is too broad for the current data.

The Lean scratch file now has the vocabulary needed to state that boundary:
`IsMaximalSidonInRange`, `NearExtremalSidonInRange`,
`prefixResidualAfterGeneralDrift`, and `densityAdjustedMassDeviation`. These
definitions compile, but they are not yet a compatibility theorem.
The first exact-surface wrapper,
`maximal_sidon_in_range_superfloor_prefix_mass_joint_envelope_external`, also
compiles and simply reuses the existing joint envelope under
`IsMaximalSidonInRange`.

The current Lean shadow of that story is intentionally weaker. The scratch file
now contains `sidon_in_range_superfloor_prefix_mass_joint_envelope_external`,
which packages the super-floor coarse prefix theorem and the density-adjusted
mass theorem into one statement about the same set. That is the nearest honest
formal analogue of the observable-split picture available right now: one object,
two observables, simultaneous envelopes. It is not yet a tradeoff theorem, not
an incompatibility law, and not a formal statement about distinct optimizer
selection.

**Working statement**

Let `A = {a_1 < ... < a_k} ⊆ [n]` be a dense Sidon set. Then the right linear
profile is

`a_m + 1 ≈ (m n) / k`

or, in division-free form,

`m n ≈ k (a_m + 1)`.

Since `k = sqrt(n) - L`, this is equivalent to

`m n ≈ (sqrt(n) - L) (a_m + 1)`.

**What it would mean**

This turns coarse occupancy into pointwise rigidity. A dense Sidon set is not merely "spread out"; its ordered elements lie near a deterministic linear profile with slope `n / |A|`. In the physics reading, this is the discrete analogue of a constrained ground-state profile: once the geometry is saturated, there is little freedom left in where the particles can sit.

**Why this matters**

This is the first theorem that starts to look classification-like. If true in a strong enough form, it says that any near-extremizer is trapped near an almost affine template with the correct density-adjusted slope. That sharply reduces the search space for construction experiments.

**Execution mode:** `Imported theorem interface first; Lean consequences after`

This is no longer the first local proof target. The honest state is that the ordered-element theorem is already available as an external theorem interface. The Lean work now is to extract local corollaries and connect it to the rest of the rigidity program.

**Likely proof ingredients**

- external ordered-element theorem interface,
- Candidate Lemma 1 if we still want a local discrepancy theorem,
- monotonicity from the ordered list of elements,
- a translation from interval count to inverse CDF / quantile control,
- careful separation between the true `n / |A|` slope and the naive `sqrt(n)` heuristic.

## Candidate Lemma 3 — Mass Balance / Sum-of-Elements Rigidity

**Working statement**

For a dense Sidon set `A ⊆ [n]`,

`sum_{a in A} a = (1/2) n^(3/2) + error(n, L)`.

**What it would mean**

This says that even the total mass of the configuration is almost forced. Once the set is close to maximal, it cannot only have the right density; it must place that density with nearly the right center of mass.

**Why this matters**

This is the cleanest bridge to a geometry-constrained physical picture. The set does not just satisfy a forbidden-sums rule. It pays a global placement cost. If the ordered elements are nearly linear and the total mass is nearly fixed, then the object behaves like a constrained low-energy configuration rather than a loose combinatorial gadget.

**Execution mode:** `Lean-first or computation-supported`

Once Lemma 2 exists, this becomes a sum of nearly linear positions. Without Lemma 2, it can still be tested computationally on exact dense Sidon sets and candidate near-extremizers.

**Likely proof ingredients**

- Candidate Lemma 2,
- summation of the ordered-element estimate,
- or direct prefix-discrepancy summation if that route is cleaner.

## Construction branch that depends on these lemmas

These lemmas do not directly solve `#30`. Their value is that they define a **restricted model class** for computation.

If the three rigidity statements are even approximately right, then construction search should stop ranging over all Sidon sets and instead range over:

- near-equispaced templates,
- perturbations of affine `(m n) / |A|` profiles,
- faithful algebraic constructions plus small structured perturbations.

That is the only construction branch that still looks worth serious compute. The old unconstrained branches already look weak:

- spectral greedy underperformed standard greedy in local experiments,
- toy Singer / toy Ruzsa proxies were not faithful enough to the actual algebraic mechanism,
- the morphism program currently classifies `#30` as a geometry-dominated generic control, not a scoped breakthrough case.

## Lean vs computation split

### Lean-first

- Candidate Lemma 1: prefix discrepancy stability
- consequences of the external ordered-element theorem at cutpoints and nearby prefixes
- Candidate Lemma 3: sum-of-elements rigidity

These are theorem statements about finite geometry and counting, which is where the current package already has traction.

### Computation-first

- exact dense-Sidon data generation,
- fitting the error term scales in Lemmas 1–3,
- searching only inside the rigid model class suggested by the lemmas,
- stress-testing whether extremizers really cluster near affine `(m n) / |A|` templates.

Computation should support conjecture sharpening, not replace the theorem program.

## Best next move

The shortest useful next step is:

1. Split the narrative explicitly into two regimes: the current floor-`sqrt(n)`
   dense corridor and the true extremal regime.
2. Keep Candidate Lemma 2 as the imported external theorem on the dense corridor.
3. Decide whether Candidate Lemma 1 should be reproved for that corridor, or
   rewritten against a different density parameter that can see maximizers.
4. Treat the exact `10 ≤ n ≤ 50` observable-split packet as a finite Pareto
   note, not as evidence that one universal affine center is about to emerge.

My current recommendation is still to treat **Candidate Lemma 2** as an imported external interface and keep mining honest consequences from it. The bridge to `∑ a ∈ A, a`, the explicit arithmetic-center rewrite, the super-floor coarse prefix theorem, the super-floor coarse cutpoint-interval theorem, the full coarse boundary split, and now the super-floor coarse mass theorem are all compiled. So the next local breakpoint is no longer missing structure in the theorem ladder. It is quantitative: whether the coarse constants can be tightened, and whether the mass side can be made less visibly looser than the prefix/interval side without claiming more than the exact maximizer packet actually supports.

## References driving this menu

- Kevin O'Bryant, *On the size of finite Sidon sets* (2024): finite `n^{1/4}` coefficient improved to `0.99703`.
- Daniel Carter, Zach Hunter, Kevin O'Bryant, *On the diameter of finite Sidon sets* (2025): finite `n^{1/4}` coefficient improved to `0.98183`, with substantial computer assistance.
- R. Balasubramanian, Sayan Dutta, *The m-th element of a Sidon set* (Journal of Number Theory, 2026 issue / 2025 DOI): ordered-element and sum-of-elements rigidity for dense Sidon sets.

## Non-goals

This note does **not** claim:

- a new bound for `h(n)`,
- a route to the infinite problem by itself,
- any scoped Leg-4 upgrade in the morphism program.

It is a theorem menu for the finite rigidity branch only.
