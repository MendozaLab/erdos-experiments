# Formal Verification of Sidon Set Upper Bounds in Lean 4

## Erdos Problem #30 — $1,000 Prize (OPEN)

**Question:** Is h(N) = N^{1/2} + O_epsilon(N^epsilon) for every epsilon > 0?

This repository contains a Lean 4 formalization of the two principal upper bounds
for Sidon sets (B_2 sets), covering the main known results on the upper bound side
of Erdos Problem #30.

The current research extension also includes a dense-finite rigidity lane in
`scratch/Erdos30_IntervalOccupancyTarget.lean`. That lane imports the published
Balasubramanian-Dutta ordered-element theorem and derives honest prefix and mass
consequences from it. It now has both a floor-`sqrt(n)` corridor interface and
a maximizer-friendly interface parameterized by `max(0, sqrt n - |A|)`. It also
has a compiled super-floor prefix refinement: once `floor(sqrt(n)) ≤ |A|`, the
set-dependent endpoint drift can be replaced by a pure ambient correction
`((Nat.sqrt (Nat.sqrt n) : ℝ) + 1) * sqrt(n)`, and then further collapsed into
a coarse single-scale `n^(7/8)` prefix theorem on the same corridor. It is a
rigidity interface, not a new bound for `h(N)`.

## What is formalized

| Result | Statement | File | Status |
|--------|-----------|------|--------|
| **Shared definition** | IsSidonSet (B₂ property) | `Erdos30_Sidon_Defs.lean` | Canonical definition |
| **Erdos-Turan (1941)** | k(k-1) <= 2N | `Erdos30_Complete.lean` | Fully proved, 0 sorry, 0 axioms |
| **Lindstrom parametric (1969)** | k^2 <= 2Nt + kt | `Erdos30_Lindstrom.lean` | Fully proved |
| **Lindstrom quadratic (1969)** | l(2k-l-1)^2 <= 4(l+1)N | `Erdos30_Lindstrom.lean` | Proved from axiom |
| **Lindstrom weak (1969)** | k <= sqrt(2N) + 1 | `Erdos30_Lindstrom.lean` | Fully proved |
| **Lindstrom full (1969)** | k <= sqrt(N) + N^{1/4} + 1 | `Erdos30_Lindstrom.lean` | Axiom (R->N step) |
| **BFR sum counting (2023)** | \|distinctSums\| = k(k+1)/2 | `Erdos30_BFR.lean` | Fully proved |
| **BFR Cauchy-Schwarz (2023)** | Variance decomposition | `Erdos30_BFR.lean` | Fully proved |
| **BFR bound (2023)** | 1000k <= 1000*sqrt(N) + 998*sqrt(sqrt(N)) + 1000 | `Erdos30_BFR.lean` | Proved from axiom |
| **Singer lower bound (1938)** | h(N) >= (1-o(1))*sqrt(N) | `Erdos30_Singer.lean` | Partial (q=2,3,5 verified) |

## Architecture

All files import `IsSidonSet` from the shared `Erdos30_Sidon_Defs.lean`, eliminating
duplicate definitions across the formalization. The three main files (Lindstrom, BFR,
Singer) compile against Mathlib via `lake build`.

## Verification status (last verified 2026-05-02 against Mathlib v4.27.0)

- **Zero sorry stubs** in all 5 main files (`Erdos30_Sidon_Defs.lean`, `Erdos30_Complete.lean`, `Erdos30_Lindstrom.lean`, `Erdos30_BFR.lean`, `Erdos30_Singer.lean`)
- **6 axioms** total in the main package, all with full bibliographic references — see [AXIOM_INVENTORY.md](AXIOM_INVENTORY.md) for the canonical catalog with closure plans:
  - 3 in Lindstrom: `order_diff_counting` (T2), `lindstrom_bound` (T2), and previously `sidon_elem_bound` (CLOSED 2026-05-02 via bridge to `sidon_difference_count`)
  - 1 in BFR: `bfr_core_bound` (T3)
  - 1 in Singer: `singer_sidon_exists` (T3)
  - 2 in the Erdős #755 sidecar (`Erdos755_*`): `singer_b2g_exists`, `lindstrom_sieve`
- **14+ theorems** fully machine-checked
- **Build:** `lake build Erdos30_Sidon_Defs Erdos30_Complete Erdos30_Lindstrom Erdos30_BFR Erdos30_Singer` PASSES (5 of 5 targets)

## Build

Requires: Lean 4 v4.24.0, Mathlib at commit `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`.

```bash
# From the lakefile.lean that imports Mathlib:
lake build Erdos30_Sidon_Defs   # shared IsSidonSet definition
lake build Erdos30_Lindstrom    # should produce 0 errors, 0 sorry
lake build Erdos30_BFR          # should produce 0 errors, 0 sorry
lake build Erdos30_Singer       # should produce 0 errors, 0 sorry
```

## Axiom inventory (canonical: see [AXIOM_INVENTORY.md](AXIOM_INVENTORY.md))

| Axiom | Statement | Closure tier |
|-------|-----------|--------------|
| ~~`sidon_elem_bound`~~ | ~~k(k-1) <= 2M~~ | ✅ **CLOSED 2026-05-02** — bridge to `Erdos30_Complete.sidon_difference_count` |
| `order_diff_counting` | Sorted diffs of orders 1..ℓ with bounded sum | T2 — `orderEmbOfFin` + telescoping (~100 lines, next round) |
| `lindstrom_bound` | k <= floor(sqrt(N)) + floor(N^{1/4}) + 1 | T2 — needs `Nat.sqrt` rounding lemma after `order_diff_counting` |
| `bfr_core_bound` | 1000(k-1) <= 1000*floor(sqrt(N)) + 998*floor(N^{1/4}) | T3 — full BFR §2–4 (~15 lemmas), research-tier |
| `singer_sidon_exists` | ∃ Sidon A ⊆ [0,q²+q] with \|A\|=q+1 for prime q | T3 — GaloisField + Singer cycle, novel formalization |


## Scratch / supplementary files

The `scratch/` directory contains exploratory files not part of the main package and not imported by any core file:

| File | Notes |
|------|-------|
| `Erdos30_difference_counting.lean` | TIER 3 — textbook k(k−1) ≤ 2N bound via difference counting. Duplicate of `Erdos30_Complete.lean`, kept for reference. |
| `Sidon_SumCount_Fix.lean` | Sum-counting helper generated during development. Referenced only in comments in BFR and Singer; not imported. |

These files compile independently but are excluded from the core package boundary.

## Exact computation artifact

The current dense-rigidity lane is cross-checked by the exact computation packet
`erdos-experiments/results/erdos-30/EXP-MM-030-DENSE-RIGIDITY-SMALLN-2026-04-22_*`.
That scan enumerates every Sidon set `A ⊆ [0,n]` with `|A| = floor(sqrt(n))`
for `10 ≤ n ≤ 50`, totaling 25,189,976 exact dense sets.

The main strategic lesson is scope: in every scanned case, the true extremal
size `h(n)` is strictly larger than `floor(sqrt(n))`. So the current Lean
package `DenseSidonAtScale` is studying a dense but sub-extremal corridor. The
same scan also shows that on this corridor the affine profile, prefix, and mass
templates from the external Balasubramanian-Dutta interface fit the exact data
cleanly enough to justify continuing the rigidity program. The scratch file now
also contains a wider prefix/mass package for `SidonInRange`, where the
terminal regime is handled by symmetric `|(A.card : ℝ) - sqrt(n)|` bookkeeping,
and the resulting all-prefix wrapper compresses to a single
`max(|(A.card : ℝ) - sqrt(n)|, 1) * sqrt(n)` drift term rather than the older
floor-deficiency bookkeeping.

That widened package is now cross-checked by a second exact packet,
`erdos-experiments/results/erdos-30/EXP-MM-030-MAXIMIZER-RIGIDITY-SMALLN-2026-04-22_*`,
which enumerates every exact maximizer `A ⊆ [0,n]` with `|A| = h(n)` for the
same window `10 ≤ n ≤ 50`, totaling 76,368 maximizers. The main lesson from
that run is that the new general prefix theorem is aimed at the right regime:
after subtracting the theorem-aligned drift
`max(|(A.card : ℝ) - sqrt(n)|, 1) * sqrt(n)`, the worst observed prefix
residual stayed below `0.5304 * n^(7/8)` in the scanned window. The explicit
mass center remains honest but looser, with worst observed deviation below
`0.9505 * n^(11/8)`. So the live open question is no longer scope mismatch; it
is whether the maximizer-friendly endpoint bookkeeping can be tightened. The
April 22 and April 23 exact packets now sharpen that sentence into a split. For
prefixes, even the raw discrepancy stayed below `0.8928 * n^(7/8)` on
`10 ≤ n ≤ 50`, but the first simple affine endpoint recentering candidate did
not help (`0.8856 * n^(7/8)`), and the stronger density-adjusted affine prefix
test was much worse (`1.5981 * n^(7/8)`). So the next prefix theorem target is
not naturally "add more endpoint drift," but it is also not "just recenter the
affine profile." For mass, the story goes the other way: the density-adjusted
mass center from
`EXP-MM-030-DENSITY-ADJUSTED-MAXIMIZER-DIAGNOSTICS-2026-04-23_*` lowers the
worst observed ratio from `0.9505 * n^(11/8)` to `0.5904 * n^(11/8)`. So the
right next target is now split: density-adjusted mass looks promising, while
prefix still wants a different structural idea.

The formal side has now taken the first honest absorption step as well. In the
same scratch file, `sidon_in_range_superfloor_prefix_external` shows that on the
super-floor corridor `floor(sqrt(n)) ≤ |A|`, the old set-dependent endpoint term
can be replaced by the pure ambient correction
`((Nat.sqrt (Nat.sqrt n) : ℝ) + 1) * sqrt(n)`. That does not yet give a
no-drift prefix theorem, but it means the next local target is genuinely about
absorbing a fourth-root ambient term into the `n^(7/8)` literature scale, not
about fixing a scope mismatch or a theorem stated against the wrong class of
sets.

That next absorption step is now compiled as well. The theorem
`sidon_in_range_superfloor_prefix_coarse_external` packages the entire
super-floor prefix bound at a single `n^(7/8)` scale, with no explicit
set-dependent drift and no separate fourth-root ambient correction term. This
is still a coarse theorem and still depends on the imported
Balasubramanian-Dutta interface, but it is the first point in the current lane
where the prefix statement itself looks structurally like a real discrepancy
theorem rather than a bookkeeping wrapper.

The cutpoint/interval layer is now in the same state. The scratch file contains
`sidon_in_range_index_difference_external`,
`sidon_in_range_cutpoint_interval_external`, and the super-floor coarse theorem
`sidon_in_range_superfloor_cutpoint_interval_coarse_external`, which together
say that on the super-floor corridor the displacement between ordered cutpoints
is controlled at the same single `n^(7/8)` scale. That upgrade has now been
pushed one step further as well: `sidon_in_range_superfloor_internal_prefix_coarse_external`
gives an index-free internal-prefix theorem at the same single `n^(7/8)` scale,
and the boundary cases are now isolated cleanly too through
`sidon_in_range_superfloor_empty_prefix_coarse_external`,
`sidon_in_range_superfloor_nonterminal_prefix_coarse_external`, and
`sidon_in_range_superfloor_terminal_prefix_coarse_external`. So the next honest
local question is no longer about prefix geometry. It is whether the mass side
can be tightened quantitatively. The scratch file now also contains
`sidon_in_range_superfloor_finset_mass_coarse_external`, so the mass side is
packaged at the honest single `n^(11/8)` scale on the same super-floor
corridor. That matches the exact maximizer packet much better than any
prefix-style rigidity claim would, but it remains visibly coarser than the
prefix and cutpoint-interval side. The new honest next step is now compiled
too: `sidon_in_range_superfloor_finset_mass_density_adjusted_external` rewrites
that same super-floor mass envelope around the density-adjusted center
`n * (|A| + 1) / 2`, still at the single `n^(11/8)` scale. So the formal side
now matches the April 23 exact split directly: density adjustment is a good
theorem move for mass, but not for prefix.

The newest exact packet then turns that split into a more structural one.
`EXP-MM-030-OBSERVABLE-SPLIT-MAXIMIZER-DIAGNOSTICS-2026-04-23_*` asks whether
the same exact maximizers optimize both the best prefix observable and the
density-adjusted mass observable. On `10 ≤ n ≤ 50`, the answer is mostly no:
the prefix-best and mass-best witnesses coincide for only `8` of the `41`
scanned values of `n`. So the best current reading of `#30` is no longer "find
one better affine center for everything." It is that prefix and mass may be
genuinely different observables on the same near-extremal Sidon object, and the
next theorem target should respect that split rather than trying to erase it by
better bookkeeping.

The April 29 PMF transfer-operator lane makes that split executable as a
lattice-gas state machine. The Rust crate
`erdos-experiments/Erdos30/rust-transfer-operator/` represents a state as an
occupied-site bitset plus used positive-difference memory, and the corrected
packet `EXP-MM-030-PMF-TRANSFER-PARITY-10-30-V3-2026-04-29_*` reproduces the
current exact Rust reference `h(n)` and maximizer counts for every `10 ≤ n ≤
30`. The field-frontier packet
`EXP-MM-030-PMF-TRANSFER-FIELD-FRONTIER-56-58-2026-04-29_*` also reproduces the
exact scanner's top-1 prefix, mass, and joint witnesses over the suspected
handoff window. The first 56-58 spectral probe
`EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29_*` says the suspected `57`
pinch is not a raw cardinality-gap singularity: the gap remains `1`, while
ground degeneracy and near-ground counts rise smoothly across `56 -> 57 -> 58`.
The pruned scale packets
`EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-69-2026-04-29_*`,
`EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-70-2026-04-29_*`, and
`EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29_*` then reproduce the
exact top-1 prefix, mass, and joint frontier witnesses at the current 69-71
frontier without enumerating low-cardinality terminal states. The live
interpretation is therefore narrower and better. The `d = 1` near-ground
packets for 69-71 retain tens of millions of first-excited states, and that
layer grows smoothly rather than spiking; the visible signal is the
zero-temperature field response on the exact maximizer face. Finite exact
maximizer faces show lattice-gas ground-state degeneracy with isolated
field-sensitive handoffs. It is not a proof of Sidon.

The first cross-problem control is now in place too. The sum-free binary
`rust-transfer-operator/src/bin/sumfree_transfer.rs` ports the same layer-pruned
state-machine API to a changed local rule, and
`EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-20-40-2026-04-29_*` plus
`EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-41-60-2026-04-29_*` match
`h(n)=ceil(n/2)` for all `20 ≤ n ≤ 60`. The frontier split count is `0`: prefix,
mass, and joint fields all choose the same odd-set witness. That is a useful
negative control. The API ports, but the Sidon handoff signature is not a
generic artifact of the machinery.

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

The first bounded beyond-50 probe preserved that exact-surface reading rather
than dissolving it. In
`EXP-MM-030-COMPATIBILITY-BEYOND50-MAXIMIZER-PROBE-2026-04-24`, the same direct
comparison holds for all `5/5` scanned values `51 ≤ n ≤ 55`; the joint witness
is mass-best in `4/5`. This does not improve the headline #30 bound, but it
does make the compatibility ridge less likely to be a `10 ≤ n ≤ 50` accident.

The Rust exact-maximizer port extends the stress test without changing the
packet contract. A replay of `51 ≤ n ≤ 55` matched the Python packet on checked
counts, witnesses, and compatibility diagnostics, reducing runtime from about
`193.3s` to about `9.93s`. The first Rust-only extension,
`EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24`, keeps the direct
comparison in `4/5` values, with the exception at `n = 57`. So the ridge
persists, but the local pinches are real.

The next Rust packet,
`EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-61-65-2026-04-24`, rebounds to `5/5`
on the direct comparison. At that point `n = 57` looked isolated rather than
the start of a failure band. The frontier was still mixed, though: the joint
witness is a third witness at `n = 63` and `n = 65`.

Combined across the bounded exact-surface probes `51 ≤ n ≤ 65`, the direct
comparison holds in `14/15` values across `49,708` exact maximizers, with only
`n = 57` failing. That is a strong finite ridge, not a headline SOTA result.

The bounded `66 ≤ n ≤ 68` Rust probe keeps the direct comparison in `3/3`,
bringing the `51 ≤ n ≤ 68` roll-up to `17/18` across `115,354` exact
maximizers. The only failure through `68` is `n = 57`. The cost boundary is now real:
`EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-66-68-2026-04-24` took about `76.9s`,
and `n = 68` has `36,234` maximizers.

The instrumented `n = 69` packet keeps the direct comparison too, bringing the
`51 ≤ n ≤ 69` roll-up to `18/19` across `181,766` exact maximizers. It also
shows why progress telemetry matters: `n = 69` has `66,412` maximizers and the
search visited about `1.52B` recursive nodes. Use `--progress` for any next
single-n packet, and treat `--seed-depth 2` as experimental because it preserved
parity but was slower on `66 ≤ n ≤ 68`.

The first-hit `n = 70` packet adds a second direct-comparison failure. It has
`117,202` maximizers, visits about `1.82B` recursive nodes, and fails the
single-witness comparison. The current first-hit `51 ≤ n ≤ 70` roll-up is
`18/20` across `298,968` exact maximizers, with failures at `n = 57` and
`n = 70`.

The pinch comparison is recorded in `PINCH_ANALYSIS_57_70_2026-04-29.md`.
The current top-k reading is that `n = 57` is the remaining small
witness-handoff pinch, while `n = 70` is a first-hit tie-selection artifact:
the exact face contains joint witnesses with zero prefix residual and zero
density-adjusted mass deviation.

Run the Rust scanner from
`erdos-experiments/Erdos30/rust-exact-maximizer/`:

```bash
cargo run --release -- --n-min 56 --n-max 60 \
  --experiment-id EXP-MM-030-RUST-COMPATIBILITY-MAXIMIZER-56-60-2026-04-24 \
  --output-dir /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30
```

For future exact-surface packets, enable bounded frontier witnesses with
`--frontier-k 5` or higher. The first-hit fields preserve compatibility with
older packets, but the `top_k_frontier` block is now the safer diagnostic for
distinguishing a genuine face-level pinch from tied optimizer selection.
New packets also emit `first_hit_vs_face_aware_diagnostics`, and the Markdown
report prints the first-hit and face-aware verdicts side by side before the
per-`n` table.

The full top-k backfill is recorded in
`FACE_AWARE_FRONTIER_BACKFILL_51_70_2026-04-29.md` and the packet
`EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`. It preserves the old first-hit
`18/20` result, but corrects the interpretation: `n = 70` is a first-hit
tie-selection artifact, while `n = 57` remains the meaningful handoff pinch.

The PMF/collider continuation is now framed in
`PMF_TRANSFER_OPERATOR_ATTACK_2026-04-29.md`. The short version: exact Sidon
maximizers are ground states of a lattice gas, top-k frontier witnesses sample
the low-energy face, `n = 57` is the current defect/handoff candidate, and
`n = 70` is degeneracy/tie-selection rather than a real defect.

The first #755 transfer-state cross-problem gate is recorded in
`B2G_TRANSFER_CROSSPROBLEM_2026-04-29.md`. The `B_2[2]` exact finite scout
uses bounded ordered-sum-count memory and restores persistent field-sensitive
face selection in `28/29` rows over `12 <= n <= 40`. The `B_2[3]` deformation
splits in `18/19` rows over `12 <= n <= 30`, and the first `B_2[2]`
and `B_2[3]` near-ground D2 passes show large smooth `h-1/h-2` layers. That
separates #755 from the rigid #166 sum-free negative control and makes it the
closer Sidon-adjacent mountain for the PMF lane.

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

The Lean side now carries the honest formal shadow of that idea too, but only
at envelope level. The scratch file contains
`sidon_in_range_superfloor_prefix_mass_joint_envelope_external`, which packages
the existing super-floor coarse prefix theorem together with the
density-adjusted mass theorem on the same set. That is the right formal
statement today: one object, two observables, simultaneous envelopes. The
stronger asymmetric Pareto story still lives only in the exact maximizer packet,
not in Lean.
## References

- Balogh, Furedi, Roy (2023). "An upper bound on the size of Sidon sets." Amer. Math. Monthly 130(5). arXiv: 2103.15850
- Lindstrom (1969). "An inequality for B_2-sequences." J. Combinatorial Theory 6(2).
- Erdos, Turan (1941). "On a problem of Sidon in additive number theory." J. London Math. Soc. 16.
- Singer (1938). "A theorem in finite projective geometry." Trans. Amer. Math. Soc. 43.

## License

Apache 2.0

## Author

K. Mendoza, Mendoza Laboratory (ken@mendozalab.io)
