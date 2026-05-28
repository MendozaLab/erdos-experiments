# Round 20 — Next Local Gate

After Round 20's two methodology resolutions, the route status is:

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL — working basis pivots from unrescaled monomial to **Chebyshev-rescaled within the canonical hyperelliptic family**. Basis-change is well-conditioned under equilibration at g=24 (cond(P) ~ 1.18e7 << 1e10 threshold, 9 orders below weighted-QR demote cond).
- C3 confound (f64-only sampling) remains open.
- Six receipts still absent.

Recommended next-local-gate packet:

```
EXP-MATH-ERDOS1038-PHI-K-HYPERELLIPTIC-CANONICAL-BASIS-INTERVAL-SEED-20260527-01
```

Same packet ID as the Round 19 assimilation's recommendation. G1 + G2 (PC-shaped) are now CLOSED. G3 + G4 remain (local-only), plus one new follow-up surfaced by Round 20's METHODOLOGY_NOTES.md §2:

## G2.5 — Chebyshev-rescaled period matrix conditioning (CAN BE PC OR LOCAL)

PC's G2 established the **basis change** is well-conditioned but did NOT compute the conditioning of the **Chebyshev-rescaled period matrix** `M_T` directly. This is a tractable f64 computation parallel to G1, can be PC-shaped if dispatched in a future round, or local-shaped if run inline.

**Specification:**

Compute `M_T[i, j-1] = ∫_{gap_i} T_{j-1}(t_i(x)) · y^{-1} dx` where `t_i(x) = (2x - (a_i + b_i)) / (b_i - a_i)` is the Chebyshev rescaling on the i-th gap. Same sweep axes as G1: `g ∈ {2, 4, 8, 12, 16, 20, 24}`, `cluster_eps ∈ {1e-1, ..., 1e-4}`, `jitter ∈ {0, ±1e-4}`. Report `cond(M_T)` per (g, eps, jitter).

**Decision rule:**

- If `cond(M_T) << cond(M)` at g=24 (e.g., `O(10³)` per Trefethen 2019's general principle): defensive pivot is a real remediation. Proceed to G3 with `M_T` as the candidate period matrix.
- If `cond(M_T) ≈ cond(M)` at g=24: Chebyshev rescaling is a relabeling, not a fix. Pivot further to Option (c) row-space/cycle-space equivalence (this becomes local-only because it requires private receipt data to fix alternate gap-row choices).

## G3 — Interval-arithmetic re-implementation (LOCAL-ONLY)

Convert the canonical-basis seed (whichever basis form survives G2.5 — likely Chebyshev-rescaled) from f64 sampling to interval arithmetic. Required additions:

- Outward-rounded interval enclosures on every Chebyshev quadrature node
- Verified upper bound on quadrature truncation error
- Verified condition number bounds (not just point estimates)
- Same endpoint-safe substitution `x = mid + half · cos(θ)`

Output (per the Round 19 recommendation, now updated for Chebyshev-rescaled):

```
CANONICAL_PERIOD_MATRIX_INTERVALS.json
CANONICAL_BASIS_CONDITION_CERTIFICATE.json
CANONICAL_BASIS_SEED_RESULTS.json
```

Scope tag advances: `F64_SAMPLED_ONLY` → `F64_INTERVAL_CONSTRAINED` (or further depending on what's achievable).

## G4 — Endpoint-limit source vector expression (LOCAL-ONLY)

Given the canonical Chebyshev-rescaled basis on the 24-row scaffold, give the explicit construction of how the endpoint-limit source kernel projects onto it. Without this, the route cannot connect to the endpoint-limit gate (which itself remains summit-level open).

Requires private structure (the actual #1038 endpoint payload). Not PC-shaped.

## What stays absent (unchanged from Round 19)

- Six private receipts (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS). Future local-agent deliverables. Independent of the canonical-basis route work.
- No theorem advance is implied by completing G2.5 / G3 / G4. The receipts would promote the canonical-basis route from "starting proposal with extended toy-scale evidence + actionable defensive pivot" to "interval-certified seed in the Chebyshev-rescaled working basis." They do NOT prove #1038 or compose into the global reduction.

## If G2.5 + G3 + G4 complete cleanly

Route status becomes:

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL with **receipt-backed Chebyshev-rescaled seed** (claim level ≤ 2 depending on what scope-tag is reached — likely `FIXED_PROJECTION_DIRECTED_INTERVAL` per the public-repo tag taxonomy)
- `weighted_qr_basis`: stays DEMOTED
- `dependent_vieta_consumer`: still BLOCKED on six absent receipts (unchanged)
- Open summit-level blockers (endpoint-limit kernel, KKT/strict-slack, global reduction): unchanged

## Round 21 dispatch question

Two reasonable Round 21 PC-shaped scopes:

**Option α — G2.5 only.** Tight, focused, builds directly on Round 20 with no re-derivation. PC computes `cond(M_T)` across the same sweep grid. Outcome decides whether the route's working basis is Chebyshev-rescaled or needs to pivot further. Likely 1-2 days of PC compute.

**Option β — Skip G2.5 to PC, do it locally; dispatch PC on a different summit-level question.** For example, PC could be sent on a literature scout for the endpoint-limit source kernel (Task C of the original Round 19 dispatch, deferred). This frees PC from the conditioning-sweep loop and asks for a different kind of work product.

Option α is the more conservative continuation. Option β is the more strategically interesting move if the local agent can run G2.5 itself in a few hours. Either is defensible; the choice belongs to the user.
