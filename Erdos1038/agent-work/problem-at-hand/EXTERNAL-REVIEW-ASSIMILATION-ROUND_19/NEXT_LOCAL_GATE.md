# Round 19 — Next Local Gate

Recommended next-local-gate packet ID (matches PC's recommendation):

```
EXP-MATH-ERDOS1038-PHI-K-HYPERELLIPTIC-CANONICAL-BASIS-INTERVAL-SEED-20260527-01
```

This is local work — not a PC dispatch. The receipts and conditioning evidence the gate produces must come from local interval-certified compute against the real 25-component private payload.

## Gate requirements (must all PASS to promote canonical-basis route to receipt-backed level)

### G1. Genus-growth conditioning probe (f64 sufficient at first)

Extend Round 19's sweep through `genus ∈ {2, 4, 8, 12, 16, 20, 24}`. Same endpoint-clustering stressor. Plot `cond(M)` vs genus. Report the genus at which `cond(M)` crosses `1e10`.

- If crossover happens at g ≤ 24: the unrescaled monomial canonical basis fails for #1038. Pivot to PC's defensive alternatives — Chebyshev-rescaled numerators or row-space/cycle-space equivalence.
- If crossover happens at g > 24: proceed to G2 with the unrescaled monomial basis.

This is f64-sufficient because we're triaging, not certifying. The interval re-implementation comes at G3.

### G2. Non-tautological transform diagnostic

Replace `cond(R) from QR(M)` (tautological with `cond(M)`) with a real transform stressor. Options:

- **(a) Weighted-coordinate transform.** Apply the same kind of weighted coordinate transform that produced weighted-QR's 4.7e16 — but using the canonical basis as the underlying matrix. Compute the transform's own condition. If the transform is well-conditioned (≪ 1e10), the route survives the same stressor that demoted weighted-QR.
- **(b) Basis-change matrix.** Compute `cond(P)` where `P` is the change-of-basis matrix from monomials `{x^{j-1}}` to Chebyshev-rescaled `{T_{j-1}(scaled_x)}`. If `cond(P) ≪ 1e10`, the two basis families are numerically interchangeable.
- **(c) Cycle-space equivalence test.** Verify that two different choices of certified gap rows yield matrices related by a well-conditioned similarity — i.e., the row-space invariant PC named is preserved.

Any of (a)/(b)/(c) is a real diagnostic. The current QR-from-M is not.

### G3. Interval-arithmetic re-implementation

Convert the sweep from f64 to interval arithmetic. The existing endpoint-safe Chebyshev quadrature scheme (`x = mid + half · cos(θ)` with square-root endpoint factored) translates directly. Required additions:
- Outward-rounded interval enclosures on every quadrature node
- Bounded quadrature truncation error
- Verified condition number bounds (not just point estimates)

Output the three files PC recommended:

```
CANONICAL_PERIOD_MATRIX_INTERVALS.json
CANONICAL_BASIS_CONDITION_CERTIFICATE.json
CANONICAL_BASIS_SEED_RESULTS.json
```

Once these exist with interval-certified bounds, the scope tag advances from `F64_SAMPLED_ONLY` to `F64_INTERVAL_CONSTRAINED` (or further).

### G4. Endpoint-limit source vector expression

PC's accept criteria include "endpoint-limit source vector expressed in same basis." This is not yet specified. The local gate must produce the explicit construction: given the canonical basis on the 24-row scaffold, how does the endpoint-limit source kernel project onto it? Without this, the route cannot connect to the endpoint-limit gate (Task C from the dispatch — deferred from Round 19 but still a real future blocker).

## What stays absent

- Six private receipts (ROOT_BOX.json etc.) remain absent. Nothing in this local gate changes that. They are still a future-round local-agent deliverable for the dependent-Vieta consumer path, independent of the canonical-basis route work above.
- No theorem advance is implied by completing G1-G4. The receipt would promote the canonical-basis route from "starting proposal" to "interval-certified seed." It does NOT prove #1038 or compose into the global reduction.

## If this gate completes cleanly

Route status becomes:

- **canonical_hyperelliptic_basis**: PRIMARY PARALLEL with receipt-backed seed (claim level ≤ 2 depending on what scope-tag is reached)
- **weighted_qr_basis**: stays DEMOTED to diagnostic
- **dependent_vieta_consumer**: still BLOCKED on six absent receipts (unchanged)
- Open summit-level blockers (endpoint-limit kernel, KKT/strict-slack, global reduction): unchanged

## Round 20 dispatch question

Whether Round 20 goes back to PC depends on whether G1-G4 are PC-shaped or local-shaped:

- G1 is PC-shaped (literature-curated decision on rescaling vs. unrescaled, plus extended sweep). Could dispatch.
- G2 is borderline PC-shaped. The diagnostic-design question can be PC; the implementation is local.
- G3 is local-shaped. Interval arithmetic against private endpoint payload — PC cannot access private receipts.
- G4 is local-shaped (algebraic derivation on private structure).

A reasonable Round 20 dispatch would be: PC produces G1 (conditioning growth curve at genus 8–24 in f64) and G2 (real transform diagnostic design), leaving G3 and G4 for local Codex/Claude work. That's a tight, well-scoped follow-up that builds directly on Round 19 without re-derivation.
