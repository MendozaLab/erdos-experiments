# Round 21 — Next Local Gate

The Round 20 G2.5 carve-out is now CLOSED with a negative verdict. After Round 21, the route status is:

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL — working basis is **UNCERTAIN**. Neither unrescaled monomial (FAIL at g=24, Round 20 G1, cond ~2.12e11) nor Chebyshev-rescaled (FAIL at g=24, Round 21 G2.5, cond ~3.12e20 with R20 branch_points) is viable in f64. Round 20's within-family pivot is walked back.
- C3 confound (f64-only sampling) remains open.
- Six receipts still absent.

What has changed since Round 20: the route is no longer on an optimistic forward trajectory within the monomial/Chebyshev family. Two basis candidates have been tested and both fail the Higham 1e10 conditioning threshold at g=24. The next candidate is different in kind, not just a rescaling of the same family.

## Option (c) — Row-space/cycle-space equivalence (LOCAL-ONLY)

This is the natural next pivot given that per-gap basis representations (both monomial and Chebyshev-rescaled) produce ill-conditioned period matrices at genus 24. The row-space/cycle-space approach reformulates the period matrix problem by choosing gap-rows that span the same cycle space as the standard gap basis but with numerically favorable row-column structure — essentially, looking for a basis in which the conditioning arises from global geometry (the arrangement of all 2g+2 branch points) rather than from per-gap normalization artifacts.

This is local-only because the gap-row selection requires knowing which rows produce near-dependent integrals (related to the specific branch-point positions in the private #1038 endpoint payload). PC's sandbox doesn't have access to that data; the computation needs the actual #1038 endpoint structure to be meaningful rather than a scaffolded toy sweep.

No PC dispatch is warranted for Option (c) until the local agent has a candidate row-selection strategy ready for PC to test.

## G3 — Interval-arithmetic re-implementation (LOCAL-ONLY, gate for any certified seed)

G3 remains the gate that converts any working basis from a triage probe into a certified seed. Even if Option (c) produces a well-conditioned basis in f64, G3 is still required before any claim level advances. The scope:

- Outward-rounded interval enclosures on every quadrature node
- Verified upper bound on quadrature truncation error
- Verified condition number bounds (not just point estimates)
- Same endpoint-safe substitution `x = mid + half · cos(θ)`

Output: `CANONICAL_PERIOD_MATRIX_INTERVALS.json`, `CANONICAL_BASIS_CONDITION_CERTIFICATE.json`, `CANONICAL_BASIS_SEED_RESULTS.json`. Scope tag advances: `F64_SAMPLED_ONLY` → `F64_INTERVAL_CONSTRAINED`.

This requires the private endpoint payload (the actual #1038 branch-point data) and is not PC-shaped.

## G4 — Endpoint-limit source vector expression (LOCAL-ONLY)

Given a canonical working basis (whichever survives Option (c) conditioning probe), give the explicit construction of how the endpoint-limit source kernel projects onto it. Without this, the canonical-basis route cannot connect to the endpoint-limit gate. Requires the private #1038 endpoint structure.

## What stays absent (unchanged from Round 19, 20)

Six receipts (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS) remain absent. The receipt bottleneck is the binding altitude constraint. Completing G3 + G4 would advance the canonical-basis route to `F64_INTERVAL_CONSTRAINED` (or better), but the dependent-Vieta consumer remains blocked on the six receipts regardless. The altitude-moving path is the local-Codex receipt sprint, not the canonical-basis route.

## On further PC dispatch for the canonical-basis route

Three PC rounds have now explored the conditioning of the canonical hyperelliptic period matrix at g=24 (Round 19 raised the question; Round 20 G1 quantified the monomial FAIL and G2 tested the basis-change; Round 21 G2.5 tested M_T directly). The sweep-based phase of this investigation is complete:

- Monomial: FAIL at g=24 (cond ~2e11, threshold 1e10, first crossing at g=24)
- Chebyshev-rescaled: FAIL at g=24, worse (cond ~1e20, first crossing at g=8)
- Basis-change between the two: well-conditioned (eq cond ~1.18e7), but this does not help the period matrix itself

The next question — whether Option (c) produces a viable working basis — requires private receipt data that PC's sandbox cannot access. A PC dispatch scoped to "probe row-space/cycle-space equivalence" would either need to work on a synthetic toy (non-informative without the private structure) or would be blocked at the branch-point data step.

The strategically useful PC dispatch at this stage is not on the canonical-basis working-basis question but on other parts of the route where PC's literature access adds value and private data is not required. Rounds 22 and 23 (endpoint-limit literature scout and global-reduction scout, noted in PC's MANIFEST.json as "in flight independently") are the natural continuation.

## If Option (c) + G3 + G4 complete cleanly

Route status becomes:

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL with **receipt-backed row-space/cycle-space seed** (claim level likely ≤ 2 depending on scope-tag reached)
- `weighted_qr_basis`: stays DEMOTED
- `dependent_vieta_consumer`: still BLOCKED on six absent receipts
- Open summit-level blockers (endpoint-limit kernel, KKT/strict-slack, global reduction): unchanged

No theorem advance is implied. The receipts remain the only altitude-moving path.
