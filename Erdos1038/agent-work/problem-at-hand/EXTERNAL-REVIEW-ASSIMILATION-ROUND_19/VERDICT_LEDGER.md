# Round 19 — Verdict Ledger

External agent: Perplexity Computer
Dispatch: Linear KEN-5 comment `348613ea-3a75-44bd-9a68-84cbf391113a` (2026-05-28T03:31:02Z)
Assimilation status: ASSIMILATED — claim level remains 0 (route management; no theorem advance)

## Bundle

PC returned a substrate-bundle (output rule preference #2; no git format-patch produced). Four files, SHA256-verified, landed cleanly in `Research-Hub/perplexity-substrate/projects/erdos-1038/jobs/ROUND_19/`. All four match Downloads source byte-for-byte.

## Per-artifact verdicts

### CANONICAL_BASIS_SEED.md — ACCEPT_AS_STARTING_PROPOSAL

PC proposed the classical holomorphic differential basis `ω_j = x^(j-1) dx/y, j=1..g` with `g=24` for #1038's 25-component cloud (24 gap rows, normalization excluded). The proposal is well-grounded: cites Frauendiener-Klein 2014 (arXiv:1408.2201) for the basis form and Clenshaw-Curtis period computation, references Mumford's theta-function treatment and Cantor's hyperelliptic Jacobian algorithms as motivating literature with the explicit caveat that they "motivate but do not prove the #1038 route." That literature/motivation separation matches the skill's discipline requirements.

Concrete components are present and acceptable:
- Construction recipe (per-row interval integral over certified gap/cycle)
- Endpoint-safe Chebyshev quadrature (`x = mid + half · cos(θ)`, factor square-root endpoints, interval-enclose smooth product, outward rounding)
- Small-genus specialization (g=2,3,4) is mathematically correct (monomial-power numerators)
- Falsifiable accept/reject thresholds (cond_hi < 1e10, inverse residual norm < 1, smallest singular surrogate > 0, endpoint-limit source expressible in same basis)
- Risk section honestly names the central concern: "The monomial canonical basis may itself be ill-conditioned at genus 24"
- Two defensive alternatives proposed: Chebyshev-rescaled numerators (`T_{j-1}(scaled_x) dx/y`) and row-space/cycle-space equivalence as the stable invariant if condition number fails

Accepted as the **starting proposal** for the canonical-basis primary parallel route. Not promoted to a route certificate — that requires the next-local-gate work described in §NEXT_LOCAL_GATE.

### FALSIFIER_SWEEP_DESIGN.md — ACCEPT_SCOPE_LIMITED

Target A (canonical basis) sweep: genus ∈ {2,3,4}, cluster_eps ∈ {1e-1 down to 1e-4} (7 values), jitter ∈ {0, ±1e-4}. Decision rule: FAIL_TOY at cond(M) ≥ 1e10 or QR-transform diagnostic ≥ 1e10. 63 trials. Well-scoped as a toy-genus stressor.

Target B (dependent-Vieta scaffold envelope): correctly identifies that with the six receipts absent the consumer cannot be tested at theorem level; envelope-only checks (root-count consistency, scaling-convention recorded, etc.) with "any absent or inconsistent field → fail closed." This is the right scaffold posture.

Accepted as scope-limited: the design is sound *for what it tests*. Two scope notes (carried into METHODOLOGY_NOTES.md):
1. Genus range stops at 4. The actual #1038 case is g=24; the toy range does not probe whether the monomial basis breaks at high genus.
2. The "QR-transform diagnostic" needs methodology review (see below).

### NULL_FALSIFIER_REPORT.json — ACCEPT_SCOPE_LIMITED

63 trials, max condition 20.81, max transform_condition 20.81, min singular 1.73, status `NO_FALSIFIER_FOUND_IN_TOY_RANGE`. claim_level 0 and claim_scope explicitly state "toy small-genus stress test only; null result is not proof." Honest scoping.

What this null result actually establishes: at small genus (2-4) and endpoint clustering down to ε=1e-4, the monomial canonical basis remains well-conditioned (cond ~ 20, nowhere near the 1e10 threshold, nowhere near weighted-QR's 4.7e16 failure regime). It does **not** establish anything about g=24.

Accepted as scope-limited supporting evidence for Track B's starting proposal at toy scale.

### round19_canonical_falsifier_sweep.py — ACCEPT_AS_F64_SCAFFOLD

Code is clean, deterministic, numpy-only dependencies. Endpoint quadrature uses the Chebyshev `x = mid + half · cos(θ)` substitution. Branch-point construction with cluster_eps stressor + small alternating jitter is sound. The script honestly declares itself a "public-safe scaffold script" that "does not consume private #1038 endpoint receipts and therefore cannot certify the real 25-component period matrix."

Scope: f64 throughout. Not interval arithmetic. Carries scope tag `F64_SAMPLED_ONLY` per the public repo's tag taxonomy — the lowest tier, not interval-certified.

Reproducibility verified locally — see `REPRODUCIBILITY_CHECK.md`.

## Aggregate verdict

Round 19 is **honestly scoped, well-grounded, and useful as a starting proposal for the canonical-basis primary parallel route**. PC respected every constraint in the dispatch:

- Mode 2 INTENT_ONLY framing throughout
- Claim level 0 on every artifact
- Did not invent any of the six absent receipts
- Did not chain into Round 20
- Did not push or open PRs
- Did not treat ErdosAtlas / Collider / own output as evidence
- Cited public literature with arXiv IDs (verifiable)

The route status from the demote stays put: weighted-QR DEMOTED, canonical hyperelliptic basis PRIMARY PARALLEL, six receipts BLOCKING the dependent-Vieta consumer path. Round 19 gives the canonical basis route a concrete starting proposal but does not move altitude on #1038.

Two methodology flags surfaced (next document).
