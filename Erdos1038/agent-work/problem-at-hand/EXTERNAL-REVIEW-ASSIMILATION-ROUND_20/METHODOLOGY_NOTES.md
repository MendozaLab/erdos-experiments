# Round 20 — Methodology Notes

Unlike Round 19's assimilation (which surfaced two methodology gaps requiring follow-up), Round 20's methodology is genuinely strong and the notes below are positive-with-followups. PC closed both C1 and C2 with rigorous methods, not just adequate ones.

## §1 — Notably strong methodology choices

### Fraction (exact rational) basis-change construction

PC built the monomial-to-Chebyshev basis-change matrix `P` in Python `Fraction` arithmetic (arbitrary-precision exact rationals), then down-cast to f64 only at the end before SVD. At `g=24, eps=1e-4`, the binomial expansion of `(mid + half·t)^n` for `n=23` involves rationals with numerators and denominators that span hundreds of orders of magnitude. Building in f64 directly produces catastrophic cancellation (PC notes ≥ 100 bits of precision loss, leading to `inf` in the subsequent SVD).

The exact-rational construction sidesteps this entirely. The **raw** cond(P) of 4.6e118 at g=24 is a genuine deterministic feature of the linear map, not f64 noise.

This is the right call and matters specifically for the high-genus regime that #1038's g=24 case sits in. A weaker methodology (f64 throughout) would have produced a misleading `inf` and obscured the actual verdict.

### Iterated Van der Sluis equilibration

The raw cond of P is dominated by the trivial scale gap between the natural norms of the two basis families on a thin gap interval (`(|mid|+half)^(j-1)` for monomials vs. `O(1)` for Chebyshev `T_{j-1}(t)`). This is column-scaling, not linear-map ill-conditioning — a diagonal rescaling makes it disappear without changing what linear maps are representable.

PC computes both raw and equilibrated conds, with 6 iterations of column/row 2-norm equilibration (Van der Sluis 1969). Reporting both is the right move: the raw cond shows the scale gap exists; the equilibrated cond shows what conditioning remains after the trivial gap is removed.

The 6-iteration choice is well within Van der Sluis's convergence regime. Single-pass equilibration is within `√min(m,n)` of optimal; iteration converges fully. At g=24, single-pass would have residual scale imbalance of order `√24 ≈ 5`; six iterations brings the diagonal-scaled cond to its scale-invariant minimum.

### Genuinely non-tautological diagnostic

The structural argument for why `cond(P)` differs from `cond(M)` is correct and explicit: P is constructed from gap geometry alone (`mid = (a+b)/2`, `half = (b-a)/2`), never touches M's cycle integrals or quadrature samples, and is not obtained from M by any orthogonal factorization. This is a verifiable structural fact, not a stipulation.

The empirical confirmation (cond(M) ~ 2.1e11 vs. equilibrated cond(P) ~ 1.18e7 at g=24, four orders of magnitude separation) closes the C1 tautology gap that Round 19's assimilation flagged.

## §2 — Carve-outs PC named (follow-up surface)

PC's `TRANSFORM_DIAGNOSTIC_RATIONALE.md` is unusually explicit about what the verdict does NOT establish. Two carve-outs deserve carry-forward into the next-local-gate:

### Carve-out 1: Basis-change well-conditioned ≠ Chebyshev-rescaled period matrix well-conditioned

The G2 diagnostic establishes that the linear map from monomial coordinates to Chebyshev-rescaled coordinates is well-conditioned under equilibration. It does NOT establish that the period matrix `M_T` (built directly with `T_{j-1}(scaled_x) dx/y` numerators instead of `x^(j-1) dx/y`) has `cond(M_T) < 1e10` at g=24.

The next-local-gate must compute `M_T` directly at g=24 across the same (eps, jitter) sweep and report `cond(M_T)`. Two outcomes are possible:

1. **`cond(M_T)` is comparable to `cond(M)` (~ 2e11) at g=24.** The Chebyshev rescaling is a relabeling, not a fix; the route needs to move further (Option (c) row-space/cycle-space equivalence — local-only).
2. **`cond(M_T)` is dramatically lower** (e.g., `O(10³)` per Trefethen 2019's general principle for switching from monomial to Chebyshev). The defensive pivot is a genuine remediation; G3 interval re-implementation proceeds with `M_T`.

This is a tractable f64 computation, well within PC scope if dispatched in a future round. Alternatively, the local agent can run it directly using PC's existing scaffold.

### Carve-out 2: Conditioning ≠ Certified seed

Both G1 and G2 are `F64_SAMPLED_ONLY` triage probes, not interval-certified seeds. The C3 confound (f64-only sampling) remains open. Even if `cond(M_T)` comes back well-conditioned in carve-out 1, the route still needs G3 interval-arithmetic re-implementation to produce a real receipt.

PC cannot do G3 — it requires the private endpoint payload that PC's sandbox doesn't have access to. G3 stays local-only.

## §3 — One small note on the rationale doc

PC's `TRANSFORM_DIAGNOSTIC_RATIONALE.md` says "Iterated Van der Sluis equilibration ... computes diagonal D_L, D_R that approximately minimize cond(D_L P D_R) over diagonal pairs." This is correct, but technically Van der Sluis 1969 proves the bound for the *one-pass* algorithm; the *iterated* form is a heuristic extension that converges in practice but doesn't have the same closed-form bound. Doesn't affect the diagnostic — PC's iteration count (6) is empirically well above the convergence threshold for this problem — but worth knowing if the methodology is later cited in a more formal context.

## §4 — What's NOT in scope this round (carry forward)

- **G3 — interval-arithmetic re-implementation.** Local-only. Required for any certified seed.
- **G4 — endpoint-limit source vector expression in canonical basis.** Local-only. Required to connect the canonical basis route to the endpoint-limit gate.
- **Six receipts (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS)** — still absent; Round 20 didn't touch the dependent-Vieta consumer path.
- **Coefficient-box theorem, endpoint-limit source kernel, KKT/strict-slack, global reduction** — none touched. Status unchanged.

## §5 — Score-card

| confound | Round 19 status | Round 20 outcome |
|---|---|---|
| C1 (tautological QR diagnostic) | Open | **RESOLVED** (non-tautology proven empirically + structurally) |
| C2 (unprobed g=24) | Open | **RESOLVED** (sweep complete through g=24; FAIL outcome quantified) |
| C3 (f64-only sampling) | Open | Still open (G3 local-only future work) |
| C4 (silent main-fallback) | New (Round 20 first attempt) | RESOLVED at template level (commit `a4b5bc6`); corrective dispatch honored |

Two of the three math confounds Round 19 left open are now closed. C3 is the remaining math gap; closing it is local-only work.

Route working-basis pivots from unrescaled monomial to Chebyshev-rescaled within the canonical hyperelliptic family. No altitude movement. No claim promotion. Six receipts still absent.
