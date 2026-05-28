# Round 19 — Methodology Notes

Two issues surfaced during assimilation that affect interpretation of the null falsifier result. Neither is fatal to Round 19's work product — both are gaps that the next-local-gate or Round 20 should close.

## Flag 1 — The "transform-condition" diagnostic is tautological

The script computes a "weighted-QR-style transform diagnostic" as `cond(R)` where `R` comes from `QR(M)`. This is intended to play the role weighted-QR's coordinate transform played when it produced the 4.7e16 conditioning that triggered the demote.

For real square matrices M, the QR decomposition has Q orthogonal and R upper-triangular with `M = QR`. Since Q is orthogonal it preserves all singular values: `cond(R) = cond(QR) = cond(M)` exactly.

This is confirmed empirically in the report itself: `max_condition = 20.81295453370498` vs. `max_transform_condition = 20.81295453370497`. The two values agree to last-digit roundoff. Independent verification (random 5×5 example): `cond(R)/cond(M) = 1.000000000000`, `||Q^T Q - I|| ≈ 9e-16`.

**Implication:** the "transform stress" diagnostic in this sweep is not actually a separate check — it tests the same quantity as `cond(M)`. PC's null falsifier therefore tests one thing (conditioning of the canonical gap matrix at small genus under endpoint clustering), not two.

This does **not** invalidate the Track B canonical-basis seed proposal. It does narrow the evidence Track A actually provides: the route survives small-genus toy conditioning, full stop. The route's behavior under the kind of weighted-coordinate-transform that destroyed weighted-QR has not been tested by this sweep.

**Next-local-gate must include:** a transform diagnostic that is NOT a tautology — e.g., apply a weighted-coordinate transform analogous to the weighted-QR setup that produced 4.7e16, then compare conditioning before and after. Or compute condition of an explicit basis-change matrix (e.g., from monomials to Chebyshev rescaled, or from canonical to row-space/cycle-space coordinates) and check whether that transform is itself well-conditioned.

## Flag 2 — Genus range stops at 4; the real concern is genus 24

PC's risk section names the central concern correctly: "The monomial canonical basis may itself be ill-conditioned at genus 24." The Vandermonde-like structure of `[x^0, x^1, ..., x^{n-1}]` is famously ill-conditioned at high degree — condition number grows roughly exponentially in n on equispaced nodes, and even on Chebyshev nodes the situation degrades at degree ~20+.

The sweep tests g=2, 3, 4. Nothing about g=24 is established. The toy null result is consistent with both "the canonical basis works at high genus" and "the canonical basis catastrophically fails at high genus" — the sweep simply doesn't probe the regime where the question matters.

This is an honest, named gap (PC flagged it in the risk section), not a failure of Round 19. But it is the question that determines whether the canonical-basis primary parallel route is viable at all.

**Next-local-gate must include:** a direct conditioning probe at g=24, or at least a continuation of the sweep through g=8, 12, 16, 20, 24 to characterize the conditioning growth curve. If the curve crosses 1e10 well below g=24, the route needs the Chebyshev-rescaled or row-space-equivalence variant PC named as defensive alternatives.

## Other notes (minor)

- f64 sampling without interval arithmetic means the entire sweep is `F64_SAMPLED_ONLY` per the public scope tags. Even if Flag 1 and Flag 2 were addressed, the sweep would not produce a route certificate — only an f64 sanity check that informs whether to invest in an interval-arithmetic re-implementation.
- The `cluster_eps` minimum of 1e-4 is reasonable for f64 (well above the ~1e-16 round-off floor). For interval re-implementation, smaller clustering can be probed without precision concerns.
- The sweep uses 160 Legendre quadrature nodes per integral. Quadrature error is not measured. For interval certification, error bounds on the quadrature are required.

## Scope of these flags

These are gaps in what Round 19's evidence supports, not flaws in PC's output. PC was honest about every limitation (claim_level 0, scope statements throughout, risk section naming the genus-24 issue). The flags belong on the **route's evidence record**, not on PC's adherence to the dispatch.

Round 19 advances the canonical-basis route from "named primary parallel" to "starting proposal with toy-scale null check." It does not yet provide a route certificate. The next-local-gate document specifies what it would take to do so.
