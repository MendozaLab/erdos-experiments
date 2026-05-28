# Round 22 — Next Local Gate

After Round 22's Mode 1 literature scout + schema proposal, the route status is:

- `canonical_hyperelliptic_basis`: PRIMARY PARALLEL — Chebyshev-rescaled working basis (post-Round 20 pivot). G2.5 (cond(M_T) sweep) is in flight via Round 21.
- `endpoint_limit_source_kernel`: SUMMIT-LEVEL OPEN — now has a literature foundation and a schema target for the certificate. Five putative morphisms with falsifiable local tests. Gate still requires local execution.
- C3 confound (f64-only sampling) remains open.
- C5 confound (boundary-component rank-drop risk) is newly open.
- Six receipts still absent.

Four recommended local actions, in priority order:

---

## Action 1 — Run T5 (Eichinger-Lukić / Vieta bound) toy test IMMEDIATELY

This is the cheapest, most immediate test from the Round 22 morphism suite. For f(x) = x^{24}, the cloud is `{x : |x^{24}| < 1}` = `{x : |x| < 1}`, which is a single interval `[-1, 1]`. The equilibrium measure of a single interval is classical; the Robin constant equals log 2 ≈ 0.693 (not log 24 as PC's T5 preamble suggests — this needs verification since the single-interval equilibrium measure is the arcsine measure with Robin constant log(cap([-1,1])) = log(1/2) = −log 2, meaning γ = log 2 for this convention). Run T5's numerical script against the x^{24} case and verify whether the Vieta bound is valid.

**Why immediate:** T5 requires no private data, runs in under a minute, and gives the first concrete data point on the M5 morphism (Eichinger-Lukić Robin constant → FIXED_CLOUD_BOUND). If T5 falsifies M5 on the toy case, that morphism should be dropped before the local agent invests further in it.

**What to check in T5's script:** PC's `robin_constant_from_vieta` function uses `log(2) − (1/n) · Σ log|rⱼ|` for roots near 0. For f(x) = x^{24} all roots are exactly 0, so the formula divides by log|0| which is undefined. The script avoids this by substituting `[0.001] * 24` as near-zero roots. Before trusting the output, confirm that the classical Robin constant for the x^{24} case is what you expect (the equilibrium measure of `{|x^n| < 1}` is the arcsine measure on `[-1,1]` rescaled by `n^{1/n}`, with Robin constant `log(n^{1/n}) = (log n)/n` — this is a different quantity from `log n`). Pin down the correct normalization before comparing the Vieta bound output to a benchmark.

---

## Action 2 — Adopt the schema as the target structure for FIXED_CLOUD_BOUND_CERTIFICATE.json

The schema proposal in `FIXED_CLOUD_BOUND_CERTIFICATE_SCHEMA_PROPOSAL.json` defines the structure that the local agent should populate when private receipt data becomes available. No code change is needed immediately; the adoption is architectural.

Concretely:

1. When the Rust backend eventually produces a FIXED_CLOUD_BOUND_CERTIFICATE.json, that file should validate against the proposed schema.
2. The 10 cross-field rules (VLD-01 through VLD-10) should be embedded as assertions in the backend certificate writer — not just as schema annotations. VLD-04 (six receipts required for FULL_INTERVAL_CERTIFIED) is the most important: the backend should refuse to emit `certificate_status: "FULL_INTERVAL_CERTIFIED"` unless all six receipt references are populated.
3. The `arithmetic_mode` field throughout (`F64_SAMPLED` / `INTERVAL_ARITHMETIC` / `EXACT_RATIONAL`) is the honest-accounting mechanism. The current backend produces F64_SAMPLED results; the schema correctly labels this as non-certifying. When G3 interval arithmetic is added, the mode upgrades to INTERVAL_ARITHMETIC and the certificate can advance to PARTIAL_F64_SAMPLED and eventually FULL_INTERVAL_CERTIFIED.

No schema changes are needed now. If the local agent discovers a needed field (e.g., the T4 b-period data discussed in METHODOLOGY_NOTES §2), propose a schema amendment at that point via a v1.1.0 bump.

---

## Action 3 — Add C5 to the confound tracking list going forward

C5 (boundary-component rank-drop risk at g=24) is a new open confound surfaced by PC's Bogatyrev morphism M3. It is directly relevant to the Round 21 G2.5 result (cond(M_T) sweep, in flight): if the boundary-near component in the actual #1038 cloud has a narrow support interval, cond(M_T) may be high even under Chebyshev rescaling.

When the G2.5 result arrives:
- If cond(M_T) < 1e7 across the (eps, jitter) sweep: C5 does not materialize; the Chebyshev-rescaled pivot is a genuine remediation; G3 proceeds.
- If cond(M_T) spikes for small eps values (small boundary-component widths): C5 is active. The spike pattern will indicate whether it follows M3's predicted `~ 1/w` or `~ log(1/w)` growth (T3 toy test distinguishes these). The remediation would then require either Option (c) row-space/cycle-space equivalence (local-only) or an explicit exclusion of the endpoint-limit source from the g=24 period matrix.

Run T3's toy test (varying boundary-component width from 0.5 down to 0.01, genus g=3) independently of G2.5 to understand the growth regime before the G2.5 result arrives.

---

## Action 4 — When ready to attack the endpoint-limit gate locally, invest in M1 and M2 first

The five morphisms M1–M5 are not equally leveraged for the endpoint-limit gate. When the local agent is ready to run `EXP-MATH-ERDOS1038-PHI-K-ENDPOINT-LIMIT-SOURCE-KERNEL-GATE-20260527-01`, the recommended investment order is:

**M1 (Lubinsky-Bessel) first.** The density-vanishing-exponent fit (T1) is the most information-efficient test. It classifies the boundary component as either Bessel (hard edge, density approaches a positive value or diverges at the endpoint) or Airy (soft edge, density vanishes like square root). This classification determines which local parametrix applies at the boundary branch point, which in turn governs whether the period contribution from that component is finite. The test is runnable in Python in minutes on toy configurations and adaptable to the actual cloud once the ordered root intervals receipt is available.

**M2 (Kuijlaars-Mo Cauchy kernel) second.** The Cauchy kernel residue convergence test (T2) directly probes whether the source kernel for the boundary atom is well-defined as the branch point approaches the physical boundary. If the residue converges to a finite nonzero value as ε→0, the Kuijlaars-Mo construction gives an explicit source-kernel formula for the endpoint atom. If it diverges, the construction fails and the morphism requires revision.

**M3 (Bogatyrev rank-drop)** runs in parallel with G2.5 — it answers the conditioning question and feeds directly into C5 resolution.

**M4 (Deift theta-function)** requires the full symplectic period matrix Ω (see METHODOLOGY_NOTES §2) for a complete test. Defer until the a-period matrix conditioning from G2.5 is known.

**M5 (Eichinger-Lukić Robin)** is addressed by Action 1 (T5 toy test). If T5 holds, M5 provides the cleanest closed-form connection between Vieta coefficients and the cloud bound — potentially the fastest path to a partial FIXED_CLOUD_BOUND certificate without requiring the full period-matrix machinery.

---

## Before / after G2.5 result (from Round 21)

**Before G2.5 (now):** Run Actions 1 and 3. The T5 toy test is immediately executable. T3 toy test can run in parallel. C5 tracked as OPEN.

**After G2.5 arrives:**
- If cond(M_T) < 1e7 at g=24: C5 resolves, C3 becomes the primary remaining math gap, G3 interval re-implementation proceeds with the Chebyshev-rescaled M_T as candidate period matrix.
- If cond(M_T) >= 1e10 at g=24 with narrow boundary component: C5 is active, the route needs Option (c) row-space/cycle-space equivalence before G3 can proceed. This remains local-only.

## What stays absent (unchanged from Round 20)

Six private receipts (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS) — still absent. No theorem advance is implied by completing Actions 1–4. Completing the endpoint-limit gate locally would advance the route toward a receipt, not toward a proof.

## Experiment packet ID (unchanged)

```
EXP-MATH-ERDOS1038-PHI-K-HYPERELLIPTIC-CANONICAL-BASIS-INTERVAL-SEED-20260527-01
```

G3 (interval-arithmetic re-implementation) and G4 (endpoint-limit source vector expression in canonical basis) are the local components of this packet. Both remain local-only.
