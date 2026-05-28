# Round 21 — Methodology Notes

Round 21's methodology is sound and, in one respect, exemplary: PC was honest about a construction caveat before the local agent needed to ask, and built the sweep to be self-validating despite that caveat. The notes below are mostly positive, with one minor flag on the script design.

## §1 — Notably strong methodology choices

### Honest construction caveat with explicit sanity-follow-up recommendation

PC's `chebyshev_m_t_sweep.py` branch_points function is a public-safe scaffold reconstruction from the documented spec — not Round 20's byte-identical substrate code, which lives inside Research-Hub and was not accessible from PC's sandbox. Rather than silently proceeding as if the constructions were identical, PC surfaced this in three places: the MANIFEST.json `construction_caveats` field, the RATIONALE.md "What this does NOT establish" section, and the decision-rule exposition. The recommended sanity follow-up — "when Research-Hub is again sandbox-accessible, replace branch_points() with Round 20's exact construction" — is precisely what the local reproducibility check executed (see REPRODUCIBILITY_CHECK.md).

This is the right behavior. A less careful PC run could have reported a verdict against the scaffold's own anchor and declared the job done. Instead, PC flagged the drift, built the side-by-side M_monomial recomputation to make the ratio construction-invariant, and handed the local agent an explicit sanity path. The qualitative verdict survived that path.

### Side-by-side M_monomial recomputation on the same construction

This is the methodological move that makes the decision rule robust. PC recomputed cond(M_monomial) on the same branch-point construction used for M_T, then applied the decision rule to the ratio cond(M_T)/cond(M_monomial) rather than against Round 20's published cond(M_monomial) anchor (2.12e11) alone.

This matters because the two constructions differ by ~100× at g=24 (Round 20 G1 monomial anchor: 2.12e11; this sweep's same-construction monomial: 1.70e9). If PC had compared cond(M_T) = 2.28e19 against Round 20's monomial anchor of 2.12e11, the ratio would be 1.08e8 — still decisively in the failure regime, but in a way that mixes two different constructions. Comparing within-construction (ratio 1.34e10) is strictly cleaner. The verdict is the same either way, but the methodology is sounder.

### Trefethen-2019 structural explanation

PC's rationale for why M_T is worse than M_monomial (rather than better, as the Trefethen 2019 Vandermonde-with-Arnoldi framing might suggest) is correct and genuinely insightful. The key observation is three-part:

1. Per-gap normalization (`t_i(x)` specific to each gap interval) breaks the global Chebyshev orthogonality structure across rows. Trefethen's result is about evaluation matrices on a single interval where orthogonality holds globally.
2. Integration against `1/y(x)` introduces a smooth per-gap Jacobian. The Chebyshev approximation theorem says smooth-function Chebyshev expansion coefficients decay geometrically — and this is the good news for approximation theory, not here. Geometric column-magnitude decay makes M_T near-rank-1 in each row's Chebyshev expansion, which pushes small singular values toward zero.
3. By contrast, `x^{j-1}` on a gap centered at `c_i` does not decay geometrically in `j` (for `|c_i| ≈ O(1)`, the monomial `|c_i|^{j-1}` grows or stays bounded), so the monomial M keeps non-trivial column magnitudes across all `j`, making it better-conditioned than M_T even though its Vandermonde structure is "theoretically worse" by classical criteria.

This is the failure mode Round 20 METHODOLOGY_NOTES.md §2 Carve-out 1 was implicitly worried about. PC could not preempt it without computing M_T directly — which is exactly what G2.5 did. The structural explanation closes the loop: the experiment result is not a surprise once you understand the mechanism.

## §2 — The question this round answers (in context of Round 20 §2 Carve-out 1)

Round 20 METHODOLOGY_NOTES.md §2 Carve-out 1 read:

> The G2 diagnostic establishes that the linear map from monomial coordinates to
> Chebyshev-rescaled coordinates is well-conditioned under equilibration. It does
> NOT establish that the period matrix M_T ... has cond(M_T) < 1e10 at g=24.

G2.5 is the direct answer to that carve-out. The answer is: no, M_T does not have cond(M_T) < 1e10 at g=24. It has cond(M_T) ~2.28e19 (PC scaffold) or ~3.12e20 (Round 20 branch_points) — eight to nine orders of magnitude above the threshold. The basis-change being well-conditioned (P: equilibrated cond ~1.18e7) does not carry over to M_T's own conditioning, for the structural reasons PC correctly identifies. Round 20's carve-out is now closed with a negative verdict.

## §3 — One minor flag: frozen dataclass import compatibility

`chebyshev_m_t_sweep.py` uses `@dataclass(frozen=True)` on the `SweepResult` container. When the script is invoked as `python3 chebyshev_m_t_sweep.py`, this works without issue. When a downstream consumer tries to import the module programmatically (e.g., `from chebyshev_m_t_sweep import branch_points`), Python 3.9's handling of frozen dataclasses in dynamically-loaded modules can raise reconstruction errors that are not present in Python 3.10+.

This is a minor ergonomics issue, not a methodological flaw — the script is clearly designed to be run as `__main__`, not as an imported library. The recommendation is either drop `frozen=True` (replacing with a regular dataclass or a plain dict) or provide a thin function-only entry point (`def run_sweep() -> dict`) so downstream consumers can call the sweep logic without needing the dataclass type. Neither change affects the numerical output.

The local reproducibility check worked around this by extracting the `branch_points()` function body directly — see REPRODUCIBILITY_CHECK.md for the specific workaround.

## §4 — Confound score-card

| confound | Round 19 status | Round 20 outcome | Round 21 outcome |
|---|---|---|---|
| C1 (tautological QR diagnostic) | Open | **RESOLVED** | Unchanged — RESOLVED |
| C2 (unprobed g=24 conditioning) | Open | **RESOLVED** | Unchanged — RESOLVED |
| C3 (f64-only sampling) | Open | Still open (G3 local-only) | Unchanged — still open |
| C4 (silent main-fallback) | New | RESOLVED at template level | Unchanged — RESOLVED |
| G2.5 (M_T own conditioning) | Named carve-out in Round 20 | Carve-out created — not yet answered | **RESOLVED** with negative verdict |

Three of the four math-relevant tracks (C1, C2, G2.5) are now closed. C3 is the remaining math gap; closing it requires G3 interval re-implementation, which is local-only.

## §5 — What's NOT in scope this round (carry forward)

- **G3 — interval-arithmetic re-implementation.** Local-only. Required for any certified seed regardless of working-basis choice.
- **G4 — endpoint-limit source vector expression in canonical basis.** Local-only. Required to connect the canonical basis route to the endpoint-limit gate.
- **Option (c) row-space/cycle-space equivalence.** Next working-basis candidate, local-only, requires private receipt data.
- **Six receipts** (ROOT_BOX, ROOT_MULTIPLICITY_LEDGER, ORDERED_ROOT_INTERVALS, SCALED_VIETA_IMAGE_CONTRACT, FIXED_CLOUD_BOUND_CERTIFICATE, ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS) — still absent.
- **Coefficient-box theorem, endpoint-limit source kernel, KKT/strict-slack, global reduction** — none touched. Status unchanged.
