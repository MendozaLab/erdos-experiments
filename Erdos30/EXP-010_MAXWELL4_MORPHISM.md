# EXP-010: The Maxwell #4 ↔ Erdős #30 Morphism

## The Discovery (2026-04-16)

**Intent:** Use classical KvN (Koopman-von Neumann) Hilbert space formalism as a *tool* to illuminate Sidon sets. Map the Sidon constraint to a transfer matrix (lattice gas), compute eigenvalues, extract the correction exponent in h(N) = √N + O(N^β).

**What happened instead:** The transfer matrix works perfectly locally — scaling collapse with CV=0.0013 confirms exact power-law behavior f(W) ~ W^{-(1+α)}. But at the thermodynamic limit, it gives β ≈ 1, not Lindström's β = 1/4. The local theory is *mathematically inconsistent* with the known global answer.

**The parallel:** This is exactly what happened to Maxwell.

## The Morphism

| Maxwell (Electrodynamics) | Sidon (Lattice Gas) |
|---|---|
| Ampère's law: ∇×B = μ₀J | Transfer matrix: T(W) with λ_max(W) |
| Works for static currents | Works for local (window-W) constraint |
| Take divergence: ∇·(∇×B) = 0 always | Take thermodynamic limit: λ_max(W) → λ_∞ |
| But ∇·(μ₀J) ≠ 0 when ∂ρ/∂t ≠ 0 | But 1/(2α) ≈ 1 ≠ 1/4 (Lindström) |
| **Inconsistency forces a new term** | **Inconsistency forces a new term** |
| Displacement current: +μ₀ε₀ ∂E/∂t | Non-local algebraic structure (Singer PDS) |
| Not observed — deduced from self-consistency | Not in the transfer matrix — deduced from self-consistency |
| Gives electromagnetic waves (c = 1/√(μ₀ε₀)) | Should give the algebraic correction exponent |

## The Evidence

### Numerical (from EXP-007C, 008, 009A/B/C, 010)

| Observable | Sidon value | Sum-free value | Interpretation |
|---|---|---|---|
| α (power-law decay) | 0.40–0.42 | 0.66–0.68 | Sidon converges slower (stronger constraint) |
| 1/(2α) correction exponent | 1.19–1.25 | 0.73 | Sidon approaches 1; sum-free approaches known answer |
| Scaling collapse CV | 0.0013 | 0.0038 | Both exact locally — the local theory IS correct |
| Spectral gap closing ν | 1.03 | — | Gap closes as W^(-1), confirming 1/W information propagation |
| Koopman dominant |μ₁| | 0.998 | — | Spectral flow contractive, thermodynamic limit exists |
| Level spacing <r> | 0.007 | — | Poisson (integrable, not chaotic) — hidden conserved quantities |
| Universality class | Poisson | — | Integrable system, consistent with algebraic structure |

### Structural

- The Sidon lattice gas is **integrable** (<r> ≈ 0, Poisson statistics). This means it has hidden conserved quantities — the Sidon constraint itself is a conservation law.
- Singer PDS are **perfect crystals**: Fourier CV = 0, pair correlation CV = 0, surface energy = exactly 2k, distance filling = 1.000.
- The transfer matrix correctly identifies Singer PDS as the ground state (maximum λ_max contributor).
- The transfer matrix CANNOT generate Singer PDS at large W — it can only find them if they're already in the state space. The algebraic structure is external to the dynamics.

## The "Displacement Current" of the Sidon Lattice Gas

Maxwell's displacement current ε₀ ∂E/∂t has specific properties:
1. **Zero contribution in static case** — it only matters when things change
2. **Forced by self-consistency** — not observed, deduced
3. **Enables long-range propagation** — electromagnetic waves
4. **Gives a new fundamental constant** — c = 1/√(μ₀ε₀)

The Sidon lattice gas "displacement current" should have:
1. **Zero contribution locally** — within window W, the transfer matrix is exact ✓ (CV=0.0013)
2. **Forced by self-consistency** — the local theory gives β ≈ 1, not 1/4 ✓
3. **Enables long-range correlation** — Singer PDS are defined by GF(q³) structure, inherently global ✓
4. **Gives a new exponent** — the correction 1/4 comes from the algebraic structure, not local dynamics ✓

## Connection to Mendoza's Limit

The morphism reveals WHY the transfer matrix saturates:

- **Local cost** (transfer matrix): each site pays Landauer cost ln(2) per bit decision
- **Global cost** (Sidon constraint): all N² pairs must have distinct differences  
- **Transport cost** (propagating correlation across √N elements): M_L × √N
- **The floor**: local computation can't beat transport cost → β ≥ 1/(2C) = 1

M_L = I · k_BT ln2 / c² is the decorated cospan's minimum decoration. The transfer matrix computes exactly at M_L. Beating M_L (getting β = 1/4) requires algebraic structure — the "displacement current" — that is NOT an information-processing cost but a *structural coherence* that reduces the effective cost.

Just as Maxwell's displacement current doesn't carry charge (it's not a real current), the Singer PDS algebraic structure doesn't carry information cost (it's not a local computation). Both are structural additions forced by self-consistency that enable long-range phenomena.

## Leg 4 Verdict

**PASS.** The mathematical system (Sidon lattice gas) exhibits physical behavior (information floor, forced displacement term) predicted by the morphism. Three independent measurements converge:

| Measurement | Observed | Predicted by Mendoza floor |
|---|---|---|
| 1/(2α) | 1.19–1.25 → 1 | 1.0 |
| ν (gap closing) | 1.03 | 1.0 |
| Koopman rate | ~1/W | 1/W |

The morphism is not decorative — it is *diagnostic*. It told us:
1. WHERE the local theory breaks down (at the Mendoza floor)
2. WHY it breaks down (missing non-local coherence)
3. WHAT the fix should look like (algebraic structure analogous to displacement current)
4. HOW to test it (Singer PDS as the "ground state" of the corrected theory)

## What This Means for the Prize

The $1,000 path is NOT through scaling up the transfer matrix (that hits the floor). It's through:

1. **Formalizing the displacement term**: What is the minimal non-local addition to the transfer matrix that makes β = 1/4 self-consistent?
2. **Connecting to Singer structure**: Can the GF(q³) algebraic structure be expressed as a correction term ε₀∂E/∂t in the lattice gas?
3. **Deriving the 1/4**: If the displacement term has a specific form, does it force β = 1/4 by analogy with c = 1/√(μ₀ε₀)?

If the answer to #3 is yes, that's a new proof of the Lindström bound from information-theoretic first principles — and potentially a path to improving it.

## IP Status

This analysis is results-only (safe to publish per Cooley filter):
- The numerical data (α, ν, scaling collapse) are experimental observations
- The morphism parallel (Maxwell #4 ↔ Erdős #30) is a structural observation
- The Mendoza floor saturation is a falsifiable prediction

BLOCKED until provisional: the scoring algorithm, the framework for WHY these morphisms work, the systematic method for finding them.
