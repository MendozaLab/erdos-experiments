# EXP-MM-1038-ENERGY-STABILITY — Experiment Report

**Date:** 2026-04-10  
**Problem:** Erdős #1038 (Extremal Polynomial Sublevel Sets)  
**Certificate:** Upper bound 1.83649 (N=120, atom at +1 w=0.826, cloud in [-1,-0.8])  
**Experiments:** Thomson/Riesz Energy Mapping + Stability Analysis  
**Status:** COMPLETE

---

## Configuration

The N=120 extremal configuration has atom+cloud architecture: a point mass at x=+1 (weight 0.826) and 120 Chebyshev-spaced cloud points in [-1.0, -0.800] (total cloud weight 0.174). For computational tractability, experiments used a reduced representation: atom + 20 cloud points (every 6th Chebyshev node), with weights renormalized to preserve total mass.

---

## Experiment 1: Thomson/Riesz Energy Mapping

**Question:** Does the #1038 extremal configuration match the log-energy minimizer? If so, this connects the sublevel measure optimization (Erdős's problem) to weighted potential theory (Saff-Totik framework), giving a principled explanation for *why* the atom+cloud structure is optimal.

### Results

| Quantity | #1038 Extremal | Log-Energy Minimizer | Δ |
|----------|---------------|---------------------|---|
| Atom weight | 0.826 | 0.809 | 0.017 (2.1%) |
| Cloud support | [-1.000, -0.800] | [-1.000, -0.417] | Right edge differs |
| Cloud center of mass | -0.904 | -0.817 | 0.087 |
| Log-energy | -0.0521 | -0.0566 | 8.7% lower (as expected) |

**Key finding:** The log-energy minimizer reproduces the atom weight within 2.1%. This is striking — two completely different optimization objectives (minimize sublevel measure vs. minimize logarithmic energy) converge to nearly the same atom mass. The cloud support differs: the energy minimizer spreads further right (to -0.417 vs -0.800), which makes physical sense — log-energy penalizes close packing more, so the minimizer spreads mass to reduce interaction energy.

### Riesz Energy Spectrum

| Kernel (s) | Energy |
|-----------|--------|
| 0.25 | 0.1523 |
| 0.50 | 0.1715 |
| 0.75 | 0.2541 |
| 1.00 (Coulomb) | 0.5261 |
| 1.50 | 4.704 |
| 2.00 | 70.46 |
| 3.00 | 31,058 |

The energy grows rapidly with s, confirming the configuration is not optimized for short-range repulsion. The Riesz spectrum shows no anomalous crossings — the extremal measure behaves like a classical potential theory object across all kernel exponents.

### Saff-Totik External Field Analysis

The atom at +1 creates an external field Q(x) = −0.826 · ln|x − 1|. In Saff-Totik's weighted energy framework, the equilibrium measure minimizes I[μ] + 2∫Q dμ, and its support concentrates near the effective potential minimum.

| Quantity | Value |
|----------|-------|
| Effective potential minimum | x = -0.561 |
| Cloud center of mass | x = -0.904 |
| Gap | 0.343 |

The cloud sits *further left* than the effective potential minimum. This is because the cloud isn't minimizing weighted energy — it's minimizing sublevel measure, which has a different cost function. The displacement from the Saff-Totik equilibrium point (0.343 units) quantifies the difference between the two optimization problems.

**Interpretation:** The atom+cloud structure is a *shared architecture* between sublevel minimization and log-energy minimization, but the two objectives tune the parameters differently. The atom weight is nearly universal (2.1% difference), while the cloud geometry is objective-specific. This is analogous to how different Thomson problems on different manifolds share the same local packing motifs but differ in global arrangement.

---

## Experiment 2: Stability Analysis

**Question:** How robust is the 1.83649 bound to perturbations? Is there a single basin, or multiple local optima? What's the curvature (stiffness) of the basin?

### 1D Perturbations

**Atom weight sweep (w ∈ [0.5, 0.95]):**

The sublevel measure |S| shows a clear single minimum. Below w ≈ 0.675 there's a flat plateau at |S| ≈ 2.76 — the atom is too weak to create a deep well, so the sublevel region spans most of [-3, 3]. Above w ≈ 0.675, |S| drops sharply, reaching a minimum at **w = 0.850** (|S| = 1.860) before rising again as the atom dominates and the cloud becomes too light to constrain the polynomial. The #1038 certificate uses w = 0.826, which is within the basin but not at its exact center — the reduced-config optimum shifts slightly.

**Cloud width sweep (δ ∈ [0.05, 0.80]):**

Minimum at **δ = 0.300** (|S| = 1.889). The #1038 certificate uses δ = 0.200, which is near-optimal. Very narrow clouds (δ < 0.1) degrade because they concentrate too much mass at one point, while wide clouds (δ > 0.5) spread mass too thin to constrain the sublevel region.

**Cloud center sweep (center ∈ [-0.95, 0]):**

Optimal at **center = -0.950** (|S| = 1.922), with a strong gradient — moving the cloud rightward monotonically increases |S|. The extremal pushes the cloud as far left as possible (against x = -1), maximizing the distance from the atom at +1. This is the physical intuition: the atom and cloud repel each other to opposite ends of [-1, 1].

### 2D Basin: Atom Weight × Cloud Width

| Parameter | 2D Optimum | #1038 Certificate |
|-----------|-----------|-------------------|
| Atom weight | 0.850 | 0.826 |
| Cloud width δ | 0.150 | 0.200 |
| Sublevel |S| | 1.854 | 1.915 (reduced) |

**Hessian analysis at 2D minimum:**

| Quantity | Value | Interpretation |
|----------|-------|---------------|
| ∂²S/∂w² | 245.76 | Very stiff in atom weight |
| ∂²S/∂δ² | 7.08 | Gentle in cloud width |
| ∂²S/∂w∂δ | 3.78 | Weak coupling |
| Determinant | 1,725.69 | > 0 → stable minimum |
| Trace | 252.84 | > 0 → confirmed minimum |

**The basin is a confirmed stable minimum** (positive-definite Hessian). The condition number (ratio of eigenvalues ≈ 245.76/7.08 ≈ 34.7) reveals highly anisotropic curvature: the system is **35× stiffer in atom weight than cloud width**. This means:

1. Small changes to the atom weight dramatically affect the bound (high sensitivity)
2. The cloud width can vary substantially without much penalty (robust)
3. The atom weight is the "control parameter" — it's the lever that matters most for bound improvement

The weak cross-coupling (∂²S/∂w∂δ = 3.78 vs diagonal terms) means the two parameters are nearly independent — optimizing one barely affects the optimal value of the other.

---

## Synthesis

### Three findings with publication value:

**Finding 1: Atom weight universality.** The log-energy minimizer and sublevel minimizer converge to atom weights within 2.1% of each other (0.809 vs 0.826). This suggests the atom weight ≈ 0.82 is a near-universal feature of the problem geometry, not an artifact of a specific optimization method. The connection to weighted potential theory (Saff-Totik) provides theoretical grounding.

**Finding 2: Single stable basin with anisotropic curvature.** The extremal configuration sits in a single convex basin with 35:1 stiffness ratio (weight vs. width). No local optima detected. The configuration is structurally robust to cloud width perturbations but sensitive to atom weight changes.

**Finding 3: Cloud-potential displacement.** The cloud sits 0.343 units left of the Saff-Totik equilibrium point, quantifying the difference between sublevel minimization and weighted energy minimization. This displacement is a measurable signature of the problem's specific cost function.

### What this does NOT prove:

- This does not improve the 1.83649 bound (the reduced config gives a looser bound)
- This does not prove the conjectured exact value 11/6
- The energy mapping is numerical, not a formal theorem

### Publishable as:

An **addendum** to the existing #1038 erdosproblems.com post or Zenodo deposit — "numerical evidence for atom weight universality and basin structure." This strengthens the certificate by showing the configuration isn't fragile and connects to established potential theory.

---

## Artifact Contract

| File | Status |
|------|--------|
| EXP-MM-1038-ENERGY-STABILITY_RESULTS.json | ✅ Written |
| EXP-MM-1038-ENERGY-STABILITY_RESULTS.sha256 | ✅ Written |
| EXP-MM-1038-ENERGY-STABILITY_REPORT.md | ✅ This file |

**SHA-256:** `333bfbb9023e7bf3c4b2918c987d769a207ea887ad4a7d694c7f7f81681931fc`

### AI Disclosure
Used Claude for experiment design, implementation, and analysis. All computations executed in Python (numpy + scipy).
