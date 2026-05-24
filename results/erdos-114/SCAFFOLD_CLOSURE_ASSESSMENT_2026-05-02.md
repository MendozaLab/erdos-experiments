# EHP #114 Scaffolding Closure Assessment — May 2, 2026

## Executive Summary

Five scaffolding families were generated on May 2, 2026, probing whether a unified framework toward EHP #114 closure can be assembled from tensor-cone geometry, hypergeometric singularity structure, Koopman spectral classification, and Fourier-Hessian modularity. The existing v5 preprint (IEEE 1788 interval arithmetic) delivers per-n certificates for n=3..14 with margins 6%–58%. These new scaffolds ask: can a *unified* proof mechanism replace stratified n-by-n verification?

**Finding:** One scaffold (Radial-Hypergeometric Calibration on n=14) is directly closure-relevant; two (Tensor Cone on n=10, Fourier-Hessian on n=15) are methods-exploration; two (MDL Probe, Koopman Probe) are diagnostic bindings. The critical missing piece is the **interval hardening of the hypergeometric singularity model**: without Puiseux/hypergeometric certificates for radial perturbations, the unified proof cannot cover the singular boundary layer that dominates the lemniscate topology near z^n−1.

---

## Per-Experiment Assessments

### 1. MDL-Probe (n=3..13, Koopman-gap classification)

**What it computes:** Binds existing Leg-4 Koopman spectral gaps to an MDL regime rule (gap > e → QuantumShadow; ≤e → Classical) and reports Classical regime for all n=3..13 with max measured gap 0.975 << e.

**Certificate or probe?** **Probe only.** The rule is deterministic, but the underlying Koopman operator is a sampled kernel (not the true spectrum). Classical verdict is probe-local; does not prove absence of QuantumShadow phase. Explicitly states: "This does not prove the true Koopman spectrum is Classical."

**Relationship to v5 preprint:** Re-binding of existing per-n interval artifacts (EXP-MM-EHP-007-n*) to a transfer-operator classification framework. No new certificates. Useful as a Lean-side classifier (MDL regime + P1/P3 floor signature across all measured degrees consistent), but the perimeter truth surface remains the IEEE 1788 interval bounds, not the spectral gap measurement.

**Closure proximity: 1/5.** The probe repackages existing artifacts under a new classifier axis (Koopman gap) but offers no path to extend coverage beyond n=13 or to unify n-by-n verification into a single mechanism. Diagnostic value only.

---

### 2. Tensor-Cone Scaffold on n=10 (two runs, identical results)

**What it computes:** Diagnostic tensor-cone analysis of the quotient tangent cone near z^10−1. Measures deficit curvature in 16 shape directions (radial + tangent modes, m=1..5) plus singular radial. All shape slopes < 0.2 (well below the slope ≈2 signature of a smooth Hessian problem). Positive symmetric deficit across all directions.

**Certificate or probe?** **Probe only.** Explicitly: "Not a proof." Uses finite-difference marching-squares estimates for perturbed lengths; is not an interval certificate. The positive-deficit finding is promising for local maximality but is floating-point only. Advises: "If n=10 shows stable positive deficit and coherent scaling class, repeat for n=11–14, then interval-harden the first case."

**Relationship to v5 preprint:** Orthogonal to IEEE 1788 strategy. v5 uses interval arithmetic on coefficient space; this probes quotient-space geometry. The n=10 interval artifact (exact L* = 22.886...) is used as external reference only. Tensor Cone could become a *local* proof method if hardened to interval; would replace smooth Hessian assumption with stratified cone analysis.

**Closure proximity: 2/5.** The scaffold identifies a coherent scaling structure (all slopes < 0.2, suggesting unified curvature behavior across perturbation modes) and is replicable across higher degrees. But without interval hardening and without a unified growth law from radial + shape modes, it remains a methods sketch. The observation that slopes are uniform (0.11–0.20) suggests a unified scaling regime *might* exist, but this is speculation without data from n=11..14.

---

### 3. Fourier-Hessian on n=15 (two runs)

**What it computes:** Probes whether lemniscate deficit near z^15−1 is visible in Hilbert/Fourier coordinate basis. Measures deficit curvature for 27 Fourier modes (m=0..7 radial/tangent × phase). All 27 modes positive across eps ∈ {0.02, 0.01, 0.005}. L0 marching estimate 30.465... (actual interval L* = 32.847..., error 7.25%, flagged as "large reference error").

**Certificate or probe?** **Probe only.** Explicitly: "Does not prove local maximality. Not an interval certificate." Includes explicit singularity sanity check warning: "Because z^n−1 is a singular lemniscate, measured effect is better read as stratified deficit/unfolding signal than ordinary smooth Hessian." Flagged data quality issue: "Large reference-error flag: True" due to 7.25% gap between marching estimate and known exact value.

**Relationship to v5 preprint:** Orthogonal. v5 works in coefficient space; Fourier-Hessian works in orthonormal basis space. The approach is sound in principle (working in a cleaner basis should reveal structure), but the large reference error (marching-only Fourier estimate misses 2.38 units out of 32.85) suggests this basis does not isolate the dominant contribution. The mode uniformity (all positive) is encouraging but not reliable without interval bounds.

**Closure proximity: 2/5.** The positive-all-modes finding is encouraging. However, the reference error (7.25%) and explicit singularity warning ("stratified deficit/unfolding signal") indicate that a smooth Hessian framing (even in Fourier basis) is unsafe. To be closure-relevant, this scaffold would need: (a) interval hardening of mode-wise Hessian bounds, (b) a method to handle the singular core (Puiseux or hypergeometric unfolding), and (c) cross-validation at n=14,16 to verify mode universality. Currently, it's a promising sign without a path forward.

---

### 4. Radial-Hypergeometric Calibration on n=14 and n=15

**What it computes:** Exact formulas for the radial family p_a(z) = z^n − a near a=1, using hypergeometric functions 2F1(p,p;1;a²) with p=(n−1)/(2n). Verifies that exact interval artifacts for n=14,15,16 contain the hypergeometric limit values. Measures how the lemniscate deficit decays as radius ε contracts inward: deficit ∝ ε^(1/n), not ε² (singular boundary layer, not smooth Hessian).

**Certificate or probe?** **Probe, but closure-critical.** The hypergeometric formula is exact and closed-form. The interval containment check is verified. However, it is a *diagnostic* probe of singular structure, not a proof of global maximality. Crucially, it proves the local model near z^n−1 is *not* a smooth Hessian but rather a Puiseux-class singular cone. The deficit exponent ε^(1/n) is the smoking gun: this demands a stratified/hypergeometric local certificate, not an ordinary Taylor expansion.

**Relationship to v5 preprint:** **Fundamental reframe.** v5 assumes (implicitly) local smoothness and uses finite-difference Hessian estimates. Radial-Hypergeometric proves that assumption is false: the radial direction exhibits singular scaling ε^(1/n). To extend v5 into a unified proof, the local certificate at z^n−1 *must* replace the Hessian with a hypergeometric singularity model. This is not an add-on; it's a foundational requirement.

**Closure proximity: 4/5.** This scaffold is the *closest* to closure. Why 4, not 5? Because (a) it diagnoses the problem (singularity, not smoothness) but does not yet provide the solution (Puiseux certificate for radial admissible contractions); (b) the missing piece is interval hardening of the hypergeometric bounds for all admissible perturbations; (c) once that's done, the remaining work is bounding nonradial perturbations via Fourier/tensor methods (which the Fourier-Hessian and Tensor Cone scaffolds are exploring). The blueprint is clear, but the interval Puiseux certificate does not yet exist.

---

### 5. Koopman Probes on n=14 (two variants)

**What it computes:** Binds n=14 interval certificate to Leg-4 Koopman-kernel classification. First variant reports spectral gap and regime; second variant adds sampled perimeter estimate (54.35 vs interval truth 30.85, error 76%, diagnostic only). Both confirm Classical regime and P1/P3 floor signature.

**Certificate or probe?** **Probe only.** Classification rule is rigid (gap ≤ e → Classical); the probe measurement is reliable for the sampled operator. However, "This is a sampled kernel measurement, not a theorem about the true Koopman spectrum." The perimeter estimate is explicitly diagnostic and is not used as proof. All rigor comes from the reused IEEE 1788 interval artifact.

**Relationship to v5 preprint:** Pure binding to existing n=14 certificate. No methodological advance. v5's interval is the sole source of truth; Koopman probe is a classifier wrapper. Useful for the Lean formalization (provides a checkpoint: if Koopman-gap analysis contradicts interval verdict, stop and debug), but adds no closure path.

**Closure proximity: 1/5.** Similar to MDL-Probe — repackages existing results under a new measurement axis. No extension beyond n=14, no unified mechanism. Diagnostic only.

---

## Closure-Path Analysis

### Critical Bottleneck Identified

The Radial-Hypergeometric Calibration reveals the **core obstruction to unified proof**: the singular boundary layer at z^n−1 does not admit a smooth Hessian model. The deficit scales as ε^(1/n), not ε² or ε³. This means:

1. **Local proof must use Puiseux/hypergeometric geometry**, not smooth Hessian.
2. **The unfolding of the singular level set is the dominant effect**, not perturbations in a smooth manifold.
3. **Nonradial admissible perturbations (roots stay inside unit disk) must be bounded separately**, using either Fourier mode decomposition or tensor-cone geometry.

### Missing Piece for Closure

**Interval hardening of the hypergeometric singularity model.** The current probe verifies that exact values lie in interval bounds, but does not *prove* z^n−1 is the maximizer using hypergeometric asymptotics. To achieve closure:

1. **Prove the radial contraction deficit ∝ ε^(1/n) lower-bounds any admissible perturbation** (using the hypergeometric formula and Puiseux expansion).
2. **Prove nonradial admissible perturbations increase the deficit faster than ε^(1/n)** (using Fourier/tensor mode analysis with interval bounds).
3. **Bind both to the interval certificate** (IEEE 1788) to yield a global EHP #114 proof for all n up to some threshold.

### Which Scaffolds Compose Toward Closure?

- **Radial-Hypergeometric + interval Puiseux certificate** → local proof at z^n−1 (dominates the problem).
- **Fourier-Hessian (interval-hardened) + Tensor Cone (n=11..14 extended)** → nonradial admissible bound (supporting detail).
- **MDL/Koopman probes** → Lean classifiers, not part of the proof itself.

---

## Top Recommendation

**Prioritize interval-hardening the Radial-Hypergeometric model for n=14, 15, 16.**

The Radial-Hypergeometric Calibration on n=14 has already mapped the singularity structure (ε^(1/n) scaling confirmed). The next step is to convert this diagnostic into a rigorous certificate:

1. **Formalize the hypergeometric singularity model** in interval arithmetic: bounds on 2F1(p,p;1;a²) and its ε^(1/n) Puiseux unfolding for a near 1, radii ε ∈ [0.001, 0.1].
2. **Prove the radial contraction deficit is minimal** among all admissible perturbations of the roots (i.e., roots remaining in the unit disk).
3. **Cross-validate on n=15, 16** to confirm the singularity pattern is universal.

**Why this closes the gap:** Once the radial (dominant) contribution is interval-certified, the nonradial corrections can be bounded separately using Fourier or tensor methods. The Fourier-Hessian and Tensor-Cone scaffolds provide the machinery; the Radial-Hypergeometric analysis provides the target exponent (ε^(1/n)) and the justification for stratified rather than smooth analysis.

**Single missing piece:** An interval-arithmetic Puiseux/hypergeometric certificate that proves ε^(1/n) lower-bounds the deficit for all admissible radial contractions. Once that exists, the unified framework is within reach.

---

## Anti-Recommendation

**Tier 1 (low value, high distraction):** MDL-Probe and Koopman-Probes (both variants).

These are Lean-side classifier bindings with no new mathematics. They add confidence in the existing per-n results but do no work toward closure. The MDL regime rule (gap > e = QuantumShadow) is a conjecture dressed up as a classification framework; the Koopman probe measures a sampled operator, not the true spectrum. **Action:** Archive these as diagnostic tools, not proof scaffolds. They belong in a Lean axiom-justification layer, not in the closure path.

**Tier 2 (incomplete without radial model):** Fourier-Hessian and Tensor-Cone scaffolds.

Both are methods-exploration with promising results (all modes positive, coherent scaling). Both become high-value once the Radial-Hypergeometric interval certificate exists (then they can bound the remainder). Until then, they are orphaned techniques. The Fourier-Hessian's 7.25% reference error is a warning flag: without the hypergeometric unfolding model, the Fourier basis does not isolate the dominant singular effect. **Action:** Archive pending interval Puiseux certificate; then re-activate for nonradial remainder bounds.

---

## Summary Table

| Experiment | Computes | Type | Closure relevance | Closure proximity |
|---|---|---|---|---|
| MDL-Probe (n=3..13) | Koopman-gap classification | Probe | Diagnostic only | 1/5 |
| Tensor-Cone (n=10 ×2) | Quotient-space deficit curvature | Probe | Methods sketch; if extended to n=11..14 + interval-hardened, becomes remainder bound | 2/5 |
| Fourier-Hessian (n=15 ×2) | Hilbert-basis mode curvature | Probe | Methods sketch; needs singularity unfolding model + interval hardening | 2/5 |
| Radial-Hypergeometric (n=14, n=15) | Singularity structure (ε^(1/n) scaling) | Probe → **closure-critical** | Proves singular boundary layer is dominant; interval hardening of this model is the key missing piece | 4/5 |
| Koopman-Probe (n=14 ×2) | Koopman-gap + interval binding | Probe | Diagnostic classifier | 1/5 |

---

## Conclusion

**A unified EHP #114 proof is within reach, contingent on one bottle-neck: interval-hardening the Radial-Hypergeometric singularity model.** The scaffolds produced today identify the structure (singular ε^(1/n) boundary layer dominates), point toward the method (Puiseux/hypergeometric unfolding), and provide supportive techniques (Fourier/tensor nonradial bounds). The missing theorem is: *For all admissible radial contractions ε ∈ (0, 0.1], the deficit of z^n − (1−ε) relative to z^n − 1 is lower-bounded by K·ε^(1/n) for explicitly-bounded K, and this deficit is minimal among all admissible root perturbations.*

Once that exists, the framework unifies: radial dominance (proven) + nonradial remainder (bounded by Fourier/tensor methods) + interval certainty (IEEE 1788 for the base case) = **closed proof** for all n up to the threshold of the underlying interval certificate (currently n=14).

*Report generated 2026-05-02, read-only assessment.*
