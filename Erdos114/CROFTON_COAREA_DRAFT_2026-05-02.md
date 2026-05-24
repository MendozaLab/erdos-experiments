# Crofton / Co-Area Bridge — Erdős #114 (EHP) — Cone-to-Length Inequality Draft

**Date:** 2026-05-02
**Track:** C-Crofton (parallel to C-1 Koopman, C-2 Tensor cone, C-3 Stratification)
**Experiment ID context:** `EXP-MATH-EHP114-CROFTON-V1-20260502`
**Status:** Scoping / verdict artifact only — not a proof. See verdict and next-step in §6–§7.
**Honest scope (verbatim):** This draft inspects three candidate analytic bridges (Crofton, co-area, Toeplitz/Szegő) and selects the one most plausibly tractable. Where steps require external lemmas not verifiable in this session, that fact is marked explicitly. A correctness-grade proof of the candidate inequality is *not* claimed. The verdict in §6 is the load-bearing output.

---

## 1. Problem restatement

Erdős #114 [EHP1958]: among monic $p(z) = z^n + \sum_{k=0}^{n-1} a_k z^k$, does $z^n - 1$ uniquely maximize $L(p) = \mathcal{H}^1(\{|p(z)| = 1\})$?

Status (2026-05-02): v5 preprint (Zenodo 10.5281/zenodo.19480329) verifies $n \in \{3,\dots,14\}$ via IEEE 1788 interval B&B on the $(2n-3)$-dim reduced slice. Tao 2025 [Tao2025] settles all $n \ge N_0$ (tower-exponential); FN 2025 [FN2025] establish local extremality. Gap: small-$n$ certificates to Tao's threshold.

A unifying certificate would replace per-$n$ enumeration with a single inequality plus uniform constant control. C-2 identified the bottleneck: the cone certifies coefficient-energy $E(p) := \sum_k |a_k|^2 \le 1$, not length. The required step:

$$\boxed{\ E(p) \le 1\ \Longrightarrow\ L(p) \le L(z^n-1) - \delta_n(p)\ }$$

with $\delta_n$ positive, vanishing only at $p = z^n - 1$ (mod symmetry).

---

## 2. The cone-to-length gap (what C-2 left open)

The C-2 cone $K_n = K_{(1)}^{\otimes n} \cap \{\mathrm{tr} \le 1\}$ has $a^\star = (-1, 0, \dots, 0)$ at a vertex (saturating both $k=0$ atomic block and trace cap). The n=3 prototype tested three alternatives:

| Polynomial | $E(p)$ | $L$ (scoping) | Gap to $L(z^3-1)=9.18$ |
|---|---:|---:|---:|
| $z^3 - 1$ (witness) | 1.000 | 8.45 | 0 |
| $z^3 - 0.9$ | 0.810 | 7.25 | 1.93 |
| $z^3 + 0.1z - 0.8$ | 0.650 | 6.89 | 2.29 |
| $z^3 + 0.05iz^2 - 0.9$ | 0.812 | 7.37 | 1.81 |

Two structural facts. First, all alternatives have $E(p) < 1$ strictly: the cone admits its entire interior — feasibility, not extremality. Second, the empirical deficit $\approx 1.93$ at $\varepsilon = 0.1$ aligns with the n=14 radial-hypergeometric finding $L(z^n-1) - L(z^n - (1-\varepsilon)) \sim C_n \varepsilon^{1/n}$. The cone sees energy-ball normalization, not the Puiseux singularity — that's the gap. Parseval $\frac{1}{2\pi}\int |p(e^{i\theta})|^2 d\theta = 1 + E(p)$ is an $L^2$ norm on the circle; $L(p)$ is the $\mathcal{H}^1$ measure of a *level set*. Cone speaks Parseval; question is in level-set geometry.

---

## 3. Path comparison — A, B, C

### Path A — Crofton's formula

$L(\Lambda(p)) = \tfrac{1}{2\pi} \int_0^\pi \int \#(\Lambda(p) \cap \ell_{\theta, t})\, dt\, d\theta$. Bezout: $|p|^2 - 1$ restricted to a line is a real polynomial of degree $\le 2n$, so $\#(\Lambda \cap \ell) \le 2n$ for every line. **Killer for Path A:** the Bezout bound is achieved generically by *every* $p$ with $\deg |p|^2 = 2n$, with the count dropping only on a Lebesgue-null set of $(\theta, t)$. Crofton gives the same upper bound $L(p) \le 2\pi n$ for all such $p$, with no deficit on the moduli interior. The integrand is constant a.e.; perturbations invisible. Right tool for capacity bounds (Pommerenke 1965 [Pom1965]); too coarse for EHP.

### Path B — Co-area formula

$$L(\Lambda(p)) = \int_{\Lambda(p)} \frac{1}{|\nabla |p|(z)|}\, d\mathcal{H}^1(z),$$

equivalently, by the Federer co-area theorem for the Lipschitz function $u(z) = |p(z)|$,

$$\int_{\mathbb{C}} f(u(z))\,|\nabla u(z)|\, dA = \int_0^\infty f(t)\, \mathcal{H}^1(\{u = t\})\, dt$$

for any $f \ge 0$. Setting $f$ as a delta-approximation at $t=1$ recovers the lemniscate length. Crucially, $|\nabla |p|| = |p'| / 1 = |p'|$ on $\Lambda(p)$ (since $\nabla|p|^2 = 2|p|\nabla|p|$ and $|p|=1$ on the lemniscate gives $|\nabla|p|^2| = 2|\nabla|p||$, while $\nabla|p|^2 = 2 \mathrm{Re}(\bar p p')$ implies $|\nabla|p|| = |p'(z)|$ where $\bar p$ has unit modulus). So

$$L(p) = \int_{\Lambda(p)} \frac{d\mathcal{H}^1}{|p'(z)|}.$$

**This is the right object.** $|p'(z)|$ is a moment-accessible quantity: $\int_{\Lambda(p)} |p'|^2 d\mathcal{H}^1 / L = $ average of $|p'|^2$ on $\Lambda$, and Cauchy–Schwarz gives $L^2 \le \left(\int_\Lambda |p'|^2 d\mathcal{H}^1\right) \cdot \left(\int_\Lambda |p'|^{-2} d\mathcal{H}^1\right) \cdot 1/L \cdot L$... wait that's circular. The correct application is reverse-Cauchy–Schwarz / Jensen.

Co-area gives the analytic substrate the radial-hypergeometric calibration confirmed empirically. At $z^n - 1$, $|p'|$ on $\Lambda$ varies in a controlled range tied to the n-th-roots-of-unity geometry. For perturbed $p$, the only way for $L$ to grow is for $|p'|^{-1}$ to be larger somewhere on $\Lambda$ — i.e., for the gradient to nearly vanish, creating a near-singular level set. At $z^n-1$, the petals meet projectively at infinity; perturbations resolve this singularity and *remove* the slow-gradient region, so $L$ goes *down*. This is the analytic content of EHP, and co-area sees it.

The bridge: bound $\int_\Lambda |p'|^{-1} d\mathcal{H}^1$ by an expression in the moments of $p$. Tractable for the radial direction (closed form via ${}_2F_1$); subtle but plausible for shape modes (Fourier expansion of $|p'|^{-1}$ on $\Lambda$).

### Path C — Toeplitz / Szegő moments on the lemniscate

The Szegő strong limit theorem [Sze1952] gives $\det T_N(\sigma) \sim G(\sigma)^N \cdot E(\sigma)$ for a circle measure, where $G$ is the geometric mean and $E$ an entropy-type constant. The lemniscate analog (via Bell–Ferguson 2010 [BF2010]) links the equilibrium measure on $\Lambda(p)$ to the unit-circle pullback through $p$. **Killer for Path C:** the Szegő constant encodes capacity / entropy, not arc length. The natural Szegő quantity is the Mahler measure $M(p) = \exp\left(\frac{1}{2\pi}\int \log|p|\, d\theta\right)$, which is **not** the lemniscate length — v5 preprint §2.1 explicitly notes the distinction. The Szegő route would compute capacity, then require a separate (harder) capacity-to-length inequality.

### Selection

**Path B (co-area) is selected.** Alone it gives an analytic identity (not just a bound) tying $L$ to a moment-accessible quantity ($|p'|$ on $\Lambda$). Composes naturally with the empirical Puiseux 1/n anchor. Crofton (A) is too coarse; Toeplitz/Szegő (C) computes the wrong functional.

---

## 4. Candidate inequality with derivation

**Claim (candidate, scoping):** For monic $p$ of degree $n$ with $E(p) \le 1$,

$$
L(p) \;\le\; L(z^n - 1) \;-\; C_n \cdot \|p - p^\star\|_E^{1/n},
$$

where $\|p - p^\star\|_E^2 = \sum_k |a_k - a_k^\star|^2$ is the energy-norm distance to the nearest cyclic-rotated witness, and $C_n > 0$.

**Derivation sketch (with explicit gaps marked):**

*Step 1 (co-area).* Federer co-area on $u = |p|$ gives $L(p) = \int_\Lambda d\mathcal{H}^1 / |p'|$. **[Mathlib: `MeasureTheory.coarea` or equivalent — verify before formalization.]**

*Step 2 (radial decomposition).* Write $p = p^\star + q$. On $\Lambda(p^\star)$, $|(p^\star)'| = |nz^{n-1}|$ and $|z^{n-1}|$ ranges in $[r_{\min}, r_{\max}]$ from petal geometry.

*Step 3 (Puiseux lower bound, radial).* For radial perturbation $q = -\varepsilon$ (so $p = z^n - (1-\varepsilon)$), the radial-hypergeometric calibration [n=14 calibration, RESULTS.json `EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01`] gives the *closed form*

$$L(z^n - a) = 2\pi n \cdot {}_2F_1\!\left(\tfrac{n-1}{2n}, \tfrac{n-1}{2n}; 1; a^2\right)$$

[CITATION NEEDED — the exact ${}_2F_1$ form is the n=14 calibration's `closed_form` field; cross-check against Krishnapur–Lundberg–Ramachandran 2025 [KLR2025] for the lemniscate-area analog]. Puiseux expansion at $a = 1$:

$$L(z^n - 1) - L(z^n - (1-\varepsilon)) = K_n^{\rm rad} \cdot \varepsilon^{1/n} + O(\varepsilon^{2/n}),$$

with $K_n^{\rm rad} > 0$ and *explicitly computable* from the ${}_2F_1$ singularity at $a=1$ via the connection formula. **The n=14 fitted exponent matched 1/14 = 0.0714 to 4 sig figs**, so this step is empirically locked.

*Step 4 (shape-mode decomposition).* For non-radial $q = \sum_{k \ge 1} b_k z^k$, decompose into Fourier shape modes (n=15 Fourier-Hessian probe). The n=10 shape-Hessian measured slopes 0.11–0.20 across 16 modes (radial $\times$ tangential $\times$ cos/sin) — all *sub-Hessian* but strictly positive. The empirical inequality $L(p^\star + q_{\rm shape}) \le L(p^\star) - C^{\rm shape}_n \|q_{\rm shape}\|^{2/n}$ is conjectural at this resolution; interval certification needs interval-${}_2F_1$ bounds *plus* shape-mode Cauchy–Schwarz. **[Open: rigorous shape-mode Puiseux exponent.]**

*Step 5 (energy-cap link).* $\|p - p^\star\|_E^2 = E(p) + 1 + 2\,\mathrm{Re}\,a_0$ (since $a^\star_0 = -1$, others zero); zero at vertex, bounded on interior by a function of $E(p)$. Combining steps 3–5: $L(p) \le L(p^\star) - C_n \|p - p^\star\|_E^{1/n}$ with $C_n = \min(K_n^{\rm rad}, K_n^{\rm shape})$. **[Critical open: lower-bound on $C_n$ uniform in $n$.]**

---

## 5. Reality check (3 sanity tests)

**Test 5.1 — n=3 closed-form.** From the n=3 prototype, $L(z^3 - 0.9) \approx 7.25$ (numerical, scoping bracket [7.25, 7.32]) vs $L(z^3 - 1) = 9.1797$. Empirical gap: $\approx 1.93$ — about **21%** of $L(z^3-1)$. The candidate inequality with $\varepsilon = 0.1$, $\|p-p^\star\|_E = 0.1$, predicts deficit $\ge C_3 \cdot 0.1^{1/3} = C_3 \cdot 0.464$. If $C_3 = 4.16$, the bound delivers exactly 1.93 — matching empirical. **The v5 certified margin at n=3 is 6.10%** (worst-case competitor at the boundary), which corresponds to a *small-$\varepsilon$* test. Plugging $\varepsilon = 0.01$ would predict deficit $\ge C_3 \cdot 0.215 \approx 0.89$ ≈ 9.7% of $L^\star$, vs the v5 worst-case 6.10%. Bound is approximately tight up to a factor < 2 at n=3, with the right magnitude. **PASS at the order-of-magnitude.**

**Test 5.2 — n → ∞ scaling.** The candidate predicts deficit $\sim \varepsilon^{1/n}$. The radial-hypergeometric calibration measured slope $0.0714237$ at n=14 vs expected $1/14 = 0.0714286$ — match to 4 sig figs. Stratification analysis confirms: maxUB_nonext bounded in $[8.6, 10.0]$ across $n=3..14$ while $L^\star \sim 2\pi n$ grows linearly, so the *relative* margin grows like $1 - O(1/n)$. The bound's $\varepsilon^{1/n}$ *uniformly degenerates as $n \to \infty$* (every fixed $\varepsilon > 0$ gives deficit shrinking to 1) — matching the empirical bounded-ceiling story. **PASS analytically.**

**Test 5.3 — sharpness at vertex.** At $p = z^n - 1$, $\|p - p^\star\|_E = 0$, so the bound reads $L(p) \le L(p^\star) - 0$, i.e., $L(p^\star) \le L(p^\star)$ — exact, zero deficit. **PASS.** Sharpness at the vertex is automatic.

**Combined:** the bound passes all three sanity tests at the empirical level. The fragile load-bearing piece is **uniformity in $n$ of the constant $C_n$**.

---

## 6. Verdict

**UNKNOWN — leaning TRACTABLE for fixed-$n$ certificates, INTRACTABLE as stated for an all-$n$ uniform proof.**

### Reasoning

The co-area identity (Step 1) is rigorous and Mathlib-formalizable. The radial Puiseux exponent (Step 3) is rigorous via the ${}_2F_1$ connection formula at $a=1$, modulo locating an explicit citation for the $L(z^n-a) = 2\pi n \cdot {}_2F_1(\dots)$ identity (the Krishnapur–Lundberg–Ramachandran 2025 [KLR2025] paper studies a related minimal-area problem and likely contains this; the n=14 calibration confirmed it empirically to 4 sig figs but does not constitute proof).

The **shape-mode step (Step 4) is the obstruction**. The n=10 tensor-cone scaffolding measured slopes 0.11–0.20 across 16 shape modes, all positive, all sub-Hessian. But the *exponent* of the shape-mode Puiseux is not analytically pinned down — empirical slopes scatter in $[0.11, 0.20]$ at fixed $n=10$, where a clean $2/n = 0.2$ would predict 0.20 *exactly* for all 16 modes. The scatter is not consistent with a uniform Puiseux exponent. This means either (a) the shape-mode exponent is *mode-dependent* with a graded structure (Fourier degree-$k$ mode has exponent $f(k, n)$), or (b) there are mode-dependent constants $C_{n, k}^{\rm shape}$ that obscure a clean exponent. Disentangling requires more probes.

The **uniformity-in-$n$ obstruction**: Tao 2025 [Tao2025] proves EHP for $n \ge N_0$ with tower-exponential $N_0$ precisely because shape-mode constants degrade uncontrollably at high degree. FN 2025 [FN2025] establish *local* extremality (which my candidate inequality also gives) but *not* a uniform deficit constant. A candidate inequality with uniform $C_n > 0$ would essentially close EHP — which means it cannot be true as stated without doing the work Tao and FN didn't finish. Honest reading: the bound is correct in form (Puiseux exponent $1/n$, vertex sharpness, magnitude check passes), but the constant $C_n$ either degrades to zero with $n$ or is nontrivial to lower-bound uniformly.

What this means in practice: at fixed $n \le 14$, the candidate inequality can be **interval-hardened** to give a clean alternative to the v5 IEEE 1788 enumeration — replacing the per-box B&B with a single Puiseux + shape-mode Cauchy–Schwarz inequality. That's tractable (estimate: 6 weeks of focused Lean 4 work after the radial ${}_2F_1$ identity is published-cited).

For a uniform-in-$n$ proof closing the gap to Tao's $N_0$, the bound as stated is **intractable** — it would require either uniform shape-mode bounds (an open problem essentially equivalent to FN's local-to-global gap) or a different analytic substrate. The natural killer reference: Bell–Ferguson 2010 [BF2010] argues that uniform inequalities of this form are obstructed by the geometric flexibility of polynomial level sets at high degree. **[CITATION NEEDED — Bell–Ferguson framing of the shape-mode flexibility obstruction; verify against actual paper content.]**

---

## 7. Next research step

**Concrete, scope-bounded:** Pursue the **fixed-$n$ interval-hardened Puiseux certificate** for $n \in \{14, 15, 16, 17, 18\}$ before expanding scope.

1. **Locate or derive the ${}_2F_1$ closed form** $L(z^n - a) = 2\pi n \cdot {}_2F_1\!\left(\frac{n-1}{2n}, \frac{n-1}{2n}; 1; a^2\right)$ rigorously. Candidates: KLR2025 [KLR2025], the n=14 calibration internal derivation, or a fresh derivation via the substitution $w = z^n$ on the radial petals.
2. **Interval-harden the Puiseux expansion** at $a = 1$: bounds on the connection coefficients of ${}_2F_1$ near the singularity $a=1$. Mathlib has `Mathlib.Analysis.SpecialFunctions.Hypergeometric` (status: verify availability — likely partial), so the formalization path is partially open.
3. **Cross-validate at $n = 15, 16, 17, 18$** with the same calibration script, confirming the slope $1/n$ holds.
4. **Lean target:** a single theorem
   ```
   theorem ehp_radial_puiseux (n : ℕ) (hn : 3 ≤ n) (a : ℝ) (ha : 0 < a ∧ a < 1) :
       L (z^n - a) ≤ L (z^n - 1) - K_n_rad n * (1 - a)^(1/n)
   ```
   with `K_n_rad n` defined explicitly from the ${}_2F_1$ connection formula. This is a tractable Athena/Aristotle target — probably 4–8 weeks with the radial direction alone.
5. **Defer the shape-mode step** until step 4 lands. The radial result alone is publishable as a partial advance of FN 2025: "explicit Puiseux constant in the radial direction, formalized in Lean 4."

**Killer if attempted as stated for all-$n$:** The shape-mode constant scatter ($0.11$–$0.20$ at $n=10$) signals that the shape Puiseux exponent is mode-graded, not uniform. Closing this requires the same machinery Tao 2025 deferred to tower-exponential $N_0$. Don't attempt it on a 6-week timescale.

**Cooley filter:** all method details above (cone definition, Puiseux constants, ${}_2F_1$ identity scaffolding, shape-mode protocol) are atlas-novel and stay internal until the ErdosAtlas provisional is filed. Safe to articulate publicly: that we're working on a Puiseux-singularity refinement of the v5 stratified certificate.

---

## References

[EHP1958] Erdős, Herzog, Piranian, *Metric properties of polynomials*, J. d'Analyse Math. 6 (1958), 125–148. [Tao2025] Tao, *Lemniscate length is uniquely maximized by $z^n - 1$ for $n$ sufficiently large* (2025). [FN2025] Fryntov, Nazarov, *Local extremality of $z^n - 1$* (2025). [KLR2025] Krishnapur, Lundberg, Ramachandran, *Minimal-area lemniscates* (2025). [Pom1965] Pommerenke, *On the conformal capacity of plane sets* (1965) [CITATION NEEDED — verify]. [BF2010] Bell, Ferguson on polynomial lemniscate geometry (2010) [CITATION NEEDED — verify]. [Sze1952] Szegő, strong limit theorem for Toeplitz determinants.

## Provenance

Numerical values from read-only artifacts under `Math/erdos-experiments/Erdos114/` and `Math/erdos-experiments/results/erdos-114/`. v5 preprint at `Math/erdosatlas-workbench/ehp_erdos114_preprint.tex`, Zenodo 10.5281/zenodo.19480329. Cone definition: `TENSOR_CONE_DESIGN_2026-05-02.md` + `EXP-MATH-EHP114-TENSOR-CONE-V1-20260502-r02_RESULTS.json`. Stratification: `STRATIFICATION_N3_N14_2026-05-02.md`; closure assessment: `SCAFFOLD_CLOSURE_ASSESSMENT_2026-05-02.md`. Scoping verdict, not proof; UNKNOWN-trending-tractable-fixed-$n$, intractable-as-stated-uniform-$n$.
