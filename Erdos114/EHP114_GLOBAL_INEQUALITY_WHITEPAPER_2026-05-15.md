# Toward a Global Inequality for Erdős #114

**Internal white paper draft**  
**Date:** 2026-05-15  
**Status:** Research program / theorem-target memo. Not a proof. Not public-cleared.  
**Problem:** Erdős-Herzog-Piranian asks whether, among monic degree-\(n\) complex polynomials, the lemniscate length of \(\{|p(z)|=1\}\) is maximized by \(p(z)=z^n-1\).

## Claim Ceiling

This document does **not** claim:

- a proof of Erdős #114;
- a Lean/Mathlib proof of the finite certificates;
- a closed bridge from Tao's sufficiently-large-\(n\) theorem to all intermediate \(n\);
- that ErdosAtlas or the workbench collider has solved the problem.

The safe internal claim is narrower:

> A plausible next route for #114 is to search for a global upper-bound functional that turns lemniscate length into a monotone or stability-controlled quantity of the root configuration. The best current candidate is a log-potential / coarea / root-stability inequality.

This white paper is an internal setup for that search.

## Starting Point

Let

\[
p(z)=\prod_{j=1}^{n}(z-a_j)
\]

be monic, and define the logarithmic potential

\[
u_p(z)=\log |p(z)|=\sum_{j=1}^{n}\log |z-a_j|.
\]

The lemniscate is the level set

\[
\Gamma_p=\{z:u_p(z)=0\}.
\]

For the regular witness

\[
p_\star(z)=z^n-1,
\]

the roots are equally spaced on the unit circle. The target is a global inequality of the form

\[
L(\Gamma_p)\le L(\Gamma_{p_\star})-\mathcal{D}(a_1,\ldots,a_n),
\]

where \(\mathcal{D}\ge 0\), and \(\mathcal{D}=0\) only at the regular \(n\)-gon up to symmetry.

## Why the Previous Routes Are Not Enough

### Crofton alone is too coarse

Crofton's formula expresses length as an average of line intersections:

\[
L(\Gamma)=c\int\#(\Gamma\cap \ell_{\theta,t})\,d\theta\,dt.
\]

For \(\Gamma_p=\{|p|=1\}\), line restrictions are controlled by a real algebraic equation of degree at most \(2n\). This gives a broad degree bound but does not see the difference between \(z^n-1\) and a nearby competitor. It is useful for coarse perimeter/capacity control, not for the EHP maximum.

### Capacity alone is too weak

The compact set \(\{|p|\le 1\}\) has capacity fixed by monicity. A capacity-only inequality would be elegant:

\[
\operatorname{cap}(\{|p|\le1\})\ \text{fixed}\quad\Longrightarrow\quad L(\Gamma_p)\le L(\Gamma_{p_\star}),
\]

but capacity does not generally control perimeter sharply enough. Without polynomial-level-set structure, capacity-perimeter inequalities admit too much geometric flexibility.

### Local collar certificates do not globalize cleanly

The local hard-cell work identifies where interval Newton and wall-separation packets fail. Those artifacts are valuable, but they remain local. The global theorem needs a functional that explains why all such local failures are geometrically forced away from producing excess length.

## Main Proposal: Log-Potential Coarea Stability

The central identity should be based on the coarea formula. On \(\Gamma_p\), where \(|p|=1\),

\[
|\nabla u_p(z)|=\left|\frac{p'(z)}{p(z)}\right|=|p'(z)|.
\]

Length can be read from a delta-localized coarea expression:

\[
L(\Gamma_p)=\int_{\mathbb{C}}\delta(u_p(z))|\nabla u_p(z)|\,dA(z).
\]

Equivalently, in level-set language, length is governed by how much area lies near \(u_p=0\) and how the gradient behaves there. Long lemniscates require extended slow-gradient geometry near the unit-potential level.

This suggests the inequality target:

\[
L(\Gamma_p)
\le
L(\Gamma_{p_\star})
-
c_n\,\Phi(a_1,\ldots,a_n),
\]

where \(\Phi\) is a root-irregularity or potential-instability defect.

## Candidate Defect Functionals

### 1. Root energy defect

Define a repulsion/stability defect against the regular \(n\)-gon:

\[
\Phi_{\rm root}(a)
=
\mathcal{E}(a)-\mathcal{E}(a^\star),
\]

where \(\mathcal{E}\) may be a logarithmic, Riesz, or pairwise angular energy after normalizing translation/rotation/scale.

Desired inequality:

\[
L(\Gamma_p)\le L(\Gamma_{p_\star})-c_n\Phi_{\rm root}(a).
\]

Risk: the correct defect may not be purely root-pairwise. Lemniscate length depends on critical points and level-set geometry, not roots alone.

### 2. Gradient smallness defect

Define

\[
\Phi_{\rm grad}(p)
=
\int_{\Gamma_p} \left(|p'(z)|^{-1}-|p_\star'(z)|^{-1}_{\rm model}\right)_+\,d\mathcal{H}^1(z).
\]

The idea: any competitor that tries to gain length must create slow-gradient corridors. If those corridors can be bounded by critical-point exclusion, then excess length is blocked.

Risk: \(|p'|^{-1}\) is singular near critical points, so the theorem must first exclude or collar critical-point regions.

### 3. Pullback-measure defect

The map \(p:\Gamma_p\to S^1\) is an \(n\)-fold covering away from critical degeneracies. Pull back arc length on \(S^1\). If the witness has the most evenly distributed pullback geometry, define a defect by deviation from uniform pullback density:

\[
\Phi_{\rm pull}(p)
=
\int_{S^1}\left(\rho_p(\theta)-\rho_\star(\theta)\right)^2\,d\theta.
\]

Desired inequality:

\[
L(\Gamma_p)\le L(\Gamma_{p_\star})-c_n\Phi_{\rm pull}(p).
\]

Risk: this may be close to restating the problem unless \(\rho_p\) is bounded by an independent moment or potential estimate.

### 4. Critical-point exclusion defect

Let

\[
F_p(x,y)=|p(x+iy)|^2-1.
\]

The local certificate route wants a lower bound:

\[
|\nabla F_p(x,y)|\ge m>0
\]

on unresolved wall/collar domains. Globalize this by defining a defect that penalizes residual regions where \(|\nabla F|\) can be small:

\[
\Phi_{\rm crit}(p)=\operatorname{Mass}\{z: |F_p(z)|\le \eta,\ |\nabla F_p(z)|\le \tau\}.
\]

Desired theorem shape:

\[
\Phi_{\rm crit}(p)>0\quad\Longrightarrow\quad L(\Gamma_p)<L(\Gamma_{p_\star})-\Delta(\eta,\tau,n).
\]

Risk: this still has a finite-certificate flavor unless the constants are uniform or structurally controlled.

## Collider Translation

The #30 lesson was:

```text
coarse invariant -> gas order parameter -> collider landscape -> stability lemma
```

For #114, the analogous pipeline is:

```text
root configuration
-> log-potential field u_p
-> level-set gas around u_p = 0
-> critical/gradient particles
-> global stability inequality
```

The "gas" is not additive pair sums. It is the distribution of area, gradient, curvature, and pullback density around the lemniscate level set.

The "particles" are not individual atoms of a Sidon set. They are:

- roots \(a_j\);
- critical points of \(p\);
- wall/collar boxes in interval certificates;
- local branches of \(|p|=1\);
- slow-gradient pockets.

The collider should perturb root configurations and measure whether these particles produce monotone changes in:

- length;
- \(\int_{\Gamma_p}|p'|^{-1}\,ds\);
- critical-point proximity;
- pullback-density variance;
- residual wall-separation defect.

## Concrete White-Paper Theorem Targets

### Target A — Coarea length identity

Formal statement:

> For monic \(p\) with regular level \(|p|=1\),
> \[
> L(\Gamma_p)=\int_{\mathbb{C}}\delta(\log|p(z)|)|\nabla\log|p(z)||\,dA.
> \]

This is not enough to solve EHP, but it gives the correct analytic substrate.

### Target B — No slow-gradient excess without critical collars

Formal statement:

> If \(|\nabla F_p|\ge m\) on a collar \(C\) containing \(\Gamma_p\cap D\), then
> \[
> \mathcal{H}^1(\Gamma_p\cap D)\le m^{-1}\cdot \operatorname{AreaBudget}(C,D,p).
> \]

This turns critical-point exclusion into a length bound.

### Target C — Root-irregularity stability

Formal statement:

> After normalizing symmetries, there exists a defect \(\Phi\) such that
> \[
> L(\Gamma_p)\le L(\Gamma_{p_\star})-c_n\Phi(p)+R_n(p),
> \]
> where \(R_n\) is a certified remainder that is negative or dominated by \(c_n\Phi\) on the admissible domain.

This is the real global-inequality target.

### Target D — Fixed-\(n\) interval-hardened version

For \(3\le n\le 14\), or the next frontier \(n=15\), prove a finite version:

> On the reduced admissible coefficient domain, every non-witness region has either a positive defect lower bound or a certified critical-point exclusion collar.

This is the most realistic next step.

## Experiment Plan

### Experiment 1 — Root gas order parameter

**Proposed ID:** `EXP-MATH-EHP114-GLOBAL-ROOT-GAS-20260515-01`

Measure candidate defects over root configurations:

- root-pair energy deviation from the regular \(n\)-gon;
- pullback-density variance on \(\Gamma_p\);
- slow-gradient collar mass;
- length deficit relative to \(z^n-1\).

Acceptance signal:

```text
defect increases whenever length deficit decreases,
with witness as the unique zero-defect configuration.
```

### Experiment 2 — Critical-particle collider

**Proposed ID:** `EXP-MATH-EHP114-CRITICAL-PARTICLE-COLLIDER-20260515-01`

Perturb one root or one coefficient mode at a time. Track:

- critical point migration;
- minimum \(|\nabla F|\) near unresolved walls;
- length response;
- defect response.

Acceptance signal:

```text
every attempted length-increasing perturbation creates a detectable slow-gradient or critical-collar penalty.
```

### Experiment 3 — WS01 remainder landscape

**Proposed ID:** `EXP-MATH-EHP114-N15-WS01-REMAINDER-LANDSCAPE-20260515-01`

Use the existing WS01 residual failures as the testbed. Split the remainder into:

- quadratic contribution;
- cubic contribution;
- mixed residual;
- critical-point proximity term.

Acceptance signal:

```text
failing cells cluster by one dominant defect coordinate,
so a theorem target can route the failure class analytically.
```

## Falsifiers

This global-inequality route should be downgraded if:

1. candidate defects correlate weakly with length deficit;
2. root-regularity defects vanish on non-witness competitors;
3. critical-point exclusion constants degrade too fast in \(n\);
4. fixed-\(n\) interval hardening still requires essentially the same per-box branch-and-bound complexity;
5. the collider only rediscovers the existing n=14 finite certificates without reducing the proof surface.

## Recommended Next Move

Run the root-gas experiment first, then the critical-particle collider.

The root-gas experiment decides whether there is a scalar defect worth making into a theorem. The particle collider then tests whether local perturbations obey that defect. Without the gas stage, collider results risk becoming another finite-search report. With the gas stage, the collider can search for a stability landscape.

## Public-Safe Summary

If this later needs external phrasing:

> I am testing whether lemniscate length can be bounded by a coarea-style stability functional of the polynomial's root configuration. The current goal is not a proof of Erdős #114, but a candidate global inequality that could reduce finite certificate work to a structural estimate.

Do not disclose internal collider scoring, defect formulas, or proof-routing details until the relevant IP/publication gate is cleared.

## Status

White-paper verdict:

```text
PROMISING_GLOBAL_INEQUALITY_PROGRAM
PRIMARY_ROUTE: log-potential + coarea + root/critical stability defect
FIRST_EXPERIMENT: EXP-MATH-EHP114-GLOBAL-ROOT-GAS-20260515-01
CLAIM_CEILING: theorem-target memo only
```
