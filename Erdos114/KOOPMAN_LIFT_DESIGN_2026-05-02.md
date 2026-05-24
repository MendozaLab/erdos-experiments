# Koopman Observable Lift for Erdős #114 — Design Document

**Experiment ID:** EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502
**Date:** 2026-05-02
**Author scope:** Prototype scoping artifact, not a certified proof.
**IP gate:** ErdosAtlas — NO PROVISIONAL FILED. Cooley filter applies. Do not publish externally.

---

## Honest Scope Statement

> This is a prototype scoping artifact, not a certified proof of EHP at any n.
> The Koopman lift here uses EDMD with finite observable basis; the resulting
> spectral gap is an empirical estimate, not an interval-arithmetic certificate.
> To translate this into a closure-relevant artifact, the spectral gap would
> need to be (1) derived analytically rather than empirically estimated,
> (2) shown n-invariant or growing in n, (3) lifted to interval-arithmetic
> certificates. None of those are achieved in this scoping run.

---

## 1. Problem Statement

Erdős–Herzog–Piranian (1958), problem #114 in the Erdős compendium: among monic
polynomials $p(z) = z^n + \sum_{k=0}^{n-1} a_k z^k$ of degree $n$, does
$z^n - 1$ uniquely maximize the lemniscate length
$L(p) = \mathcal{H}^1\bigl(\{z \in \mathbb{C} : |p(z)| = 1\}\bigr)$? Tao (2025)
proved this for all sufficiently large $n$ via dispersion + Stokes, but his
threshold $N_0$ is tower-exponential and computationally inaccessible. The
v5 preprint (Mendoza, March 2026, Zenodo
[10.5281/zenodo.19480329](https://zenodo.org/records/19480329)) verifies the
conjecture for $n \in \{3, \ldots, 14\}$ via stratified IEEE 1788 interval
arithmetic with branch-and-bound, but does so n-by-n with no unified framework.
The question this design tackles: can a **single Koopman-operator construction**
provide a uniform spectral statement that implies EHP for all $n$ at once?

## 2. Coefficient State Space and Slice Convention

A monic degree-$n$ polynomial has coefficient vector $a \in \mathbb{C}^n \cong
\mathbb{R}^{2n}$. We adopt the v5 preprint's two symmetry reductions:

1. **Translation** $p(z) \mapsto p(z - c)$: kills $a_{n-1}$ (sum of roots = 0).
2. **Rotation** $p(z) \mapsto e^{i\theta} p(e^{-i\theta/n} z)$: fixes the phase
   of $a_0 \in \mathbb{R}_{\ge 0}$.

The reduced slice has real dimension $2n - 3$. For $n = 3$ this is 3 (one
complex $a_1$ and one real $a_0 \ge 0$). The EHP optimum $z^n - 1$ lives at
$x^* = (a_1, a_0) = (0, -1)$ in slice coordinates (we use signed $a_0$ here
because the prototype's L evaluator is sign-symmetric in $a_0$ and centering
on $-1$ avoids boundary effects at $a_0 = 0$).

## 3. Flow Choice — Normalized Gradient Ascent on $L$

Picked: **gradient ascent flow** $\dot x = c \, \nabla L(x)$ with normalization
constant $c = 1 / \mathrm{FLOW\_NORMALIZER}$, FLOW_NORMALIZER $\approx 5 \times 10^6$.

Rationale. Gradient ascent on $L$ makes $z^n - 1$ a *stable* fixed point: at a
strict local maximum, $\nabla L(x^*) = 0$ and the Hessian $H(x^*) \prec 0$, so
the linearized flow $\dot{\delta x} = H \delta x$ has all eigenvalues negative.
This means a Koopman gap statement of the form "the dominant non-trivial
eigenvalue $\lambda^*$ satisfies $\Re \lambda^* \le -\rho_n < 0$" directly
encodes "every initial polynomial $p \neq z^n - 1$ flows toward $z^n - 1$ with
$L$ strictly increasing along the trajectory." Combined with uniqueness of the
fixed point (which the v5 verifier separately establishes), this would give
the EHP statement.

The normalization is necessary because the v5 Level3 results report Hessian
diagonal $[-3.8\text{M}, -5.1\text{M}, -41\text{M}]$ at $n = 3$ — the flow is
extremely stiff. Without a flow-time rescaling, forward Euler would require
$h \le 10^{-7}$ to stay numerically stable, at which point the snapshot
displacement $|\Phi_h(x) - x|$ falls below the marching-squares L noise floor
($\sim 10^{-3}$ relative). The normalization $c = 1/\mathrm{NORMALIZER}$ is a
**flow-time coordinate change**: it rescales all Koopman eigenvalues uniformly
along the real axis, so the qualitative spectral signature (which eigenvalues
sit inside vs outside the unit disk) is preserved.

Alternatives considered and rejected: (a) heat-equation smoothing flow on $L$
itself — interesting but invents a free parameter (diffusion constant) without
a clear spectral interpretation; (b) parametric homotopy from a generic random
$p$ to $z^n - 1$ — would require choosing a path family, breaking the
intrinsic-geometry framing.

## 4. Observable Basis — Monomials of Degree ≤ 2 plus $L$

Picked: $K = 11$ observables for $n = 3$:

$$
\Psi(x) = \bigl[1,\; x_1,\; x_2,\; x_3,\; x_1^2,\; x_2^2,\; x_3^2,\;
                x_1 x_2,\; x_1 x_3,\; x_2 x_3,\; L(x)\bigr]^T
$$

Rationale. Three constraints. (i) $L$ itself must be in the basis — it's the
target observable; without it the EDMD operator can't see the cost surface
structure. (ii) Polynomial monomials of degree ≤ 2 span the linearization
subspace at $x^*$; the Hessian eigenvectors are linear combinations of
$\{x_1, x_2, x_3\}$ and the local L-curvature lives in the
degree-2 subspace. (iii) $K = 11$ is small enough to assemble densely (no
SVD scaling concerns) and large enough that the dominant non-trivial mode
has room to express itself away from the constant eigenfunction.

Alternatives considered: Hermite polynomials in $x$ (no domain advantage over
monomials on the bounded grid), trig basis $\sin/\cos(k \cdot \arg a_1)$ (the
natural invariance is rotation in the *full* coefficient space, but in the
post-symmetry-reduction slice, $L$ has *no* residual rotational structure —
the symmetry has already been quotiented out), radial-only basis (loses
information about which direction the Hessian is shallowest in).

## 5. Koopman Operator Definition + EDMD Approximation

Continuous-time Koopman operator: for an observable $f : \mathbb{R}^{2n-3} \to
\mathbb{R}$ and the flow $\Phi_t$, define $\mathcal{K}_t f(x) := f(\Phi_t(x))$.
For a finite step $h$, $\mathcal{K}_h$ is a bounded linear operator on the
function space (not finite-dimensional in general).

Discrete approximation: **EDMD** (Williams, Kevrekidis, Rowley, 2015). Take a
snapshot dataset $\{x_m\}_{m=1}^M$ in the slice (here, the $7^3 = 343$ regular
lattice in $[x^* - 0.4, x^* + 0.4]^3$). Advance each by one forward Euler step
with the normalized flow to get $y_m = \Phi_h(x_m)$. Assemble the $M \times K$
observable matrices $\Psi_X$ (row $m$: $\Psi(x_m)$) and $\Psi_Y$ (row $m$:
$\Psi(y_m)$). The discrete Koopman matrix is

$$
\widehat{\mathcal{K}}_h = \Psi_Y^T \Psi_X \bigl(\Psi_X^T \Psi_X + \epsilon I\bigr)^{-1}
$$

with a small Tikhonov regularizer $\epsilon = 10^{-10}$ to handle the
near-collinearity between the constant column and the $L$ column near $x^*$.
Eigenvalues of $\widehat{\mathcal{K}}_h$ approximate the discrete-time Koopman
spectrum.

## 6. Predicted Spectral Gap → Certificate Path

Suppose we could establish: **for every $n$, the dominant non-trivial eigenvalue
of the Koopman operator linearized at $z^n - 1$ has real part $\le -\rho_n$
with $\rho_n > 0$, and $\rho_n$ either constant or growing in $n$.** Then EHP
follows in two steps:

1. **Local strict optimality.** A negative-real-part dominant eigenvalue at the
   linearized Koopman around $x^*$ is equivalent to negative-definite Hessian
   of $L$ at $x^*$, since the linearized Koopman generator on the linear
   observable subspace is exactly the Hessian. Negative-definite Hessian gives
   strict local optimality. This is what the v5 preprint already proves
   numerically (their Hessian diagonals are uniformly negative for $n \le 13$).

2. **Global uniqueness.** For the *non-linearized* Koopman operator, a global
   spectral-gap statement says every initial $p \neq z^n - 1$ flows under
   ascent to $z^n - 1$ with $L$ strictly increasing. Combined with the absence
   of secondary critical points (which would manifest as additional
   eigenvalue-1 fixed-point modes in the Koopman spectrum), this rules out
   any other global maximum.

The path from spectral bound → certificate: take the EDMD-estimated bound,
compare to a *certified* Hessian eigenvalue bound (computable via interval
arithmetic, which is what the v5 preprint already does), tighten by including
higher-order observables, and finally lift to interval-arithmetic enclosure
of the *true* Koopman generator's discrete spectrum on a Galerkin truncation.
The last step is the open challenge: known results (Korda-Mezić 2018, Klus
et al. 2020) give EDMD convergence in the $M, K \to \infty$ limit but with
no quantitative rate that's tight enough for a closure-relevant interval bound.

## 7. Generalization Story — What Is $n$-Invariant and What Is Not

| Step | $n$-invariant? | Notes |
|---|---|---|
| Symmetry reduction (translation + rotation) | **Yes** — same recipe at every $n$ | Reduces $\dim = 2n \to 2n-3$ uniformly |
| Flow choice (normalized gradient ascent on $L$) | **Yes** in form | Normalizer scales with $\max\|H_{kk}(z^n-1)\|$ which grows in $n$; the *shape* of the flow is $n$-invariant |
| Stable fixed point at $z^n - 1$ | **Yes** if EHP holds | Tautology in the forward direction |
| Observable basis (monomials degree ≤ 2 + $L$) | **No** in size | Need basis size $\ge 2n-3 + 1$ minimum; for fixed degree-2 truncation, $K = (2n-3)(2n-2)/2 + (2n-3) + 2$, grows quadratically. The *form* generalizes |
| EDMD assembly procedure | **Yes** — pure linear algebra | Cost is $O(M K^2)$ where $K = O(n^2)$ for degree-2 truncation |
| Snapshot grid | **No** in size | Lattice $r^{2n-3}$ blows up; need either Monte Carlo sampling or sparse grids for $n \ge 5$ |
| Spectral gap analytic bound | **Open** | This is the bottleneck. We have no analytic lower bound on $\rho_n$ — neither from Tao's proof nor from the v5 preprint's empirical Hessian eigenvalues. **Critical-path step.** |
| Interval-arithmetic certification | **In principle yes**, in practice the bottleneck of the bottleneck | Once analytic $\rho_n$ exists, certifying it via inari-style intervals is mechanical but expensive. The cost again scales with the basis size |

The bottleneck is step 7: **converting the empirical EDMD spectral gap into an
analytic lower bound on the true Koopman gap, uniformly in $n$.**
Everything before that step is a finite refinement that scales polynomially with
$n$ (via dimension and basis-size growth) but is structurally identical at every
$n$. Everything after that step is mechanical interval-arithmetic plumbing.

## 8. Comparison to v5 Preprint at $n = 3$

The v5 preprint reports certified margin $\sim 6.1\%$ at $n=3$ — i.e., the best
non-extremal polynomial in the search space achieves $L < 0.939 \cdot L(z^3-1)$.
Our prototype's coarse-grid sample margin (the worst $L$ on the $343$-point
grid divided by $L^*$) is **4.53%** — within 2pp of the certified value, and a
positive sign that the EDMD grid actually probes the structurally relevant
neighborhood of $x^*$. The closed-form $L(z^3-1) = 9.179724\ldots$ is recovered
by our marching-squares evaluator to relative error $9.1 \times 10^{-4}$.

The Koopman spectrum, however, **does not currently resolve a positive gap**.
Across a step-size sweep $h \in \{0.1, 0.3, 1.0, 3.0, 10.0, 30.0\}$ (in
flow-time units after the FLOW_NORMALIZER rescale), the dominant non-trivial
$|\lambda_d|$ stays within $10^{-4}$ of unity, never crossing into a clearly
*decaying* regime. The interpretation: with a marching-squares L evaluator at
$\sim 10^{-3}$ relative precision and a finite degree-2 monomial basis, the
spectral gap signal is below the EDMD noise floor. The constant-function
eigenvalue *does* track 1 to ~$10^{-5}$, confirming the basis includes the
correct trivial mode; but the linear-and-quadratic monomial modes that should
encode the Hessian-driven contraction toward $x^*$ are too tightly clustered
near $\lambda = 1$ to be resolved.

To get a meaningful spectral gap from this same construction would require:
(a) replacing marching-squares with a higher-precision L evaluator (e.g.
mpmath at 30+ digits, or the inari-bound chain itself), (b) expanding the
basis to degree 3 or 4 (raises $K$ to ~30 for $n=3$), or (c) using a finer
snapshot grid combined with mean-centering of observables to suppress the
trivial eigenmode. Item (a) is the binding constraint and where the next
iteration should focus.

---

## Files Produced This Session

- `KOOPMAN_LIFT_DESIGN_2026-05-02.md` — this document
- `koopman_lift_prototype.py` — runnable Python prototype (≈400 lines)
- `EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_RESULTS.json` — eigenvalues, sweep, margin
- `EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_REPORT.md` — auto-generated terse summary
- `EXP-MATH-EHP114-KOOPMAN-LIFT-V1-20260502_RESULTS.sha256` — reproducibility hash

## Inputs Read

- `Math/erdosatlas-workbench/ehp_erdos114_preprint.tex` (v5 preprint, lines 34–245)
- `Math/erdos-experiments/results/erdos-114/EHP_N3_LEVEL3_RESULTS.json`
  (Hessian, margin, closest competitor)

## Inputs NOT Modified (per task constraints)

- The v5 preprint TeX, all `EHP_N*_*.json` files in `results/erdos-114/`,
  and the existing N14-KOOPMAN-PROBE artifacts. This V1 is a fresh experiment
  ID with no overwrite collision.
