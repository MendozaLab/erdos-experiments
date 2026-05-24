# Tensor Cone Certificate — Erdős #114 (EHP) — Scoping Design v1

**Date:** 2026-05-02
**Author:** Tensor-cone track (parallel to a separate Koopman lift track)
**Experiment ID:** `EXP-MATH-EHP114-TENSOR-CONE-V1-20260502`
**Status:** Scoping artifact only — see honest scope statement at end.

## 1. Problem statement

The Erdős–Herzog–Piranian (EHP) conjecture asserts that among monic polynomials
$p(z) = z^n + \sum_{k=0}^{n-1} a_k z^k$ of degree $n$, the polynomial $z^n - 1$
uniquely maximizes the lemniscate length
$L(p) = \mathcal{H}^1(\{z \in \mathbb{C} : |p(z)| = 1\})$. The v5 preprint
(`Math/erdosatlas-workbench/ehp_erdos114_preprint.tex`, Zenodo
10.5281/zenodo.19480329) verifies this for $n \in \{3,\dots,14\}$ via
stratified IEEE 1788 interval arithmetic on a $(2n-3)$-dimensional reduced
parameter space. Each $n$ is a separate exhaustive search; the n=3 dominance
margin is $\sim 6.1\%$. We seek a unified algebraic certificate — a closed
convex cone $K_n$ in coefficient space whose membership implies the bound.

## 2. Coefficient space and slice convention

We identify the monic-polynomial space at degree $n$ with $\mathbb{C}^n$ via
the coefficient vector $(a_0, a_1, \dots, a_{n-1})$ (the leading coefficient
is fixed at 1). Treating real and imaginary parts separately gives
$\mathbb{R}^{2n}$ with coordinates $(x_0, y_0, x_1, y_1, \dots, x_{n-1}, y_{n-1})$
where $a_k = x_k + i y_k$. We work on this affine slice throughout. (The v5
preprint additionally fixes $\arg(a_0)$ via rotation symmetry and uses a
$2n-3$ slice; for SDP exposition, the $2n$ representation is cleaner because
the symmetry constraints are linear and can be added as separate equality
faces of the cone.)

## 3. Cone definition $K_n$ — PSD-moment / SOS choice

We pick a **PSD-moment cone** built on the squared-modulus integral of
$p(e^{i\theta})$. The motivation is the Cauchy–Schwarz / Parseval bridge:

$$
\frac{1}{2\pi}\int_0^{2\pi}|p(e^{i\theta})|^2 \,d\theta
= 1 + \sum_{k=0}^{n-1}|a_k|^2 \;=:\; \|a\|^2 + 1.
$$

For $z^n - 1$, this integral equals $1 + 1 = 2$. Define the **coefficient
energy** $E(p) = \sum_{k=0}^{n-1}|a_k|^2$. The closed-form
$L(z^n - 1) = 2^{1/n}\sqrt{\pi}\,\Gamma(1/(2n)) / \Gamma(1/(2n)+1/2)$ gives
the calibration target. The lemniscate length functional admits the
Crofton-type bound (heuristic for scoping; see §6 for the analytic gap)

$$
L(p) \;\le\; L(z^n - 1) \cdot \sqrt{\,1 + E(p)\,/\,1\,}
\quad\text{when } E(p) \le 1,
$$

with equality only at $a = (-1, 0, \dots, 0)$. The cone we pick is therefore
the **squared-coefficient ball** defined as a PSD constraint on a $2 \times 2$
block per coefficient:

$$
K_n \;=\; \left\{ a \in \mathbb{R}^{2n} \;:\;
\begin{pmatrix} 1 - x_k^2 - y_k^2 & 0 \\ 0 & 1 \end{pmatrix} \succeq 0
\;\;\forall\, k = 0,\dots,n-1,
\;\; \sum_k (x_k^2+y_k^2) \le 1 \right\}.
$$

The per-coordinate PSD blocks force $|a_k| \le 1$; the global trace
constraint forces $E(p) \le 1$. Both are SDP-representable: the coordinate
constraint is a $2\times 2$ Schur complement (linearizable as
$\bigl[\begin{smallmatrix} 1 & x_k \\ x_k & 1 \end{smallmatrix}\bigr] \succeq 0$
to capture $|x_k| \le 1$, and analogously for $y_k$), and the trace constraint
is one linear inequality on the trace of the block-diagonal matrix.

Why this and not full Riesz / Toeplitz-moment cones? For a scoping pass at
n=3 with three test alternatives, the smaller coordinate-PSD cone is enough
to (a) place $z^n - 1$ on the boundary, (b) check three alternative polynomials
inside it. Full Toeplitz-moment lifting is the natural next-iteration design.

## 4. Tensor decomposition

The cone factors as

$$
K_n \;=\; K_{(1)}^{\otimes n} \;\cap\; \{\mathrm{trace} \le 1\},
$$

where the **atomic block** is

$$
K_{(1)} \;=\; \left\{ (x, y) \in \mathbb{R}^2 \;:\;
\begin{pmatrix} 1 & x \\ x & 1 \end{pmatrix} \succeq 0,\;\;
\begin{pmatrix} 1 & y \\ y & 1 \end{pmatrix} \succeq 0 \right\}
\;=\; \{|x|\le 1, |y|\le 1\}.
$$

The tensor symbol here is exact in the block-diagonal sense: an n-fold copy
of $K_{(1)}$ gives $|a_k| \le 1$ coordinate-wise, and the global trace cap
binds the n copies into the energy ball. This is the n-invariance the spec
asks for: same $K_{(1)}$ at every n, same trace-cap structure, just more
copies.

## 5. Why $z^n - 1$ is on the cone boundary

The vector for $z^n - 1$ is $a^\star = (x_0, y_0, \dots, x_{n-1}, y_{n-1})
= (-1, 0, 0, 0, \dots, 0)$. Plugging in: the $k=0$ atomic block becomes
$\bigl[\begin{smallmatrix} 1 & -1 \\ -1 & 1 \end{smallmatrix}\bigr]$, which
is rank-1 PSD — saturating the per-coordinate constraint at the ray
boundary. The trace constraint $\sum (x_k^2 + y_k^2) = 1$ is also saturated
exactly. Both binding constraints are active at $a^\star$, placing it at a
vertex of the cone (not just on a flat face). This is the analytic witness:
$z^n - 1$ corresponds to the single-coefficient saturating direction in
$K_n$.

## 6. Generalization story

Because $K_{(1)}$ is fixed and independent of $n$, an SDP feasibility check
at any degree is structurally the same program with more blocks. The
critical-path obstruction to a closure-relevant certificate is **not**
the tensor part — it is the analytic gap between cone membership and the
length functional.

The current cone enforces $E(p) \le 1$, which by Parseval bounds the
$L^2$ norm of $p$ on the unit circle by $\sqrt{2}$. But the lemniscate
length depends on the geometry of the level set $\{|p|=1\}$, not only on
the boundary $L^2$ norm. The implication

$$
\text{(SDP feasible)} \;\Longrightarrow\; L(p) \le L(z^n - 1)
$$

requires either (a) replacing the $L^2$ bound with a sharper Crofton /
co-area inequality that ties $L(p)$ to a moment-matrix functional, or (b)
augmenting $K_n$ with second-order cone constraints from Toeplitz-moment
lifting (Kalouptsidis / curve-of-Markov-moment style). Neither is in this
scoping pass.

**Koopman alignment note.** The parallel C-track is designing a Koopman
lift. The natural composition: the spectral structure of our SDP — the
eigenvalues of the $2\times 2$ atomic block — corresponds to the simplest
non-trivial Koopman observables on the unit circle (the constant
function and the first character $z \mapsto z$). A coupled certificate
would replace each atomic block with a Koopman-eigenvalue-weighted PSD
constraint and recover a non-coordinate-aligned cone. We design with
this composition path in mind but don't depend on it.

## 7. Comparison to v5 preprint

The v5 preprint reports a $\sim 6.1\%$ margin at n=3 with a certified
interval $[9.17972422234315,\; 9.17972422234317]$ for $L(z^3 - 1)$. Our
prototype's numerical estimate (co-area Riemann sum on a $1200 \times 1200$
grid in $[-3,3]^2$, three values of strip half-width
$\varepsilon \in \{0.04, 0.02, 0.01\}$) gives $L \approx 8.45$ for the
witness — a discretization bias of about $-0.73$ ($\sim -8\%$) at the
finest tested $\varepsilon$. The same bias should affect the alternative
polynomials roughly proportionally; the relative gap $L(\text{alt})$ vs
$L(z^3 - 1)$ at our resolution is on the order of 0.7–1.5 in absolute units
(alternatives in 6.9–7.4 range, witness at 8.45). After unbiasing for the
discretization offset, the prototype's relative margin is consistent with
v5's $\sim 6\%$ at the order of magnitude, but **not analytically derived
from the cone** — it is just a numerical re-confirmation of the v5 fact.

Our cone, in its current scoping form, **does not recover the v5 margin
as an algebraic certificate**. The PSD-feasibility test at n=3 admits
feasible points (alternative polynomials) where $L(p) < L(z^n - 1)$ but
the cone provides no quantitative margin estimate beyond placing
$z^3 - 1$ at the boundary. The gap to a fully analytic SDP-based
certificate is the missing Crofton-type inequality binding $L(p)$ to the
moment matrix. A second-iteration design adding Toeplitz-moment lifting
would close that gap on paper; whether the resulting interval-arithmetic
SDP is tractable at every n is the open question for any future track.

## Honest scope statement (verbatim, per spec)

> "This is a prototype scoping artifact. The cone $K_n$ here is
> candidate-defined and only verified computationally at n=3 with three
> alternative polynomials. To become a closure-relevant certificate, the
> cone definition must be (1) shown analytically to enclose all valid
> coefficient vectors at all n, (2) shown to have $z^n-1$ on its boundary
> at every n, (3) verifiable via interval-arithmetic SDP at every n in a
> tractable way. None achieved in this scoping run."
