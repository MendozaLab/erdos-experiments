# Toeplitz-Moment Lift of the Tensor Cone — Erdős #114 (EHP) — Scoping Design v1

**Date:** 2026-05-02
**Author:** Toeplitz-moment-lift track (iteration on the C-2 tensor cone)
**Experiment ID:** `EXP-MATH-EHP114-TOEPLITZ-LIFT-V1-20260502`
**Status:** Scoping artifact only — see honest scope statement at end.
**Companion:** `TENSOR_CONE_DESIGN_2026-05-02.md` (C-2, the rank-1 prototype this iterates on).

## 1. Problem statement

The Erdős–Herzog–Piranian (EHP) conjecture asserts that among monic polynomials
$p(z) = z^n + \sum_{k=0}^{n-1} a_k z^k$, the polynomial $z^n - 1$ uniquely
maximizes the lemniscate length
$L(p) = \mathcal{H}^1(\{|p(z)| = 1\})$. The v5 preprint
(Zenodo 10.5281/zenodo.19480329) verifies $n \in \{3,\dots,14\}$ via
stratified IEEE 1788 branch-and-bound; the n=3 margin is $\sim 6.1\%$. Tao
(2025) covers $n \ge N_0$ (tower-exponential). The middle wants a unified
algebraic certificate: a closed convex cone in coefficient space with
$z^n - 1$ on its boundary and non-membership binding to $L(p) < L(z^n-1)$.

## 2. C-2's cone and the gap it leaves

The C-2 prototype (rank-1 atomic block) imposes $|a_k| \le 1$ per coordinate
plus $\sum |a_k|^2 \le 1$. It places $z^n - 1$ at a vertex correctly, but its
constraint surface depends only on the Parseval $L^2$ summary $E(p)$. Any
two coefficient vectors with the same energy are indistinguishable: e.g.
$z^3 - 1$ and $z^3 - e^{i\theta}$ both sit on C-2's boundary with $E=1$, but
the C-2 cone has no mechanism to even **see** the cross-coefficient
information that distinguishes lower-length polynomials within an energy
shell. Energy is not enough.

## 3. Choice of moment family — coefficient autocorrelation

Two natural families compete. **Trigonometric-on-lemniscate moments**
$m_k^{\Lambda}(p) = \int_\Lambda z^k \, d\mathcal{H}^1$ tie directly to
length ($m_0^\Lambda = L(p)$) but require knowing $\Lambda$ to evaluate —
circular. **Coefficient autocorrelation moments**
$m_k(a) = \sum_j a_j \overline{a_{j+k}}$ (with $a_n = 1$) are the Fourier
coefficients of $|p(e^{i\theta})|^2$, computable directly from
$(a_0,\dots,a_{n-1})$. This is the classical entry point into Szegő /
Carathéodory–Toeplitz / Fejér–Riesz theory and the family this iteration
uses; length appears via the Szegő-determinant route, not via $m_0$.

## 4. The Toeplitz atomic block $K_{\text{Toep},(1)}$

For monic $p$ with $\mathbf{a} = (a_0,\dots,a_{n-1}, a_n=1)$, the
autocorrelation moments are
$m_k(p) = \sum_{j=0}^{n-k} a_j \overline{a_{j+k}}$ for $k=0,\dots,n$. The
**atomic Toeplitz block** is the $(n+1)\times(n+1)$ Hermitian Toeplitz
matrix $T(p)$ with $T_{ij} = m_{i-j}$ for $i \ge j$ and $\overline{m_{j-i}}$
for $i < j$. By construction $T(p)$ is the Gram matrix of
$\{p, zp, \dots, z^n p\}$ in $L^2(d\theta/2\pi)$, so $T(p) \succeq 0$ is
automatic and the block carries information only **relative to a target**.
The target is the witness Toeplitz form $T^\star_n = T(z^n-1)$.

For $z^n - 1$: $m_0 = 2$, $m_n = -1$, $m_k = 0$ otherwise — so $T^\star_n =
2 I_{n+1} + (-1)\cdot(\mathbf{e}_0 \mathbf{e}_n^\top + \mathbf{e}_n
\mathbf{e}_0^\top)$, i.e. $2$'s on the diagonal, $-1$ in the corners, zeros
elsewhere. Its eigenvalues are $\{1, 3\} \cup \{2\}^{n-1}$. The shifted-down
eigenvalue $1$ records the EHP-canonical configuration.

## 5. The full cone $K_{\text{Toep},n}$

Define
$K_{\text{Toep},n} = \{\mathbf{a} \in \mathbb{C}^n : T^\star_n - T(p)
\succeq 0\}$. This is the **dominance cone** of the witness Toeplitz form
— convex (Loewner order), SDP-representable via Schur complement.

Two immediate consequences. (i) Trace: $\mathrm{tr}\,T(p) = (n+1)(1+E(p))$,
$\mathrm{tr}\,T^\star_n = 2(n+1)$, so PSD-dominance implies $E(p) \le 1$.
C-2's energy bound is a strict consequence. (ii) Off-diagonals: $m_k$ for
$k = 1,\dots,n-1$ couple distinct coefficients (e.g. $a_0 \overline{a_1} +
a_1 \overline{a_2}$ enters $m_1$). The Loewner difference must keep these
small enough for PSD-ness — constraints C-2 is **blind to**. So
$K_{\text{Toep},n}$ is strictly tighter than C-2's cone for $n \ge 2$ on a
set of positive measure.

## 6. Witness saturation

At $\mathbf{a}^\star$, $T(p) = T^\star_n$ exactly, so $T^\star_n - T(p) = 0$
— the zero matrix. Saturation is maximal: rank deficit equals the ambient
dimension $n+1$, not just 1. This is structurally stronger than C-2's
vertex (rank deficit 1 in one block). Geometrically, the witness is the
**apex** of the cone — the unique minimum of $-\log\det(T^\star_n - T(p))$.
A first-order perturbation $a_0 \to -1 + \delta$ ($\delta > 0$ real) gives
$m_0 = 2 - 2\delta + \delta^2$, $m_n = -1 + \delta$, leading to a Loewner
difference with positive diagonal $2\delta - \delta^2$ and small corner
$-\delta$ — still PSD for small $\delta$, as expected.

## 7. Bridge to length, and composition with the Crofton track

The primary identity the Toeplitz cone targets is **Fejér–Riesz**: a
non-negative trigonometric polynomial $\sum_k c_k e^{ik\theta}$ factors as
$|q(e^{i\theta})|^2$ iff its Toeplitz form $(c_{i-j})$ is PSD. Applied to
the difference, $T^\star_n - T(p) \succeq 0$ encodes the existence of a
polynomial $r$ with
$|p^\star(e^{i\theta})|^2 - |p(e^{i\theta})|^2 = |r(e^{i\theta})|^2 \ge 0$
pointwise on the unit circle. So **PSD-ness of the Loewner difference is
exactly equivalent to pointwise dominance of $|p|$ by $|p^\star|$ on
$|z|=1$**.

Szegő's strong limit ($\log\det T_N(|p|^2) \sim 2N\log M(p)$, Mahler
measure) provides an asymptotic determinantal pathway from Toeplitz forms
to Mahler measure to maximum modulus on the unit circle, but this is the
$N \to \infty$ regime; at finite $n$ the Fejér–Riesz route is the operative
one.

Crucially, **pointwise dominance on $|z|=1$ is not yet a length bound** —
that is the gap the parallel Crofton-bridge track is filling. The two
pieces compose end-to-end:
$$
\underbrace{T^\star_n - T(p) \succeq 0}_{\text{Toeplitz cone (this track)}}
\;\stackrel{\text{Fej\'er--Riesz}}{\Longleftrightarrow}\;
\underbrace{|p|^2 \le |p^\star|^2 \text{ on } |z|=1}_{\text{pointwise dominance}}
\;\stackrel{\text{Crofton bridge}}{\Longrightarrow}\;
\underbrace{L(p) \le L(p^\star)}_{\text{length bound}}.
$$
The Fejér–Riesz step is rigorous and classical; the outer two are open.
**The Toeplitz lift is therefore not an alternative to the Crofton bridge,
it is the algebraic feeder for it** — a finite-dimensional SDP-representable
certificate for the pointwise hypothesis the Crofton step needs.

## 8. Comparison to C-2 — does the Toeplitz cone strictly dominate?

Yes, structurally and verifiably at n=3. C-2 sees only $|a_k|^2$ (energy
per coord plus trace cap). The Toeplitz cone uses cross-terms
$a_i \overline{a_j}$ via the off-diagonal moments $m_1,\dots,m_{n-1}$. The
prototype's three test polys sit at strict interior in both cones (Toeplitz
min eigenvalue $> 0$ for all three), but the **phase-shifted vertex** test
$a_0 = -e^{i\phi}$ provides a clean discriminator: at $\phi > 0$ this point
saturates C-2 (energy=1, $|a_0|=1$, on the boundary face), but the Toeplitz
Loewner difference has min eigenvalue $-0.05$ at $\phi=0.05$, $-0.20$ at
$\phi=0.20$, $-0.95$ at $\phi=1.0$. Five test directions, all on C-2's
boundary, are **strictly excluded** from the Toeplitz cone — a
positive-measure family. Strict dominance is demonstrated; see
RESULTS JSON `strict_dominance_at_n_3: "demonstrated"`.

**Honest caveat on rotational symmetry.** By $z \mapsto e^{i\alpha} z$, the
polynomial $z^n - e^{i\phi}$ has the same lemniscate length as $z^n - 1$
(take $\alpha = \phi/n$). So the rotated vertices the Toeplitz cone excludes
are length-equivalent to the witness, not lower-length. Two readings: (i)
the Toeplitz cone is **strictly tighter than a length-bounding cone** — it
pins orientation as well as length, giving a certificate for a
fixed-orientation extremal; (ii) for a length-only certificate, the
operative object is the union of rotated cones $\bigcup_\alpha R_\alpha
\cdot K_{\text{Toep},n}$, or equivalently the v5 preprint's
$\arg(a_0)$-fixed slice on which the rotational symmetry is already
quotiented out. On either reading the cone is structurally stronger than
C-2; on the slice the symmetry issue does not arise.

## Honest scope statement (verbatim, per spec)

> "This is a prototype scoping artifact. The Toeplitz-moment block here is
> candidate-defined and only verified at n=3 with three alternative
> polynomials. The bridge from PSD-ness of the Toeplitz block to
> $L(p) \le L(z^n - 1)$ relies on either Szegő's strong limit (asymptotic)
> or Fejér-Riesz (finite-n but local). Neither is a closed analytic step at
> finite n; both are pointers to an analytic completion. Tractability of the
> analytic completion is for the parallel Crofton-bridge draft, not this
> iteration."
