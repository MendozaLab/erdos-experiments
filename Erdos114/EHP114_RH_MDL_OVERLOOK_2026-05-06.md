# EHP114 -> RH MDL Overlook

Date: 2026-05-06
Scope: internal research-program note
Status: OVERLOOK / PROGRAM TARGET, not theorem progress

## Meaning

EHP114 does not hand over a route to a proof of the Riemann Hypothesis. The
honest bridge is narrower and more useful: EHP114 is a calibrated example where
a short, highly symmetric description (`z^n - 1`) sits at the extremal boundary
of a complex-analytic length functional. RH has several equivalent formulations
that already look like approximation, positivity, or extremal-geometry
statements. The research program is to make one of those formulations carry an
explicit description-length curve.

The clean version is:

> EHP114 is a proved calibration point for MDL extremality. Beurling-Nyman is
> the RH-shaped target where an MDL curve can be made mathematically precise.

This is not "solve RH with EHP." It is "use EHP as a ground-truth theorem in a
class of MDL extremality statements, then ask whether a Beurling-Nyman MDL rate
recovers known RH-equivalent or zero-free-region criteria."

## Claim Ceiling

Safe:

- "EHP114 suggests an MDL-extremality pattern worth testing against RH
  equivalents."
- "Beurling-Nyman is the strongest bridge because it is already an
  approximation criterion equivalent to RH."
- "Arithmetic zeta truncation level sets are lemniscate-like diagnostics, not
  literal EHP polynomial lemniscates."
- "A useful near-term result would be an implication such as: RH or a
  zero-free-region criterion forces a specific asymptotic bound on a
  Beurling-Nyman MDL curve."

Unsafe:

- Any claim that EHP114 implies RH.
- Any claim that RH is equivalent to an EHP arclength bound before a theorem is
  written.
- Any claim that truncated zeta level sets are polynomial lemniscates in the
  strict one-variable EHP sense.
- Any public claim exposing atlas scoring, transfer-operator construction, or
  unpublished MDL internals before the ErdosAtlas IP filter is satisfied.

## External Anchors Checked

- Tao, "The maximal length of the Erdos--Herzog--Piranian lemniscate length in
  high degree", arXiv:2512.12455. Tao proves the EHP conjecture for all
  sufficiently large `n`, building on Fryntov-Nazarov.
- Fryntov and Nazarov, "New estimates for the length of the
  Erdos-Herzog-Piranian lemniscate", arXiv:0808.0717. They prove local
  maximality at `z^n - 1` and the asymptotically sharp `2n + o(n)` bound.
- Alouges, Darses, and Hillion, "Polynomial approximations in a generalized
  Nyman-Beurling criterion", Journal de theorie des nombres de Bordeaux 34
  (2022), 767-785. Their abstract explicitly frames Nyman-Beurling as an
  approximation problem equivalent to RH, with coefficient-control and
  Gram/Hankel structures.

## Bridge Ranking

### Bridge A - Beurling-Nyman MDL Curve

This is the primary bridge.

The Nyman-Beurling formulation asks whether the indicator target can be
approximated in an `L^2` space by a span of normalized fractional-part
functions. That is already an MDL-shaped object: choose a dictionary, choose a
coefficient budget, and measure how much approximation gap closes per
description bit.

Define a predeclared dictionary `B_N = {rho_theta : theta in Theta_N}` and a
cost model `C(a, theta)` that charges for:

- number of active basis functions,
- precision of the chosen `theta` values,
- coefficient precision,
- and any structural restriction such as geometric or randomized `theta`
  schedules.

Then define the MDL approximation curve

```text
K(epsilon; B_N, C)
  = min C(a, theta)
    subject to ||1_(0,1] - sum_k a_k rho_{theta_k}||_2 <= epsilon.
```

The research target is not to prove RH immediately. The first publishable target
is weaker:

```text
RH or a standard zero-free-region condition
  => an explicit asymptotic upper bound on K(epsilon)
```

or conversely,

```text
an explicit asymptotic bound on K(epsilon)
  => a known zero-free region for zeta.
```

That would turn "MDL extremality" into a theorem-shaped invariant rather than a
story.

### Bridge B - Arithmetic Level Sets Of Zeta Truncations

This is secondary and must be worded carefully.

Finite Euler products and Dirichlet polynomials produce unit-level sets such as
`|zeta_N(s)| = 1` or `|sum_{n <= N} n^{-s}| = 1`. These are useful
lemniscate-like diagnostics on the critical strip, but they are not literally
the same as EHP polynomial lemniscates in one complex variable. The corrected
question is:

```text
Among arithmetic truncations with the same description budget, does the
critical-line unit-level geometry extremize a length, variation, or curvature
functional relative to nearby off-line probes?
```

The value of this bridge is experimental. It can produce numerical evidence and
counterexamples quickly, but it is weaker than Beurling-Nyman because the
functional is not already RH-equivalent.

### Bridge C - Li Positivity As MDL Slack

This is conceptually attractive but tertiary.

Li's criterion expresses RH as nonnegativity of a sequence of coefficients
derived from the xi-function. The MDL reading would interpret each coefficient
as a slack coordinate: every coordinate of the zeta compressor must remain
nonnegative if the extremal geometry is correctly balanced.

This bridge should not lead the program yet. It needs more translation work
before "description length" is more than an analogy.

## Why EHP114 Matters Here

EHP114 gives a real calibration pattern:

```text
short symmetric description
  -> extremal boundary length
  -> defect/dispersion terms measure deviation from symmetry
```

The current local #114 work adds one useful internal lesson: the smooth-Hessian
picture is too naive near the radial boundary. The better object is stratified:
radial Puiseux behavior plus nonradial shape-cone positivity plus mixed
remainder absorption. That lesson transfers to RH only as a warning: if an RH
MDL curve exists, the hard part may be a singular or coefficient-control term,
not a smooth quadratic Hessian.

## First Concrete Experiment

Build a Beurling-Nyman MDL probe with a brutally explicit cost model.

Minimum viable version:

1. Fix three dictionaries:
   `Theta_N = {1/k : 1 <= k <= N}`,
   `Theta_N = {exp(-k h_N) : 1 <= k <= N}`,
   and a seeded randomized schedule from the Alouges-Darses-Hillion style.
2. For each dictionary, solve the finite least-squares projection problem in
   high precision.
3. Quantize coefficients to a predeclared bit grid and recompute the residual.
4. Save `epsilon_N`, active support size, coefficient-bit cost, theta-bit cost,
   Gram condition number, and residual decay.
5. Report the curve:

```text
description_bits_N vs -log2(residual_N)
```

Pass signal:

- stable residual decay under coefficient quantization,
- no single dictionary artifact dominates,
- and Gram/Hankel conditioning is explicitly reported rather than hidden.

Fail signal:

- residual improvement disappears after quantization,
- the best curve depends on an unprincipled theta schedule,
- or the Gram matrix condition number explains the whole effect.

## Artifact Contract For A Future Run

If this becomes an actual experiment run, use a new immutable experiment ID,
for example:

```text
EXP-MATH-RH-BEURLING-NYMAN-MDL-PROBE-20260506-01
```

and emit:

```text
*_RESULTS.json
*_REPORT.md
*_RESULTS.sha256
```

Do not update D1 or any public status from this note. A run would be
`SUGGESTIVE` at best unless it proves a formal implication to RH-equivalent or
zero-free-region statements.

## Next Move

Write the probe only after the dictionary and cost model are frozen in the
header. The point is not to find a pretty curve. The point is to make it hard
for the curve to be an artifact of coefficient precision, basis selection, or
Gram-matrix instability.
