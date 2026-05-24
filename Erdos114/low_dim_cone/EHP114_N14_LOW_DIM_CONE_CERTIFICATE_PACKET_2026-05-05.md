# EHP114 n=14 Low-Dimensional Cone Certificate Packet

Date: 2026-05-05

## Meaning

The current n=14 route is no longer trying to preserve transported shape
positivity after radial contraction. The admissible spectral run showed real
shape softening, so the useful target is the scalar reserve theorem:

```text
D14(eps, s) >= 12 eps^(1/14)
```

on the epsilon-scaled admissible cone. This remains a shadow signature, not
universal law. It is a local n=14 certificate route, not an Erdős #114 closure.

## New Finite Artifact

This packet adds one low-dimensional finite probe:

```text
EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01
```

Files:

- `ehp114_n14_low_dim_cone_box_probe.py`
- `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_ORACLE_INPUT.json`
- `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_ORACLE_OUTPUT.json`
- `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_RESULTS.json`
- `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_REPORT.md`
- `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_RESULTS.sha256`

The run reuses:

- `ehp114_n14_interval_taylor_m14_packet.py` for interval arithmetic helpers.
- `ehp114_n14_eps_scaled_spectral_direction_search.py` for the lowest
  admissible-Taylor midpoint eigendirections.
- `probe_ehp114_tensor_cone_scaffold.py` for quotient basis and coefficient
  construction.
- The existing `ehp114_batch_interval_lengths` oracle binary for interval
  lemniscate lengths.

## Exact Subspace Tested

For each epsilon in:

```text
eps in {1e-4, 1e-3, 1e-2, 1e-1}
```

the tested subspace is:

```text
V_eps = span(u0(eps), u1(eps))
```

where `u0(eps)` and `u1(eps)` are the two lowest midpoint eigenvectors of the
24-dimensional admissible-stencil Taylor matrix from:

```text
EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01
```

The sampled points are:

```text
s = eps^(1/28) * (a u0(eps) + b u1(eps))
sqrt(a^2 + b^2) <= 0.014
a,b in {-0.014, -0.007, 0, 0.007, 0.014}
```

That is a two-dimensional coefficient-disk grid inside the lowest spectral
subspace, not a continuous box proof and not a full 24-dimensional cone proof.

## Result

The finite grid had 52 candidate points before enforcing root admissibility.
The root-admissibility filter left 16 oracle-evaluated points. Every evaluated
admissible point passed the scalar reserve target.

| eps | admissible/evaluated | inadmissible | all evaluated pass | min margin |
|---:|---:|---:|:---:|---:|
| 1e-04 | 1/13 | 12 | yes | 7.58178533831 |
| 1e-03 | 1/13 | 12 | yes | 8.94388640959 |
| 1e-02 | 1/13 | 12 | yes | 10.5029403545 |
| 1e-01 | 13/13 | 0 | yes | 12.0595154164 |

Interpretation: at `eps = 0.1`, this is a real mixed two-dimensional grid
test inside the low spectral span. At the three smaller epsilons, the proposed
eta cap is too large for this particular 2D grid under root admissibility:
only the base point remained admissible. So the packet supports the route, but
it does not yet certify the eta `0.014` two-dimensional cone uniformly over
all four epsilon rows.

## eps = 0.1 Dense Box-Candidate Follow-Up

A denser follow-up was run on the only epsilon row that behaved like a true
two-dimensional local box candidate:

```text
EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01
```

This run fixes `eps = 0.1`, keeps the same two-dimensional spectral span
`span(u0,u1)`, and densifies the coefficient disk

```text
sqrt(a^2 + b^2) <= 0.014
```

using a `17 x 17` grid plus cell centers for cells fully inside the disk.

Result:

```text
status = EPS01_LOW_DIM_BOX_CANDIDATE_PASS
evaluated oracle points = 361
failure count = 0
minimum point margin = 12.05951541635811
candidate cells whose four corners and center all pass = 164
minimum sample margin among candidate cells = 12.062197693283327
continuous certificate = false
```

This is the strongest local finite evidence so far for the interval-box lane.
It is still not a continuous certificate because the current Rust oracle
evaluates concrete coefficient points, not coefficient boxes.

## What This Certifies

It certifies only this finite statement:

```text
For the 16 root-admissible grid points emitted in the oracle input,
the interval oracle lower bound gives D14(eps,s) >= 12 eps^(1/14).
```

It does not certify every point in the 2D disk, every point in a 2D box, or
any point outside the root-admissible filtered set. It also does not certify
the full quotient cone.

## What Would Promote This To A Genuine Certificate

A genuine low-dimensional certificate needs three extra pieces:

1. Root-admissibility boxes for the `(a,b)` coefficient domain, not just
   pointwise max-root-radius tests.
2. Interval lower bounds for `D14(eps,s) - 12 eps^(1/14)` over those same
   `(a,b)` boxes.
3. A rule that either fixes the eigenspace on epsilon intervals or replaces
   epsilon-dependent eigenvectors with an analytic basis whose variation is
   bounded.

The immediate certificate upgrade is therefore not another wider point sample.
It is an interval-box proof on a smaller admissible coefficient region, starting
with the `eps = 0.1` two-dimensional span where all 13 grid points survived
admissibility.

After the dense follow-up, the next blocker is sharper:

```text
prove a rigorous cell-wise variation bound for D14 over the (u0,u1)
coefficient boxes, or upgrade the Rust oracle to accept interval coefficient
boxes directly.
```

The dense packet shows there is ample sampled margin at `eps = 0.1`; the proof
gap is no longer point coverage. The proof gap is continuous enclosure between
sampled points.

## One-Cell Variation Probe

The next step was run on the strongest sampled cell from the dense packet:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01
```

Selected coefficient rectangle:

```text
u0 in [-0.00175, 0]
u1 in [0, 0.00175]
```

Result:

```text
status = ONE_CELL_ROOT_BOX_PASS_EMPIRICAL_VARIATION_PASS
root box admissible by affine radius bound = true
max affine root-radius upper bound = 0.9936164637647286
oracle points inside selected cell = 289
sampled deficit failures = 0
minimum sampled margin = 12.101568027485637
maximum sampled margin = 12.103239066384457
max neighbor margin slope = 1.4657307798578183
required Lipschitz upper bound for this margin = 9779.543788829036
empirical margin after one cell-radius variation = 12.099754278181432
continuous deficit certificate = false
```

This is a meaningful upgrade: the root-admissibility side is enclosed over the
entire selected cell by a simple affine root-radius bound. The deficit side is
still not continuous-certified; it is a dense variation diagnostic.

The resulting next theorem is extremely concrete:

```text
For the selected eps = 0.1 coefficient rectangle, prove a rigorous Lipschitz
constant L < 9779.543788829036 for the scalar margin function, or directly
enclose the deficit by interval coefficient arithmetic.
```

The empirical local slope is around `1.47`, so the required rigorous bound is
very loose relative to the observed behavior. That makes this cell the best
first target for an actual interval/Cauchy certificate.

## Continuous-Cell Interval Attempts

The first continuous-cell attempt deliberately started with the selected
`eps = 0.1` coefficient rectangle from the one-cell variation probe:

```text
u0 in [-0.00175, 0]
u1 in [0, 0.00175]
```

The naive coefficient-interval oracle failed:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-INTERVAL-COEFF-ORACLE-20260505-01
status = ONE_CELL_INTERVAL_COEFF_ORACLE_FAIL
length upper = 2847.9689848799253
margin lower = -2827.2961888173054
active cells = 36920
definite-case cells = 0
uncertain-corner cells = 36920
```

Interpretation: this is dependency blowup, not mathematical evidence against
the cell. Passing coefficient boxes through elementary symmetric coefficients
destroys the usable enclosure.

The root-affine oracle then evaluated the polynomial directly as

```text
p(z) = prod_i (z - r_i)
```

where each root is an affine interval in the `(u0,u1)` coefficient rectangle.
That attempt still failed on the whole cell, but narrowed the interval bound by
two orders of magnitude:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-ORACLE-20260505-01
status = ONE_CELL_ROOT_AFFINE_ORACLE_FAIL
length upper = 44.5971758610202
margin lower = -23.92437979840053
active cells = 640
definite-case cells = 104
uncertain-corner cells = 536
```

This identified the viable representation: intervalize the roots, not the
coefficients, then reduce box width.

## First Continuous Grid-Oracle Cell Certificate

Subdivision of the same selected cell gives the first PASS for a continuous
root-affine interval certificate of the current marching-squares oracle
functional:

| run | subcells | status | failures | max length upper | min margin lower |
|---|---:|---|---:|---:|---:|
| `SUBDIV4` | 16 | `ONE_CELL_ROOT_AFFINE_SUBDIV4_FAIL` | 16 | 22.99013026541943 | -2.3173342027997634 |
| `SUBDIV8` | 64 | `ONE_CELL_ROOT_AFFINE_SUBDIV8_PASS` | 0 | 18.110795101362747 | 2.5620009612569206 |

The `SUBDIV8` artifact is:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01
```

Checksum verification:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.json: OK
```

This is materially stronger than point sampling. It certifies, for this one
selected coefficient cell and this current marching-squares oracle functional,
that every one of the 64 root-affine interval subcells has lower margin above
zero.

Claim ceiling:

```text
This is a continuous grid-oracle cell certificate, not exact lemniscate
certification, not a Lean theorem, and not a proof of Erdős #114.
```

## Rust Port Decision

The exploratory code can remain Python while the certificate geometry is still
moving. The hardened engine should be Rust.

Reason: the root-affine interval representation is now stable enough to
deserve a deterministic, fast, reviewable implementation. Rust should own the
next certificate layer:

```text
input: epsilon, spectral basis, coefficient box, subdivision, grid extent/resolution
output: root-admissibility interval bound, length upper bound, deficit lower bound,
        margin lower bound, per-subcell audit rows, sha256-bound result JSON
```

The Rust engine should reproduce the Python `SUBDIV8` PASS before it is allowed
to generalize to more cells, more epsilon rows, or exact lemniscate-length
enclosures.

## Rust Reproduction Gate

The Rust reproduction gate has now passed in the existing `erdos-114` Rust
crate:

```text
EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01
status = RUST_ROOT_AFFINE_SUBDIV8_REPRO_PASS
failure count = 0
rows = 64
max length upper = 18.110795101366655
min margin lower = 2.5620009612530126
reproduction gate pass = true
```

Saved Rust output:

```text
erdos-experiments/scripts/erdos-114/rust-root-affine-subdiv8-output-v2/
```

The local smoke rerun into `/tmp` also reproduced the same pass:

```text
status = RUST_ROOT_AFFINE_SUBDIV8_REPRO_PASS
failure count = 0
max length upper = 18.110795101366655
min margin lower = 2.5620009612530126
```

Meaning: the selected-cell continuous grid-oracle certificate is no longer a
Python-only exploratory artifact. The hardened Rust path reproduces it. The
claim ceiling does not change: this is still not exact lemniscate
certification, not a Lean theorem, and not a proof of Erdős #114.

## Exact-Length Lift Budget

The exact-length lift packet records the budget needed to promote this
grid-oracle cell certificate to proof-grade exact/validated length:

```text
EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01
status = BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE
uniform exact length cap = Lstar_lower - 12 * eps^(1/14)
                           = 20.672796062619668
worst marching length upper = 18.110795101362747
allowed additive exact-length error = 2.5620009612569206
```

The current bottleneck is now sharply isolated:

```text
prove L_exact(C) <= L_ms(C) + E(C)
with E(C) <= margin_lower(C) for every accepted root-affine subcell.
```

Equivalently, replace the marching-squares oracle functional by a validated
exact lemniscate-length enclosure.

The first Rust bridge diagnostic on the worst accepted subcell `(6,4)` did not
close this gap:

```text
EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
status = BRIDGE_REGULARITY_INTERVAL_UNRESOLVED
z-subdivision = 16
candidate level boxes = 95071
regularity unresolved boxes = 13560
sum normal-drift error candidate = 149.62095139919995
sum relative-length error candidate = 423.1887593268329
available exact-length bridge budget = 2.5620009612530126
```

Meaning: naive ambient interval boxes for `p'` are not sharp enough. The exact
length lift now needs curve-following charts, Bernstein/affine arithmetic, or a
coarea/Crofton bound with constants that do not require proving `p'` separation
on every ambient subbox.

## Next Theorem Target

```lean
theorem ehp114_n14_eps_scaled_cone_deficit
    (eps : Real) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk14 eps s)
    (hcone : quotientNorm s <= eta0_14 * Real.rpow eps ((1 : Real) / 28)) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= totalDeficit14 eps s := by
  -- direct total-deficit interval/Cauchy certificate target
  sorry
```

The low-dimensional bridge theorem before the full statement is:

```text
For a fixed epsilon slab and the two-dimensional spectral span V_eps,
certify D14(eps,s) >= 12 eps^(1/14) on an admissible coefficient box.
```

## Claim Ceiling

Safe claim:

```text
At n=14, the scalar epsilon-scaled route now has finite admissible evidence
along axes, along the lowest spectral rays, and on a sampled two-dimensional
grid in the lowest spectral span. The new grid run is strongest at eps = 0.1
and remains only pointwise at smaller epsilons because most proposed grid
points leave the root-admissible domain.
```

Overclaim bans:

- Do not state that Erdős #114 is closed or complete.
- Do not state that local stability is established.
- Do not state that the full cone is certified.
- Do not state that eta `0.014` is a uniform analytic radius.
- Do not state that epsilon `1/28` is the final cone law.
- Do not present sampled grid evidence as a Lean theorem.
- Do not use outreach language before an interval-box certificate exists.

## Added Verification Command

Run from
`/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/low_dim_cone`:

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-LOW-DIM-BOX-CANDIDATE-20260505-01_RESULTS.sha256
```

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-VARIATION-PROBE-20260505-01_RESULTS.sha256
```

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-INTERVAL-COEFF-ORACLE-20260505-01_RESULTS.sha256
```

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-ORACLE-20260505-01_RESULTS.sha256
```

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV4-20260505-01_RESULTS.sha256
```

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.sha256
```

## Verification Commands

Run from `/Users/kenbengoetxea/container-projects/apps/H2/Math` unless noted.

```bash
python3 erdos-experiments/Erdos114/low_dim_cone/ehp114_n14_low_dim_cone_box_probe.py
```

Run from `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/low_dim_cone`:

```bash
shasum -a 256 -c EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_RESULTS.sha256
```

Run from `/Users/kenbengoetxea/container-projects/apps/H2/Math`:

```bash
test -f erdos-experiments/Erdos114/low_dim_cone/EHP114_N14_LOW_DIM_CONE_CERTIFICATE_PACKET_2026-05-05.md
```
