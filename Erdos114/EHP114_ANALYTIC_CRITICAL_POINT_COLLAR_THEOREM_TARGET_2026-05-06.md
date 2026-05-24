# EHP114 Analytic Critical-Point Collar Theorem Target

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` branch-isolation residuals only.

Claim ceiling: this is a theorem target for local certificate work. It is not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.

## Why This Target Exists

The local length budget is still favorable:

```text
accepted length upper = 20.316451752723314
hard-cell cap         = 20.672796062619668
margin               = 0.3563443098963539
```

The blocker is not length. The blocker is branch certification: the unresolved hard-cell pieces still need a proof that each residual tube is either empty or contains a controlled graph branch with counted ownership.

The first-order axis collars, first-order normal collars, and sampled third-order collar remainder inequality all failed to certify branches. The latest third-order pilot did improve the diagnosis: 19 of 64 residual regions passed critical-point exclusion, but none passed wall separation. The first failing inequality was:

```text
sup_r |F(0,r)| + normal_remainder = 0.8371883316743337
S * lower(|F_n|)                 = 0.03309058624726584
margin                           = -0.8040977454270679
```

That says the next object should not be another box-counting run. It should be an analytic collar theorem that reduces the center-strip and wall-separation constants.

## Target Theorem Shape

Let

```text
F(x,y; u) = |p_u(x + i y)|^2 - 1
```

where `u` ranges over the root-affine parameter cell for n=14 hard subcell `(6,4)`. For a residual branch tube, choose midpoint-gradient normal/tangent coordinates:

```text
z = z0 + s n + r t
```

with normal coordinate `s`, tangent coordinate `r`, normal radius `S`, and tangent radius `R`.

For every admissible root parameter `u` and every residual tube, prove constants:

```text
m0  <= |F_n(0,0;u)|
Vn  >= sup_{|s|<=S, |r|<=R} |F_n(s,r;u) - F_n(0,0;u)|
C0  >= sup_{|r|<=R} |F(0,r;u)|
Rn  >= normal Taylor remainder from s=0 to |s|=S
Kt  >= sup_{|s|<=S, |r|<=R} |F_t/F_n|
```

such that:

```text
m0 - Vn > 0
C0 + Rn < S * (m0 - Vn)
```

Then the collar has no critical point, the two normal walls separate the zero set, and the zero set inside the collar is a certified graph over the tangent coordinate. Its length is bounded by:

```text
L <= width_tangent * sqrt(1 + Kt^2)
```

## What Must Improve Over L18

The L18 pilot used conservative interval bounds:

```text
Vn <= |F_nn| S + |F_nt| R + 0.5 T3 (S + R)^2
Rn <= 0.5 |F_nn| S^2 + (T3 / 6) S^3
T3 <= 2 |p'''| |p| + 6 |p''| |p'|
```

That was enough to pass critical exclusion on 19 regions, but wall separation failed because `C0` remained far too large. The analytic target must therefore attack one of these constants:

- reduce `C0` by proving a center-strip cancellation or one-dimensional tangent expansion instead of using a raw interval hull;
- reduce `R` or split tangent collars while preserving ownership;
- replace the whole-box root-affine interval with a root-collar constant that respects the actual residual tube geometry;
- prove a better lower bound for `m0 - Vn` that does not spend the full root-affine uncertainty budget at once.

## Proof-Facing Acceptance

A useful theorem-target artifact must state, for each processed residual region:

- the chosen collar coordinates;
- the constants `m0`, `Vn`, `C0`, `Rn`, and `Kt`;
- whether `m0 - Vn > 0`;
- whether `C0 + Rn < S * (m0 - Vn)`;
- the resulting length contribution if certified;
- the first failed inequality if not certified.

No branch may be promoted unless both critical exclusion and wall separation pass. No local hard-cell pass may be reported unless all processed residual regions are certified or excluded and the total validated length remains under `20.672796062619668`.

## Next Implementation Lane

The next diagnostic should be a theorem-target packet or Rust pilot that isolates `C0`, the center-strip term, with a one-dimensional tangent Taylor model before trying another full branch-atlas sweep. If that does not shrink the wall-separation inequality, the remaining route is a more global analytic critical-point exclusion theorem rather than further local subdivision bookkeeping.
