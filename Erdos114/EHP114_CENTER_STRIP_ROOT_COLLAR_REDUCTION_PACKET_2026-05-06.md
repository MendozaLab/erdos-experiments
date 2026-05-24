# EHP114 Center-Strip Root-Collar Reduction Packet

Date: 2026-05-06

Scope: local EHP114 / n=14 hard-cell `(6,4)` residual branch certification only.

Claim ceiling: proof-facing route design only. This is not a proof of Erdos #114, not a global n=14 certificate, not an exact lemniscate-length certificate, and not a claim upgrade. This remains a shadow signature, not universal law.

## Executive Diagnosis

The current blocker is no longer "need more compute." It is the center-strip term in the wall-separation inequality.

The L18 third-order collar pilot tested 64 normal-collar residual regions. It certified no branch, but it separated the obstacle:

```text
processed residual regions       = 64
third-order certified branches   = 0
critical-exclusion passes        = 19
wall-separation passes           = 0
critical-exclusion failures      = 45
wall-separation failures         = 19
accepted length upper            = 20.316451752723314
hard-cell cap                    = 20.672796062619668
```

For the first wall failure, the inequality was:

```text
C0 + Rn = 0.8371883316743337
S*m     = 0.03309058624726584
margin  = -0.8040977454270679
```

where:

```text
C0 = sup_r |F(0,r)|
Rn = normal Taylor remainder
m  = lower bound for |F_n|
S  = normal wall distance
```

The normal remainder was tiny there:

```text
C0 = 0.8366399140044845
Rn = 0.0005484176698492166
```

So the failure is not third-order derivative size. The first wall failure is almost entirely `C0`.

## Two Different Wall-Failure Classes

The 19 wall failures split into two classes.

### Class A: Center-Strip Reducible

There are 14 regions where:

```text
allowed_C0 = S*m - Rn > 0
```

Those regions could pass wall separation if the center-strip bound is sharpened enough.

The raw interval hull `sup_r |F(0,r)|` is too pessimistic because it spends all root-affine dependency at once. A root-centered or parameter-adapted tangent frame should replace it with a Taylor form:

```text
C0 <= E0 + |F_t(0,0)| R + 0.5 |F_tt| R^2 + (T3/6) R^3
```

The key opportunity is to make `F_t(0,0)` vanish or become small by centering at a validated branch point and using a parameter-adapted tangent direction. If `F_t(0,0)` is controlled away, the quadratic-plus-third-order bound would already pass 10 of the 19 wall failures under the L18 constants.

Concrete first-region comparison:

```text
raw C0 from L18                = 0.8366399140044845
allowed C0                    = 0.032542168577416625
quadratic C0 if F_t(0,0)=0    = 0.0010832541787435003
```

That is the proof-facing opening. The wall condition does not need a new global length method; it needs the center strip to be centered on the branch with a moving/validated tangent frame.

### Class B: Normal-Budget Failed

There are 5 regions where:

```text
allowed_C0 = S*m - Rn <= 0
```

In those regions, even perfect center-strip cancellation cannot pass the current wall inequality. They need a different repair first:

- smaller normal collar radius;
- sharper normal remainder;
- better lower bound on `m`;
- or a different analytic collar scale.

The largest offender family in this class is source ownership key `3707:1`, where the normal remainder already exceeds the wall budget in all four refined subregions.

## Critical-Exclusion Failures

The remaining 45 regions fail before wall separation:

```text
|F_n(0,0)|_lower - variation_bound(F_n) <= 0
```

These are not center-strip failures. They are derivative-variation failures. The least bad one is close:

```text
source key       = 4169:1
split path       = root/y0/y0
center |F_n|     = 0.45866683260596713
variation bound  = 0.51703377189946
margin           = -0.05836693929349285
```

The worst failures have large Hessian/third-derivative variation and should not be attacked by the same method as Class A.

## Proof-Facing Lemma To Try Next

### Lemma: Branch-Centered Moving-Frame Collar

For each residual region in Class A, prove the existence of a validated center point `z_*(u)` and a parameter-dependent frame `(n(u), t(u))` such that:

```text
F(z_*(u); u) = 0
t(u) perpendicular to grad F(z_*(u); u)
|F_n| >= m > 0 on the collar
```

Then the tangent center strip satisfies:

```text
F(0,0;u) = 0
F_t(0,0;u) = 0
C0 <= 0.5 M_tt R^2 + (T3/6) R^3
```

and the wall condition becomes:

```text
0.5 M_tt R^2 + (T3/6) R^3 + Rn < S*m
```

This is the clean theorem-shaped replacement for the raw `sup_r |F(0,r)|` interval hull.

## Why This Is Orthogonal To More Rust

Rust can certify the final constants, but Rust cannot fix the current inequality by running longer. The present interval calculation is asking a fixed midpoint frame to control a branch that should be described in a moving branch-centered frame.

The next mathematical decision is therefore:

```text
stop bounding a thick center strip;
prove a branch-centered moving-frame collar.
```

Only after that theorem shape is fixed should a new Rust diagnostic be written.

## Next Artifact Gate

The next executable artifact should not process all 4398 unresolved slab branches. It should target only the Class A wall-failure subset first and report:

```text
processed_class_a_count
validated_center_point_count
moving_frame_tangent_zero_count
quadratic_center_strip_pass_count
remaining_class_a_unresolved_count
```

Acceptance for that pilot:

```text
quadratic_center_strip_pass_count > 0
```

Failure condition:

```text
validated branch-centered moving frame cannot be established on the easiest key 2286:0 residuals
```

If even `2286:0` cannot be centered, the correct next route is not subdivision. It is a global analytic critical-point exclusion theorem for the hard-cell residual geometry.
