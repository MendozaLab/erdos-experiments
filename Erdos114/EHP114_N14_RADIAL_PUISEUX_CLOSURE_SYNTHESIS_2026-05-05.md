# EHP114 n=14 Radial Puiseux Closure Synthesis

## Meaning

The n=14 radial/Puiseux lane has moved from a numerical signature to a
two-part interval certificate for the radial family

```text
p_a(z) = z^14 - a,  a = 1 - eps.
```

The certified radial statement is:

```text
0 < eps <= 1e-1
L_14(1) - L_14(1-eps) >= 24 eps^(1/14).
```

This is a shadow signature, not universal law. It does not prove EHP #114 and
does not settle nonradial local stability. It removes the radial singularity as
the first blocker in the middle-kingdom closure strategy.

## Evidence Stack

| layer | artifact | status | meaning |
|---|---|---|---|
| finite anchor | `EXP-MM-EHP-007-n14-inari_RESULTS.json` | DOI-backed Rust/inari certificate | n=14 is the largest byte-reconciled finite proof anchor in the current corpus |
| theorem target | `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01` | `CERTIFICATE_TARGET_READY` | fixed the constant, domain split, and Lean-shaped theorem |
| compact middle | `EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01` | `COMPACT_MIDDLE_CERTIFIED` | IEEE-1788 inari interval check for `1e-4 <= eps <= 1e-1` |
| singular tail | `EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01` | `RADIAL_TAIL_CERTIFIED` | Arb/FLINT connection-formula check for `0 < eps <= 1e-4` |

The compact certificate minimum margin is `0.4683310451937004`. The tail
connection-formula coefficient lower bound is
`25.304169858286159276467082057586698228234056420266`, leaving margin
`1.3041698582861592764670820575866982282340564202659` over the required
constant `24`.

## Why This Matters

Before this packet, the Hessian-fantasy lane had a real ambiguity: the radial
direction is singular, so an ordinary quadratic Hessian story was the wrong
local model. The Puiseux exponent `1/14` is now the right local coordinate for
the radial defect, and the radial defect has a certificate path that no longer
depends on brute-force branch-and-bound.

That is the middle-kingdom pattern we wanted: the finite certificate supplies
the anchor, Tao supplies the large-n endpoint, and the shadow term supplies the
local completion mechanism.

## What Is Still Open

The next blocker is not radial. It is the shape-cone/remainder splice:

```text
show that nonradial perturbations cannot erase the certified
24 eps^(1/14) radial deficit inside the n=14 local cone.
```

The next theorem target should be:

```lean
theorem ehp114_n14_radial_shape_remainder_splice
    (eps : Real) (s : ShapeQuotient 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hs : quotientNorm s <= eta14 eps) :
    (12 : Real) * Real.rpow eps ((1 : Real) / 14)
      <= D14 (radialMode14 eps + shapeMode14 s) := by
  -- radial certificate + shape cone + mixed remainder absorption
  sorry
```

The constant `12` is deliberately half the radial constant. It leaves half the
radial deficit as budget for mixed-term absorption. It is a theorem target, not
a claimed Lean theorem.

## Next Attack Order

1. Build the n=14 shape-cone interval matrix around the quotient space
   orthogonal to rotation and scaling.
2. Prove a conservative lower bound for the shape quadratic term on that cone.
3. Bound the mixed radial/shape remainder by at most half the radial deficit.
4. Only after the n=14 local cone is closed, repeat the same scaffold at n=15
   and n=16 after reconciling their zero-evaluation provenance.

No scorecard, D1, public document, Lean file, git, or email state was changed.
