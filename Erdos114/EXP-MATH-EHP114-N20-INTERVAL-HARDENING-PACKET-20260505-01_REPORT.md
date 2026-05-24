# EXP-MATH-EHP114-N20-INTERVAL-HARDENING-PACKET-20260505-01 Report

## Status

This packet takes the successful n=20 Hessian/Puiseux triad and narrows it into
proof obligations. It is not a proof of EHP114 and not a publication packet.
The claim ceiling remains: this is a shadow signature, not universal law.

## Verdict

- Packet verdict: `CERTIFICATE_TARGET_READY`
- Source: `EXP-MATH-EHP114-N20-HESSIAN-TRIAD-20260505-01`
- Current blocker: `Turn the diagnostic constants into interval arithmetic lemmas and prove the mixed remainder absorption bound.`

The useful fantasy is now precise: do not run another ordinary Hessian sweep.
Prove the radial Puiseux interval bound, prove the quotient shape-cone lower
bound, then absorb the mixed remainder.

## Candidate Constants

| layer | diagnostic lower evidence | working constant |
|---|---:|---:|
| radial Puiseux | min sampled ratio 40.1621701231 at eps 0.1 | C20 = 32 |
| radial tail | min tail ratio 41.5922788546 for eps <= 1e-4 | C20_tail = 40 |
| shape cone | Gershgorin lower 545638.890084; floating eig min 601487.934073 | lambda20 = 250000 |

The constants are deliberately conservative. They are not certified constants
until interval arithmetic replaces the floating diagnostic input.

## Lean-Shaped Targets

### Radial

```lean
theorem ehp114_n20_radial_puiseux_interval
    (eps : Real) (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10000) :
    (32 : Real) * Real.rpow eps ((1 : Real) / 20)
      <= D20 (radialMode20 eps) := by
  -- exact hypergeometric interval bound for the radial family
  sorry
```

### Shape Cone

```lean
theorem ehp114_n20_shape_cone_interval
    (s : ShapeQuotient 20) :
    (250000 : Real) * quotientNormSq s
      <= shapeQuadratic20 s := by
  -- interval matrix lower bound for the quotient shape block
  sorry
```

### Combined Local Certificate

```lean
theorem ehp114_n20_stratified_local_certificate
    (r : Real) (s : ShapeQuotient 20)
    (hr_pos : 0 < abs r) (hr_small : abs r <= delta20)
    (hs_small : quotientNorm s <= eta20) :
    (16 : Real) * Real.rpow (abs r) ((1 : Real) / 20)
      + (125000 : Real) * quotientNormSq s
      <= D20 (radialMode20 r + shapeMode20 s) := by
  -- radial Puiseux + shape cone + mixed remainder absorption
  sorry
```

## Mixed Remainder Obligation

For |r| <= delta20 and ||s|| <= eta20, prove R20(r,s) <= 0.5 * (C20 * |r|^(1/20) + lambda20 * ||s||^2).

This is the single mathematical blocker. The diagnostics already separate the
radial and nonradial positive terms; the closure theorem needs a proof that
mixed terms cannot erase that positivity inside the local cone.

## Source Boundary

No scorecard, D1, public document, git, email, CLAUDE.md, or AGENTS.md was
changed.
