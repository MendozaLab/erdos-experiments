# EHP114 Tao Middle-Kingdom Closure Packet

Experiment: `EXP-MATH-EHP114-TAO-MIDDLE-KINGDOM-20260505-01`

## Meaning

This packet turns the gap between the DOI-backed n <= 14 certificate and Tao's sufficiently-large-n theorem into a concrete attack surface.

The target is not brute-force computation past n = 20. The target is to make Tao's asymptotic proof quantitative enough to splice with Rust/inari interval certificates, while replacing global B&B with radial, shape-cone, and remainder interval lemmas.

Preserve the book-facing phrase: shadow signature, not universal law.

## Source

- arXiv source: `https://arxiv.org/e-print/2512.12455`
- arXiv abstract: `https://arxiv.org/abs/2512.12455`
- source SHA-256: `ab792bd8fe806a985e44f692ccb8305cf386677509846979bb8ff35a276814a6`
- local DOI anchor: `10.5281/zenodo.19480329`

## What Tao Gives

Tao proves EHP for sufficiently large degree. The paper states that all implied constants are effectively computable, so the remaining problem is finite, but the paper does not optimize the resulting numerical bound.

The practical consequence is a middle range:

```text
n = 1,2        known / classical
n = 3..14      DOI-backed Rust/inari IEEE-1788 certificate
n = 15..20     plausible finite-certificate extension, not yet safe
n = 21..N0-1   middle kingdom
n >= N0        Tao, after explicit constant extraction
```

## Dependency Landmarks

| label | source lines | role | middle-gap use |
|---|---:|---|---|
| `main-thm` | 185-192 | Tao high-degree ceiling | Defines the upper end of the finite middle range, but no practical N0. |
| `erem` | 203-204 | Extremizer normalization | Required if the finite certificate is to match Tao's proof variables. |
| `stokes` | 357-370 | Area-to-length functional | This is the natural bridge for an interval-compressed certificate. |
| `out` | 612-624 | p0 reference value | This is where DOI L* intervals and Tao's asymptotic reference meet. |
| `geomcontrol` | 1190-1220 | First geometry control | First constant-chase choke point: every hidden constant enters later. |
| `x2-lem` | 1275-1278 | Small-field error term | Quantifying the 'mu sufficiently small' choice is a tractable first chase. |
| `x3-lem` | 1305-1311 | Large-psi error term | Pairs with X2 to optimize mu explicitly instead of asymptotically. |
| `lemni` | 1346-1363 | Total size bound | Separates global shape control from final local deficit. |
| `inside` | 1419-1421 | Inner annulus estimate | Ancestor of the radial Puiseux interval target. |
| `annulus` | 1471-1473 | Intermediate region estimate | Middle-region term must be made quantitative for any finite splice. |
| `outside` | 1516-1518 | Outer tip estimate | This is where endpoint/tip constants enter the final comparison. |
| `ets` | 1640-1643 | Critical-point collapse | This is the first final-section gate toward uniqueness. |
| `pots` | 1706-1708 | Total-size collapse | This makes the final split local around p0. |
| `inside-2` | 1739-1751 | Final inner deficit | Closest analytic cousin of the radial Puiseux certificate. |
| `annulus-2` | 1781-1784 | Final middle deficit | The shape-cone interval bound should target this term. |
| `outside-again` | 1831-1833 | Final outer deficit | Remainder absorption must keep this error below the inner/middle gain. |

## The Attack Order

1. **Constant chase Tao's final proof.**
   Start at `inside-2`, `annulus-2`, and `outside-again`, then walk backward through `pots`, `ets`, `inside`, `annulus`, `outside`, and `geomcontrol`. The goal is not a beautiful bound; the first goal is any explicit N0.

2. **Interval-harden the radial term first.**
   This corresponds to the inner-deficit lane. It should reuse the n = 14 DOI certificate as calibration and target a fixed-n theorem for the radial hypergeometric/Puiseux deficit.

3. **Then harden the shape cone.**
   This corresponds to the intermediate-region deficit. The theorem should prove that nonradial perturbations cannot restore the lost length after the radial term is controlled.

4. **Then absorb the mixed remainder.**
   This corresponds to keeping the outer-region error and local mixed terms below the inner plus middle gains.

5. **Only then extend finite n.**
   Reconcile the n = 15 and n = 16 zero-eval artifacts before claiming them. Treat n = 20 as a diagnostic endpoint until a new versioned Rust/inari certificate exists.

## First Concrete Theorem Target

```lean
theorem ehp114_fixed_n_radial_puiseux_interval
    (n : Nat) (eps : Real)
    (hn : n = 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10000) :
    C14 * Real.rpow eps ((1 : Real) / 14)
      <= D14 (radialMode14 eps) := by
  -- hypergeometric connection formula + interval constants
  sorry
```

This is intentionally fixed-n. If n = 14 closes, general n becomes a
parameterization problem. If n = 14 does not close, n = 20 is premature.

## Claim Ceiling

Safe: `finite certified computation n <= 14`, `Tao proves sufficiently large n`,
and `we now have a concrete middle-range closure map`.

Unsafe: saying the middle range is closed, saying n = 20 is certified, or
claiming a Lean/formal proof of the analytic bridge.

## Verification

- arXiv source fetched and parsed.
- `lemniscate.tex` landmarks extracted by label.
- No scorecard, D1, public deployment, or Lean status was mutated.

