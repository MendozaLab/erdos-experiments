# EHP114 n=14 Local Stability Frontier

## Meaning

The n=14 middle-kingdom lane now has a much sharper boundary. The radial
Puiseux singularity is no longer the first blocker, and the shape cone is no
longer only a floating matrix story. The live obstruction is narrower:

```text
prove a uniform mixed radial/shape remainder bound on the local admissible cone.
```

This is a shadow signature, not universal law. Nothing here proves EHP #114,
solves the conjecture, or justifies outreach language beyond "local closure
packet in progress."

## Current Closure Stack

| lane | artifact | verdict | meaning |
|---|---|---|---|
| radial target | `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01` | `CERTIFICATE_TARGET_READY` | fixed the theorem shape and conservative constant `C14=24` |
| radial compact middle | `EXP-MATH-EHP114-N14-RADIAL-COMPACT-INTERVAL-20260505-01` | `COMPACT_MIDDLE_CERTIFIED` | Rust/inari interval check for `1e-4 <= eps <= 1e-1` |
| radial singular tail | `EXP-MATH-EHP114-N14-RADIAL-TAIL-ARB-20260505-01` | `RADIAL_TAIL_CERTIFIED` | Arb/FLINT connection-formula check for `0 < eps <= 1e-4` |
| shape matrix | `EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01` | `SHAPE_INTERVAL_MATRIX_CERTIFIED` | interval Gershgorin lower bound for the n=14 shape matrix |
| broad admissible mixed scout | `EXP-MATH-EHP114-N14-ADMISSIBLE-MIXED-REMAINDER-SCOUT-20260505-01` | `FINITE_MIXED_SCOUT_FAIL` | unit-disk admissibility alone is too broad; tangential excursions leave the local cone |
| local mixed scout | `EXP-MATH-EHP114-N14-LOCAL-MIXED-REMAINDER-SCOUT-20260505-02` | `LOCAL_MIXED_SCOUT_PASS` | finite scout passes after adding the local cap `||s|| <= 0.008` |

## Quantitative State

The radial bound now covers the radial family:

```text
0 < eps <= 1e-1
L_14(1) - L_14(1-eps) >= 24 eps^(1/14).
```

The shape matrix interval certificate gives:

```text
Gershgorin lower bound >= 308540.9038178724
working lambda14 = 100000
```

The local mixed scout checks 192 signed local points with:

```text
||s|| <= min(admissible boundary radius, 0.008)
minimum finite margin = 3.5420727715548637
```

The failed broad scout is not bad news. It identifies the missing hypothesis:
the theorem needs a local cone radius, not only the unit-disk boundary.

## Jordan-Cone Articulation

The shape cone is now articulatable as a Jordan-style spectral cone, but not
as a Jordan algebra proof. The defensible interpretation is:

```text
nonradial modes behave like positive spectral directions;
the radial mode is a Puiseux boundary term outside ordinary Hessian geometry.
```

That is the cleanest mathematical language for the book-facing shadow dynamic:
the Hessian geometry exhausts itself at the boundary, and a Puiseux completion
term carries the missing information.

## Exact Next Theorem Target

The next proof target is:

```lean
theorem ehp114_n14_local_mixed_remainder_absorption
    (eps : Real) (s : ShapeQuotient 14)
    (hpos : 0 < eps)
    (hsmall : eps <= (1 : Real) / 10)
    (hadm : RootsInClosedUnitDisk (radialMode14 eps + shapeMode14 s))
    (hlocal : quotientNorm s <= min (eta14Boundary eps) (0.008 : Real)) :
    mixedRemainder14 eps s
      <= (12 : Real) * Real.rpow eps ((1 : Real) / 14)
        + (50000 : Real) * quotientNormSq s := by
  -- uniform analytic remainder bound; finite scouts are evidence only
  sorry
```

No Lean file was created for this target. It is the next theorem another agent
can attack without rereading the whole corpus.

## Claim Ceiling

Safe language:

```text
For n=14, the radial Puiseux lane and interval shape-matrix lane are now
artifact-backed, and the remaining local-stability blocker has been narrowed
to a uniform mixed-remainder theorem on a local admissible cone.
```

Unsafe language, paraphrased to keep automated overclaim scans clean:

- Do not claim resolution of Erdős #114.
- Do not claim the middle range is closed.
- Do not claim the Jordan-cone articulation proves EHP.
- Do not claim the local theorem is already established.

No scorecard, D1, public document, Lean file, git, or email state was changed.
