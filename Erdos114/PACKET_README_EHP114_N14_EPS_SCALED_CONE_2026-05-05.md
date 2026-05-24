# EHP114 n=14 Epsilon-Scaled Cone Packet

Date: 2026-05-05

## Meaning

This packet records the current #114 local-stability route at n=14 after the
interval Taylor / M14 split verdict. The strong theorem that transported shape
positivity remains uniformly positive after radial contraction is no longer the
right target. The evidence now points to a weaker and cleaner closure route:

```text
on an epsilon-scaled root-admissible cone, the total deficit itself stays above
12 eps^(1/14).
```

The phrase to preserve for book-facing language is: shadow signature, not
universal law.

## What Changed

The admissible-only spectral rerun confirms that negative shape-softening
directions persist even when all Taylor stencil points are root-admissible.
That removes the main worry that the earlier spectral failure was only a
non-admissible-stencil artifact.

The epsilon-scaled axis and spectral-direction searches then show that every
evaluated root-admissible point passed the scalar reserve target. The strongest
sampled eta in the mixed spectral-direction run was `0.014`.

The low-dimensional follow-up now adds the first continuous cell certificate
for the current marching-squares oracle functional. On the selected
`eps = 0.1` spectral-span coefficient cell

```text
u0 in [-0.00175, 0]
u1 in [0, 0.00175]
```

the naive coefficient-interval oracle fails by dependency blowup, but the
root-affine interval oracle passes after `8 x 8` subdivision:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01
status = ONE_CELL_ROOT_AFFINE_SUBDIV8_PASS
failure count = 0
minimum margin lower bound = 2.5620009612569206
maximum length upper bound = 18.110795101362747
```

This is a continuous grid-oracle cell certificate, not exact lemniscate
certification, not a Lean theorem, and not a proof of Erdős #114.

That cell certificate has now been reproduced by the Rust hardening path:

```text
EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01
status = RUST_ROOT_AFFINE_SUBDIV8_REPRO_PASS
failure count = 0
rows = 64
maximum length upper = 18.110795101366655
minimum margin lower = 2.5620009612530126
```

The exact-length lift packet also extracts the proof-grade bridge budget:

```text
EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01
status = BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE
allowed additive exact-length error over worst marching upper = 2.5620009612569206
```

This isolates the next bottleneck: replace the marching-squares oracle
functional by a validated exact lemniscate-length enclosure, or prove a
conservative error theorem connecting the two.

The first worst-subcell bridge diagnostic has now been run:

```text
EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01
status = BRIDGE_REGULARITY_INTERVAL_UNRESOLVED
z-subdivision = 16
regularity unresolved boxes = 13560
candidate error sums = 149.62 / 423.19, both above the 2.562 budget
```

So the next proof move is not more naive `z`-box subdivision. It is a sharper
level-set chart or Bernstein/affine/coarea bridge that avoids ambient interval
dependency in `p'`.

## Claim Ceiling

Safe claim:

```text
At n=14, the current local closure target has been redirected from transported
shape-Hessian positivity to an epsilon-scaled scalar reserve theorem, with
finite admissible axis and spectral-direction evidence through eta = 0.014.
```

Unsafe claims:

- Do not say Erdős #114 is solved.
- Do not say local stability is proved.
- Do not say the full cone is certified.
- Do not say epsilon^(1/28) is the true cone law.
- Do not say this is ready for outreach as a proof.
- Do not present the `SUBDIV8` cell certificate as exact lemniscate
  certification.

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

## Included Artifacts

- `QUORUM_SYNTHESIS_EHP114_INTERVAL_TAYLOR_M14_2026-05-05.md`
- `EHP114_N14_LOCAL_MIXED_REMAINDER_THEOREM_TARGET_2026-05-05.md`
- `EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01_REPORT.md`
- `EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01_REPORT.md`
- `EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01_REPORT.md`
- Matching `*_RESULTS.json` and `*_RESULTS.sha256` files
- Scripts used to produce the three new experiment artifacts
- `Ehp114LocalMixedRemainderScratch.lean`
- `low_dim_cone/EHP114_N14_LOW_DIM_CONE_CERTIFICATE_PACKET_2026-05-05.md`
- `low_dim_cone/EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_*`
- `cauchy_bounds/EHP114_N14_ANALYTIC_CAUCHY_BOUND_PACKET_2026-05-05.md`
- `exact_length_lift/EHP114_N14_EXACT_LENGTH_LIFT_PACKET_2026-05-05.md`
- `exact_length_lift/EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_*`
- `EHP114_N14_ROOT_AFFINE_RUST_PORT_SPEC_2026-05-05.md`
- Rust output under `erdos-experiments/scripts/erdos-114/rust-root-affine-subdiv8-output-v2/`
- Rust bridge diagnostic under
  `erdos-experiments/scripts/erdos-114/bridge-diagnostic-worst-subcell-output-z16/`

## Verification Commands

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114
shasum -a 256 -c EXP-MATH-EHP114-N14-ADMISSIBLE-SPECTRAL-TAYLOR-20260505-01_RESULTS.sha256
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS-SCALED-CONE-AXIS-SEARCH-20260505-01_RESULTS.sha256
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS-SCALED-SPECTRAL-DIRECTION-SEARCH-20260505-01_RESULTS.sha256
cd low_dim_cone
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01_RESULTS.sha256
cd ../exact_length_lift
shasum -a 256 -c EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.sha256
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114
cargo build --release --bin ehp114_n14_root_affine_cell
cd rust-root-affine-subdiv8-output-v2
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01_RESULTS.sha256
cd ../bridge-diagnostic-worst-subcell-output-z16
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01_RESULTS.sha256
```
