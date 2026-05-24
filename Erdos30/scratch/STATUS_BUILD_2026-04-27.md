# Erdos30_SpectralSidon.lean — Build Status 2026-04-27

## Toolchain
- lean-toolchain: `leanprover/lean4:v4.24.0`
- mathlib version (from lakefile): commit `f897ebcf72cd16f89ab4577d0c826cd14afaafc7`

## Sorry / Axiom Inventory (Erdos30_SpectralSidon.lean)
- sorries: **4** at lines 49, 79, 110, 142
- axioms: **1** at line 95 (`spectral_collision_30`)

Header claim ("4 explicit sorries + 1 axiom") is **ACCURATE**.

Details:
- Line 49: `sorry` — sin injectivity on `[0, π/2]` for `spectralDistanceSet_card`
- Line 79: `sorry` — finite verification (27 cases) for `spectral_sidon_small`
- Line 95: `axiom spectral_collision_30` — trig identity `sin(π/15) + sin(7π/15) = sin(2π/15) + sin(4π/15)`
- Line 110: `sorry` — instantiate `IsRealSidonSet` with four elements for `spectral_not_sidon_30`
- Line 142: `sorry` — saturation argument for `spectral_sidon_saturates`

## Build Status
- lake exe cache get: **skipped** — `.lake/` directory already present; warm incremental cache
- lake build (default targets): **PASS**
- duration: 162 seconds (incremental; Mathlib already cached)
- lake build Erdos30_SpectralSidon (explicit): **FAIL** — `error: unknown target`

## If FAIL: First 3 Compile Errors
```
error: unknown target `Erdos30_SpectralSidon`
```
No compile errors from the file itself — it was never exercised.

## Critical Question Answered
**Is `scratch/` part of the build target?**
Partially. The `lakefile.lean` defines two scratch targets:
- `Erdos30_difference_counting` (srcDir = "scratch")
- `Sidon_SumCount_Fix` (srcDir = "scratch")

`Erdos30_SpectralSidon` is **NOT** listed as any `lean_lib` target in `lakefile.lean`. It is present in `scratch/` but invisible to the lake build system entirely.

**Did `Erdos30_SpectralSidon.lean` compile (specifically)?**
**NOT EXERCISED.** No build artifact was produced for it. `lake build` (default targets) succeeded in 162 seconds but never touched this file. `lake build Erdos30_SpectralSidon` returns `error: unknown target` immediately, confirming it is not a registered lake target.

Per the H² status pipeline: **"nothing counts until COMPILED."** This file has not been compiled. It cannot receive CLEAN or COMPILED status.

## Recommendation
The file is currently invisible to the build system — it lives in `scratch/` but has no `lean_lib` entry in `lakefile.lean`. Before any status assessment is meaningful, Ken needs to add a target entry to `lakefile.lean` such as:

```lean
lean_lib Erdos30_SpectralSidon where
  srcDir := "scratch"
  roots := #[`Erdos30_SpectralSidon]
```

Only after that addition can `lake build Erdos30_SpectralSidon` run and surface real compile errors (or pass). The file contains 4 sorries and 1 axiom — all explicit and well-documented — so it is correctly labeled PROOF ARCHITECTURE, not COMPILED. It is not a candidate for promotion to `lean/` until (1) it is added to lakefile, (2) `lake build` runs without error, and (3) the 4 sorries are resolved. The axiom (`spectral_collision_30`) requires either a Mathlib trig proof or explicit axiom declaration, which it already has correctly.
