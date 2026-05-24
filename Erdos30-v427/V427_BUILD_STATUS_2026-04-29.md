# Erdos30 Lean 4.27 Build Status

Date: 2026-04-29
Workspace: `erdos-experiments/Erdos30-v427`

## Version Target

`erdosproblems.com` does not appear to impose a standalone Lean toolchain pin.
Its formalisation thread discusses whether and how Lean formalizations should
count for the site/wiki, including conditional formalizations, but it does not
publish a separate Lean version requirement.

The practical compatibility target is therefore Google DeepMind
`formal-conjectures`, because that is the public Lean corpus that explicitly
uses Erdős Problems as a source and requires contributors to run `lake build`.
As of this local build, the relevant checked-out target is:

- Lean: `leanprover/lean4:v4.27.0`
- Mathlib: `v4.27.0`
- Mathlib commit resolved by Lake: `a3a10db0e9d66acbebf76c5e6a135066525ac900`

## Local Result

Created a separate v4.27 workspace instead of mutating the legacy v4.24 packet:

`/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30-v427`

The default target is now:

`Erdos30_PublicV427`

Build command:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30-v427
lake build
```

Result:

```text
Build completed successfully (7900 jobs).
```

## Public Bundle Included

The successful v4.27 aggregate imports:

- `Erdos30_Sidon_Defs`
- `Erdos30_Lindstrom`
- `Erdos30_BFR`
- `Erdos30_Singer`
- `Erdos755_BhG`
- `Erdos755_DifferenceCount`
- `Erdos755_Lindstrom`
- `Erdos755_Singer_BhG`
- `Erdos755_Complete`
- `Erdos755_B3G`
- `Erdos755_BhG_General`
- `Erdos1_DistinctSubsetSums`
- `Erdos755_HigherOrder`
- `Erdos166_SumFree`

## Porting Changes Applied

Only the v4.27 copy was changed.

- `lean-toolchain` updated from `v4.24.0` to `v4.27.0`.
- `lakefile.lean` updated from Mathlib commit `f897ebcf72cd16f89ab4577d0c826cd14afaafc7` to `v4.27.0`.
- `Erdos30_Lindstrom.scaled_range` changed from a `< succ` proof term to the direct `Nat.div_le_div_right` inequality expected by Mathlib v4.27.
- `Erdos30_BFR` diagonal-cardinality proof now rewrites through `Finset.diag_eq_filter` before applying `A.diag_card`.
- Added `Erdos30_PublicV427.lean` as the default aggregate build target.

## Known Non-default Legacy Drift

The following copied legacy/scratch targets are not part of the default v4.27
public bundle and still need separate port work if we want the entire old
packet to compile:

- `Erdos30_Complete`
- `Erdos30_OrderedElements`
- `Erdos30_SharpDiff`
- `Erdos30_difference_counting`
- `Sidon_SumCount_Fix`

The failures are Mathlib/API/tactic drift, not evidence that the v4.27 public
bundle failed.
