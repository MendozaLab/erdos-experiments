# Weak pinned distances: an Erdős 604 proof integration

This Lean package derives the weak Erdős 604 statement from OpenAI's weak pinned
planar-distance theorem:

```text
For every ε > 0, for all sufficiently large n,
  n^(1−ε) ≤ pinnedDistinctDistances n.
```

The function minimizes, over all planar sets of `n` points, the largest number of
distinct positive distances seen from a point of that set. Its definitions are
copied from the referenced Formal Conjectures revision and placed in a separate
namespace for comparison without declaration conflicts.

The main declaration is
[`erdos604_from_openai`](Erdos604Proof/ExternalApplication.lean). It has no endpoint
assumption: it applies the imported, proved `OAI.WeakPinned.main`.

The same file proves `almost_all_pins_from_openai`: for every `ε > 0` and `δ > 0`,
one cutoff works for every sufficiently large planar configuration. More than
`(1−δ)n` of its pins each see more than `n^(1−ε)` distinct positive distances.

## Source and scope

The mathematical result is attributed to OpenAI's
[*The weak pinned planar distance theorem*](https://github.com/openai/math/blob/1581d4cbdd69a0faf34b2089614fe2a814316a6a/preprints/The-weak-pinned-planar-distance-theorem-September-23-2026/paper.pdf).
OpenAI already supplies a downstream pinned-distance consequence in `Pins.lean`.
This package contributes a compatibility port, an independently written counting
adapter, and an exact integration with the weak Erdős 604 formulation. It does not
claim an original solution of the mathematical problem.

The checked toolchain is Lean **4.33.1**, with Mathlib pinned to
`0df444a360eaa60ab8c11dca51a86af692955474`. This differs from the original OpenAI
release's Lean 4.34.1 / Mathlib
`d13f23b723b8a846827a245b89c10fc7d3f11612`.

- `OAI/Geometry/PinnedDistances/` contains the 47-module dependency closure of
  OpenAI's endpoint, with documented compatibility repairs in 26 files.
- `Erdos604Proof/` contains the finite counting and uniform asymptotic adapter,
  target definitions, application, boundary examples, and axiom audit.
- The sharper `n / sqrt(log n)` lower bound is outside this package's scope.
- The unchanged release and its official Comparator workflow have not been replayed
  by this package's local verification.

Full source pins, hashes, licensing, and changes are in [PROVENANCE.md](PROVENANCE.md).
Verification evidence is described in [VERIFICATION.md](VERIFICATION.md).

## Reproduce with Lake

With the pinned toolchain and dependencies available, run from this directory:

```sh
lake --wfail build OAI Erdos604Proof Erdos604Proof.Audit
```

The package configuration is a plain `lakefile.toml` with no custom build hooks.
The lockfile records the exact dependency revisions. Standard dependency setup
may fetch the pinned official dependencies if they are not already present.

## Reproduce with existing compiled dependencies

The optional direct compiler script uses a pre-existing pinned dependency cache.
It performs no dependency resolution, downloads, installation, or Lake execution:

```sh
python3 scripts/replay.py \
  --lean /path/to/lean-4.33.1/bin/lean \
  --packages /path/to/existing/.lake/packages \
  --output /path/to/new/replay-output
```

It checks the compiler version and dependency Git revisions/origins, compiles all
package modules into a fresh output directory, confirms unchanged source hashes,
and audits the selected endpoint and adapter declarations. Full compiler logs and
exact commands are saved in that output directory. Existing dependency caches are
read and never modified by the script.

## Attribution and assistance

OpenAI's proof and the Formal Conjectures definitions retain their source
attribution. The compatibility work, adapter, and packaging were prepared under
Ken Bengoetxea's direction with Codex assistance. Compiler checks and axiom audits
are distinguished from independent mathematical review and Comparator validation.

The source package is distributed under [Apache-2.0](LICENSE); see [NOTICE](NOTICE).
