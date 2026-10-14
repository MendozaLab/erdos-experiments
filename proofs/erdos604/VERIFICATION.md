# Verification record

The publication sources were checked on 14 October 2026 with Lean 4.33.1
(`819816b2e0a3bf405af45ae5c7af2491d8f5bee6`, arm64 macOS) and Mathlib
`0df444a360eaa60ab8c11dca51a86af692955474`. The other eight dependency revisions
are pinned in `lake-manifest.json` and recorded in the receipts.

## Fresh direct compiler replay: passed

All 58 package modules were compiled from the final, notice-bearing sources into
a fresh output directory. Previously compiled, pinned Mathlib/dependency caches
were reused. The replay checked their Git revisions and origins; it did not rebuild
Mathlib or download dependencies. All module invocations exited 0, and all source
SHA-256 hashes were unchanged at completion.

The reproducible invocation is:

```sh
python3 scripts/replay.py --lean "$LEAN" --packages "$PACKAGES" --output "$OUTPUT"
```

Here `LEAN` names the installed Lean 4.33.1 binary, `PACKAGES` names the existing
pinned dependency directory, and `OUTPUT` is a new directory. Each module is
compiled with `autoImplicit=false`, `relaxedAutoImplicit=false`, and
`warningAsError=true`. The script records every command and compiler log.

The [direct receipt](verification/direct-receipt.json) records all 58 source hashes,
dependency pins, per-module commands, exit codes, timings, and log hashes.
[Combined compiler output](verification/direct-compiler.log) contains every
module's original output under an added module heading. Private machine paths in
the recorded command arrays are replaced with the explicitly documented `$LEAN`
and `$OUTPUT` placeholders; compiler output is unchanged.

## Transitive axiom audit: passed

`Erdos604Proof/Audit.lean` prints the transitive axioms of 27 selected declarations,
including `OAI.WeakPinned.main`, the counting and asymptotic adapters,
`almost_all_pins_from_openai`, and `erdos604_from_openai`. Every report contains
only `propext`, `Classical.choice`, and `Quot.sound`.

The publication Lean sources contain no `sorry`, `admit`, new `axiom`,
`native_decide`, or `unsafe` declaration. The final theorem is derived from the
included OpenAI proof without an endpoint assumption. This is an ordinary Lean
compiler/kernel check of this compatibility port against trusted compiled
dependencies, together with an axiom audit; it does not independently validate
the compiler binary or rebuild every dependency.

## Standalone Lake build: passed

The package's plain `lakefile.toml` and exact lockfile were tested in a separate
copy with fresh package outputs and isolated writable copies of the existing
dependency caches:

```sh
lake --wfail build OAI Erdos604Proof Erdos604Proof.Audit
```

The command exited 0. It was run through the reviewed local build-slot guard.
Every publication Lean source and each configuration file matched the tested copy
by SHA-256. The package defines no custom Lake build hooks. This check exercises
Lake configuration, imports, and all package targets using pinned cached
dependencies; it is not a cold dependency download/build test.

See the [Lake receipt](verification/lake-receipt.json) and
[complete build log](verification/lake-build.log).

## Exact current Formal Conjectures target: passed

The separate [bridge source](verification/CurrentTargetBridge.lean.txt) imports the
actual `FormalConjectures.ErdosProblems.«604»` module from
[`84cf7f8adfb0ea8d14c7d00bdf145cc9032bdf88`](https://github.com/google-deepmind/formal-conjectures/tree/84cf7f8adfb0ea8d14c7d00bdf145cc9032bdf88)
alongside this package. It proves both copied counting definitions equal to the
upstream definitions by `rfl`, then obtains the exact weak Erdős 604 proposition
with `exact erdos604_from_openai`.

The bridge exited 0 under the same strict compiler flags. Its three axiom reports
contain only the same standard axioms; the upstream placeholder proof is not used.
See the [bridge receipt](verification/current-target-bridge.json) and
[compiler output](verification/current-target-bridge.log).

To repeat that separate check, build the pinned Formal Conjectures 604 module,
copy the supplied `.lean.txt` file to `CurrentTargetBridge.lean`, and set `LEAN_PATH`
to include this package's output directory, the Formal Conjectures
`.lake/build/lib/lean` directory, and the `.lake/build/lib/lean` directories of its
pinned dependencies. From the bridge file's directory, run:

```sh
"$LEAN" -DautoImplicit=false -DrelaxedAutoImplicit=false \
  -DwarningAsError=true CurrentTargetBridge.lean
```

The text-file suffix prevents this separate cross-project check from becoming a
dependency of the standalone package.

## Limits and source review

The 47 OpenAI modules and eight adapter/application modules were independently
compared with the saved source records. The compatibility patch covers every
difference from the pinned OpenAI source, including the modification notices.
The target definition bodies are unchanged apart from namespace and notation.
Per-file provenance and attribution are recorded in [PROVENANCE.md](PROVENANCE.md).

These receipts do not claim an independent mathematical rederivation of every
argument in OpenAI's paper, a replay of the unchanged Lean 4.34.1 release, or an
official Comparator run. The exact target bridge is a direct Lean equality and
application check, rather than a claim that Comparator was used.

The proof establishes the weak, fixed-ε asymptotic bound and the stated uniform
almost-all-pins consequence. It supplies neither the sharper `n / sqrt(log n)`
bound nor quantitative control as ε varies with n.
