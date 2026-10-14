# Provenance and changes

| Component | Source revision or version |
|---|---|
| OpenAI source | `1581d4cbdd69a0faf34b2089614fe2a814316a6a` |
| Original OpenAI Lean | `leanprover/lean4:v4.34.1` |
| Original OpenAI Mathlib | `d13f23b723b8a846827a245b89c10fc7d3f11612` |
| Checked compatibility port | Lean `4.33.1` |
| Checked Mathlib | `0df444a360eaa60ab8c11dca51a86af692955474` |
| Local checked-port record | `31cfe229796310c109785a29ba6174e6b39c51a7` |
| Copied Formal Conjectures target definitions | `57f2e0cf2b374563dcb42bdd9adec97aa7e27517` |

The local checked-port record identifies the pre-publication work; it is not
represented as an upstream OpenAI or Formal Conjectures commit.

## OpenAI proof

The original source is
[`lean/OAI/Geometry/PinnedDistances`](https://github.com/openai/math/tree/1581d4cbdd69a0faf34b2089614fe2a814316a6a/lean/OAI/Geometry/PinnedDistances).
This package includes the 47-module transitive closure of `Main.lean`. The source
`Pins.lean` is outside that closure and is not included.

The compatibility port changes 26 files for the older installed Mathlib APIs:

- use the available `Measure.map_smul` form;
- supply probability-measure instances and measurability arguments;
- replace conditional rewrite names with the supported `if_pos`, `if_neg`, and
  `dif_pos` forms;
- expose a negative-tail membership goal before linear arithmetic;
- prove pointwise function equalities explicitly in two quotient-limit arguments.

The original proof review found the theorem/lemma statements and data declarations
unchanged. One conditional rewrite is inside the existing `retained.disjoint`
proof field; its data fields and required proposition are unchanged. `Model.lean`
and `Main.lean` are byte-identical to the original source.

The 21 unmodified OpenAI modules are byte-identical to the original and checked
port. Each of the 26 modified modules carries a prominent compatibility notice;
removing that exact comment prefix gives byte-identical checked-port source.
The complete, portable [compatibility diff](provenance/compatibility.patch) shows
every difference from the pinned original. The
[source manifest](provenance/source-manifest.json) records upstream, checked-port,
and publication SHA-256 hashes for each file.

The original release's executable Lake configuration, dependency cloning/patching
hooks, unrelated proofs, and Comparator infrastructure are not included or run.

## Adapter and target definitions

The copied `distinctDistancesFrom` definition comes from
[`FormalConjecturesForMathlib/Geometry/Metric.lean`](https://github.com/google-deepmind/formal-conjectures/blob/57f2e0cf2b374563dcb42bdd9adec97aa7e27517/FormalConjecturesForMathlib/Geometry/Metric.lean).
The `pinnedDistinctDistances` definition comes from
[`FormalConjectures/ErdosProblems/604.lean`](https://github.com/google-deepmind/formal-conjectures/blob/57f2e0cf2b374563dcb42bdd9adec97aa7e27517/FormalConjectures/ErdosProblems/604.lean).
They exclude self-distance with `erase`, use the maximum pinned count in each
configuration, and take the natural-number infimum over all configurations of the
specified cardinality. The adapter proves that this infimum is attained.

Packaging makes only these changes to the previously checked eight adapter and
application modules:

1. Module imports gain the `Erdos604Proof.` prefix.
2. The copied target definitions move to `Erdos604Proof.Target`; references to
   them are qualified or opened accordingly. Their bodies are preserved.
3. `#print axioms` commands move to a dedicated `Erdos604Proof.Audit` module.
4. The target-definition header records its two source paths and namespace change.
   An application docstring describes the public package rather than the former
   local application harness; trailing whitespace/empty lines are normalized.
5. Two library-root import files provide standard Lake targets.

The theorem `erdos604_from_openai` retains its name, quantifiers, and proof body,
with the renamed target definition in its conclusion. Its proof applies the
actual imported OpenAI theorem. No endpoint hypothesis, placeholder proof,
additional mathematical axiom, or native decision procedure is introduced.

## Licensing

The exact Apache-2.0 license from the pinned OpenAI repository is included as
`LICENSE`. No upstream `NOTICE` file was present in the saved source snapshot;
this package's `NOTICE` identifies its sources and modifications. Per-file notices
identify the modified OpenAI files. The copied Formal Conjectures definitions
retain their copyright/license attribution and identify their namespace change.

Mathlib and its dependencies are not vendored as source or binaries. Their exact
source repositories and revisions are recorded in `lake-manifest.json`; each
retains its own license. This package does not claim to audit licenses for a
redistributed binary bundle, because no such bundle is supplied.
