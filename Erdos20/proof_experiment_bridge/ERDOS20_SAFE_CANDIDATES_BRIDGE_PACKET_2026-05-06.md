# Erdos #20 Safe Candidates Bridge Packet

**Date:** 2026-05-06
**Worker:** #20-E
**Lane:** sunflower core-closure proof/experiment bridge
**Scope:** local bridge specification only. No D1, scorecard, git, public-doc, or theorem-status update.
**Required phrase:** shadow signature, not universal law.

## Bottom Line

The live blocker is semantic, not syntactic. `Erdos20CoreClosureScratch.lean` already has the right structural slot:

```lean
safeCandidates : Finset (Finset α)
safe_subset_candidates : safeCandidates ⊆ candidateExtensions
```

But that slot is not yet bound to real enumerator data. The enumerator currently reports useful per-core counts such as `candidate_count`, `local_valid_count`, and `I_core_local_bits`; those are good diagnostics, but they are not enough for Lean to check that a concrete candidate set is safe. To turn the scratch carrier into a theorem target, the enumerator must emit a candidate-level table with explicit safety flags and unsafe witnesses.

The bridge is:

```text
enumerator family/core rows
  -> candidateExtensions as concrete Finsets
  -> SafeExtension predicate
  -> safeCandidates = candidateExtensions.filter SafeExtension
  -> real theorem target about filtered safe candidates and measured closure cost
```

This remains A0 structural/protocol work. It does not prove, improve, advance, or resolve the sunflower conjecture.

## Inputs Read

- `erdos-experiments/Erdos20/formal_core_closure/ERDOS20_FORMAL_CORE_CLOSURE_TARGETS_2026-05-05.md`
- `erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean`
- `erdos-experiments/Erdos20/experimental_deepening/ERDOS20_CORE_CLOSURE_DEEPENING_PACKET_2026-05-05.md`
- `Math-Problems/proofs/archive/Sunflower_MDL.lean`
- `erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md`
- Continuity checks against existing May 5 per-core result JSON files, read-only.

## Meaning of `safeCandidates`

Fix:

- ground set `[1..n]`;
- uniform set size `w`;
- petal count `k = 3`;
- sunflower-free family `F` of `w`-sets;
- candidate core `C` of size `s`;
- candidate extension `B`.

For the proof bridge, `B` is a candidate extension through `C` when:

```text
B subset [1..n]
|B| = w
C subset B
B notin F
```

`B` is locally safe through `C` when adding `B` does not close a forbidden 3-sunflower using core `C`. Equivalently, there do not exist distinct `A1, A2 in F` such that:

```text
C subset A1, C subset A2
A1 != A2
(A1 \ C), (A2 \ C), and (B \ C) are pairwise disjoint
A1 cap A2 = C
A1 cap B = C
A2 cap B = C
```

For `k = 3`, that witness pair is the whole local obstruction. This is the exact predicate that should populate:

```lean
safeCandidates = candidateExtensions.filter (SafeExtension C F)
```

The bridge should keep both local and global notions separate:

- `local_safe_through_core`: no 3-sunflower is closed through this exact `core`.
- `global_safe`: adding `B` preserves sunflower-freeness across all possible cores.

The Lean target for the current scratch should use local safety first. Global safety is valuable later, but it mixes all cores and makes the first bridge harder than necessary.

## Required Enumerator Schema

The next enumerator output should use schema version `erdos20.safe_candidates.v1` and emit one JSON file per run. It can include summaries, but the proof bridge requires the candidate-level rows.

```json
{
  "schema_version": "erdos20.safe_candidates.v1",
  "experiment_id": "EXP-MATH-ERDOS20-SAFE-CANDIDATES-YYYYMMDD-NN",
  "generated_at_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "runner": {
    "name": "sunflower-core-closure-rust",
    "source_path": "erdos-experiments/Erdos20/rust_core_closure/src/main.rs",
    "binary_sha256": "<sha256-or-null>",
    "command": "<exact command>",
    "deterministic": true
  },
  "parameters": {
    "n": 8,
    "w": 3,
    "k": 3,
    "ground_set": [1, 2, 3, 4, 5, 6, 7, 8],
    "ground_set_indexing": "1-based",
    "target_m": 8,
    "jamming_m_star": 8,
    "window_offset": 0
  },
  "family_rows": [
    {
      "family_id": "F-00000001",
      "family_hash_sha256": "<canonical sorted set-list hash>",
      "m": 8,
      "sets": [
        {"set_id": "S-0001", "elements": [1, 2, 3], "mask_u64": 7}
      ],
      "sunflower_free_k3": true,
      "sunflower_free_check": {
        "method": "exhaustive_triple_scan",
        "checked_triples": 56,
        "forbidden_witness": null
      },
      "core_rows": [
        {
          "core_id": "C-0001",
          "core_elements": [1, 2],
          "core_mask_u64": 3,
          "core_size_s": 2,
          "petal_size": 1,
          "candidate_count": 6,
          "local_safe_count": 4,
          "local_unsafe_count": 2,
          "global_safe_count": 3,
          "candidate_rows": [
            {
              "candidate_id": "B-0001",
              "elements": [1, 2, 5],
              "mask_u64": 19,
              "petal_elements": [5],
              "petal_mask_u64": 16,
              "contains_core": true,
              "uniform_card_w": true,
              "unused_by_family": true,
              "local_safe_through_core": false,
              "global_safe": false,
              "unsafe_reason": "closes_3_sunflower_through_core",
              "unsafe_witness_pairs": [
                {
                  "a1_set_id": "S-0002",
                  "a1_elements": [1, 2, 3],
                  "a2_set_id": "S-0003",
                  "a2_elements": [1, 2, 4],
                  "pairwise_intersection_core": true,
                  "petals_pairwise_disjoint": true
                }
              ]
            }
          ]
        }
      ]
    }
  ],
  "summary_rows": [
    {
      "n": 8,
      "w": 3,
      "k": 3,
      "m": 8,
      "core_size_s": 2,
      "family_count": 48389005,
      "core_count": 28,
      "candidate_count_total": 6968016720,
      "local_safe_count_total": 5675439840,
      "I_core_local_bits": 0.296016,
      "claim_ceiling": "shadow signature, not universal law; no theorem or lower-bound progress"
    }
  ]
}
```

For large runs, `family_rows` may be emitted as a bounded proof sample plus a separate `summary_rows` table, but the theorem bridge requires at least one complete family/core/candidate block. Counts alone are not a proof object.

## Field Requirements

The fields below are mandatory for the proof bridge.

| Field | Level | Purpose |
|---|---|---|
| `schema_version` | run | Allows future schema evolution without silently changing theorem meaning. |
| `experiment_id` | run | Immutable run identity for artifact traceability. |
| `runner.command` | run | Reproducibility hook; records exact enumerator command. |
| `runner.binary_sha256` | run | Binds output to the executable when available. |
| `parameters.n` | run | Ground-set size. |
| `parameters.w` | run | Uniform set size. |
| `parameters.k` | run | Must be `3` for the first theorem target. |
| `parameters.ground_set` | run | Concrete universe to load into Lean as `Fin n` or `Nat` values. |
| `parameters.target_m` | run | Family size being checked. |
| `parameters.jamming_m_star` | run | Experimental context; not used as a theorem hypothesis. |
| `family_id` | family | Stable row key for a concrete family. |
| `family_hash_sha256` | family | Canonical hash of sorted set masks. |
| `sets[].elements` | family | Concrete `Finset Nat` payload. |
| `sets[].mask_u64` | family | Fast consistency check against elements. |
| `sunflower_free_k3` | family | Required precondition for extension safety. |
| `sunflower_free_check.method` | family | Must state how the boolean was certified. |
| `core_elements` | core | Concrete core `C`. |
| `core_mask_u64` | core | Fast consistency check against `core_elements`. |
| `core_size_s` | core | Must equal `core_elements.length`. |
| `petal_size` | core | Must equal `w - core_size_s`. |
| `candidate_count` | core | Must equal `candidate_rows.length` for proof samples. |
| `local_safe_count` | core | Must equal rows with `local_safe_through_core = true`. |
| `candidate_rows[].elements` | candidate | Concrete extension `B`. |
| `candidate_rows[].petal_elements` | candidate | Must equal `elements \ core_elements`. |
| `contains_core` | candidate | Must be true for all candidate extensions. |
| `uniform_card_w` | candidate | Must be true for all candidate extensions. |
| `unused_by_family` | candidate | Must be true for all candidate extensions. |
| `local_safe_through_core` | candidate | The boolean that instantiates `SafeExtension C F B`. |
| `unsafe_witness_pairs` | candidate | Required nonempty witness list when `local_safe_through_core = false`. |
| `pairwise_intersection_core` | witness | Checkable certificate for the obstruction. |
| `petals_pairwise_disjoint` | witness | Checkable certificate for the obstruction. |

Optional but useful:

- `global_safe`, for future all-core extension safety.
- `automorphism_orbit_id`, if symmetry reduction is used.
- `family_weight`, if rows represent an orbit rather than one literal family.
- `excluded_reason`, if a row is omitted from proof samples.

## Invariants

The enumerator must enforce these before writing `PASS`.

1. Ground set invariant: every element in every set/core/candidate is in `[1..n]`.
2. Encoding invariant: each `mask_u64` is exactly the bit encoding of `elements` under the stated 1-based convention.
3. Uniformity invariant: every family set and candidate has cardinality `w`.
4. Family-size invariant: `m = sets.length`.
5. Family uniqueness invariant: family set masks are pairwise distinct.
6. Candidate uniqueness invariant: candidate masks are pairwise distinct within each core row.
7. Core invariant: `core_elements` is a subset of every candidate in that core row.
8. Candidate extension invariant: every candidate contains the core, has size `w`, and is not already in the family.
9. Petal invariant: `petal_elements = candidate_elements \ core_elements` and `petal_size = w - core_size_s`.
10. Baseline family invariant: `sunflower_free_k3 = true` is backed by an exhaustive triple scan or by an explicitly named certified method.
11. Local unsafe witness invariant: if `local_safe_through_core = false`, at least one witness pair is present and checkable.
12. Local safe invariant: if `local_safe_through_core = true`, exhaustive search over all distinct pairs in `F` found no witness pair through the core.
13. Count invariant: `candidate_count = local_safe_count + local_unsafe_count`.
14. Proof-sample invariant: for complete family/core samples, `candidate_count = candidate_rows.length`.
15. Summary invariant: summary totals must be computed from the same predicate as the candidate rows, not from a second implementation.

The most important anti-artifact rule is number 15: validation cannot silently reimplement a different safety predicate.

## Lean-Shaped Predicate

The current scratch can be extended by adding a predicate like this. This is a theorem target shape, not a claim that the code below has been compiled today.

```lean
namespace Erdos20FormalCoreClosure

open Finset

variable {alpha : Type*} [DecidableEq alpha]

def Petal (C A : Finset alpha) : Finset alpha :=
  A \ C

def CandidateExtension (w : Nat) (C : Finset alpha)
    (F : Finset (Finset alpha)) (B : Finset alpha) : Prop :=
  B.card = w /\ C subset B /\ B notin F

def ClosesThreeSunflowerThroughCore (C : Finset alpha)
    (F : Finset (Finset alpha)) (B : Finset alpha) : Prop :=
  exists A1 in F, exists A2 in F,
    A1 != A2 /\
    C subset A1 /\ C subset A2 /\ C subset B /\
    A1 inter A2 = C /\
    A1 inter B = C /\
    A2 inter B = C

def SafeExtension (w : Nat) (C : Finset alpha)
    (F : Finset (Finset alpha)) (B : Finset alpha) : Prop :=
  CandidateExtension w C F B /\
  not ClosesThreeSunflowerThroughCore C F B

end Erdos20FormalCoreClosure
```

In real Lean syntax, `subset`, `notin`, `exists ... in ...`, and `inter` must be replaced by Mathlib forms (`⊆`, `∉`, bounded membership hypotheses, and `∩`). They are written ASCII here to keep the packet portable across agents.

## Lean-Shaped Theorem Target

The first real theorem target should be deliberately small:

```lean
theorem safeCandidates_eq_filter_and_subset
    {alpha : Type*} [DecidableEq alpha]
    (w : Nat)
    (C : Finset alpha)
    (F : Finset (Finset alpha))
    (candidateExtensions : Finset (Finset alpha))
    (safeCandidates : Finset (Finset alpha))
    (hSafe :
      safeCandidates =
        candidateExtensions.filter
          (fun B => SafeExtension w C F B)) :
    safeCandidates subset candidateExtensions := by
  intro B hB
  rw [hSafe] at hB
  exact (Finset.mem_filter.mp hB).1
```

Once that compiles, the next theorem target binds the experiment count:

```lean
theorem coreClosureCost_from_filtered_safeCandidates
    {alpha : Type*} [DecidableEq alpha]
    (w : Nat)
    (S : CoreClosureState (alpha := alpha))
    (hSafe :
      S.safeCandidates =
        S.candidateExtensions.filter
          (fun B => SafeExtension w S.core S.family B)) :
    CoreClosureCost S =
      S.candidateExtensions.card - S.safeCandidates.card := by
  rfl
```

That second statement is intentionally bookkeeping. The real value is not the proof difficulty; it is that `safeCandidates` becomes a table generated by the same predicate used in Lean.

## Proof-to-Experiment Bridge Procedure

1. Run the enumerator with `--emit-safe-candidates-v1` on one cheap calibration target, preferably `w=3,n=5,k=3,target_m=4` before repeating `w=3,n=8`.
2. Validate schema and invariants locally.
3. Convert one complete `family_rows[*].core_rows[*]` block into Lean literals:
   - `family : Finset (Finset Nat)`
   - `core : Finset Nat`
   - `candidateExtensions : Finset (Finset Nat)`
   - `safeCandidates : Finset (Finset Nat)`
4. Prove `safeCandidates = candidateExtensions.filter (SafeExtension w core family)`.
5. Instantiate `CoreClosureState` with that equality and prove subset/cost bookkeeping.
6. Only after the literal sample compiles, consider summary-level claims about `I_core_local_bits`.

This separates three things that should not be conflated:

- the concrete proof object for one family/core sample;
- the aggregate enumerator summary over many families/cores;
- the Leg-4 classifier interpretation.

## Claim Ceiling

Safe statement:

> The #20 lane now has a specified bridge from per-core enumerator rows to the Lean `safeCandidates` slot. The next proof target is a filtered-subset theorem over concrete enumerator samples. The status remains A0: shadow signature, not universal law.

Unsafe statements:

- The sunflower conjecture has been advanced.
- A new lower bound has been proved or improved.
- The Maxwell analogy is a law for sunflowers.
- The May 5 per-core signal is Leg-4 PASS evidence.
- A Lean theorem is COMPILED unless it is built in the current Lean project.

## Files Changed

- Created: `erdos-experiments/Erdos20/proof_experiment_bridge/ERDOS20_SAFE_CANDIDATES_BRIDGE_PACKET_2026-05-06.md`

No other files were written.

## Commands Run

Read-only inspection:

```bash
rg -n "Erdos #20|Erdős #20|sunflower|safeCandidates|proof_experiment_bridge|formal_core_closure" /Users/kenbengoetxea/.codex/memories/MEMORY.md
sed -n '1,220p' /Users/kenbengoetxea/.agents/skills/math-problems-manager/SKILL.md
sed -n '1,260p' erdos-experiments/Erdos20/formal_core_closure/ERDOS20_FORMAL_CORE_CLOSURE_TARGETS_2026-05-05.md
sed -n '1,260p' erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean
sed -n '1,260p' erdos-experiments/Erdos20/experimental_deepening/ERDOS20_CORE_CLOSURE_DEEPENING_PACKET_2026-05-05.md
ls -la erdos-experiments/Erdos20/proof_experiment_bridge
sed -n '1,280p' Math-Problems/proofs/archive/Sunflower_MDL.lean
sed -n '1,260p' erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md
find erdos-experiments/Erdos20 -maxdepth 3 -type f | sort | sed -n '1,220p'
sed -n '1,220p' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json
sed -n '1,220p' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-W3N8-20260505-02_RESULTS.json
sed -n '1,260p' erdos-experiments/Erdos20/SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md
rg -n "safe|core|closure|jamming|m_star|unsafe|blocked|candidate|signature|universal" erdos-experiments/Erdos20 -g '*.{md,json,py,rs}'
```

Write command:

```bash
mkdir -p erdos-experiments/Erdos20/proof_experiment_bridge
```

Manual edit:

```text
apply_patch added this bridge packet.
```

## Next Concrete Step

Add `--emit-safe-candidates-v1` to the Rust enumerator in a future worker lane, but do not write that code from this packet unless ownership explicitly expands to `erdos-experiments/Erdos20/rust_core_closure/`. The first emitted proof sample should be small enough to translate by hand into Lean literals and compile before scaling back to `w=3,n=8`.
