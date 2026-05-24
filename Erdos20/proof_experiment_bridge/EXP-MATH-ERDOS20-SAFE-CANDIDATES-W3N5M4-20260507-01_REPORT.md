# Erdős #20 Safe-Candidates Bridge — Calibration Run W3N5M4

**Experiment ID:** EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01
**Date:** 2026-05-07
**Worker:** A (sunflower core-closure lane)
**Schema:** `erdos20.safe_candidates.v1`
**Claim ceiling:** A0. shadow signature, not universal law.

## Bottom Line

A new mode `--emit-safe-candidates-v1` was added to the existing Rust enumerator at `erdos-experiments/Erdos20/rust_core_closure/src/main.rs`. The mode emits one complete `family_rows` block conforming to the `erdos20.safe_candidates.v1` schema specified in the bridge packet. This run encodes one concrete sunflower-free family at the calibration target `w=3, n=5, k=3, target_m=4`, computes candidate-level local-safety flags with checkable witness pairs, and enforces invariants 1 through 15 inline in Rust. The default per-core summary path is unchanged.

This is a structural plumbing run. It does not prove, advance, improve, or resolve the sunflower conjecture. It only specifies the data the proof bridge needs.

## Files Created

- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/proof_experiment_bridge/EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01_RESULTS.json`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/proof_experiment_bridge/EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01_RESULTS.sha256`
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/proof_experiment_bridge/EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01_REPORT.md`

Modified (in-place enumerator extension):
- `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/rust_core_closure/src/main.rs`

## Runner

| Field | Value |
|---|---|
| Source | `erdos-experiments/Erdos20/rust_core_closure/src/main.rs` |
| Binary | `target/release/sunflower-core-closure` |
| Binary sha256 | `51361e17f96f26ae99b6452878bc14a4ff383a5e77226841e228f3dae588c1bb` |
| Build | `cargo build --release` |
| Determinism | true (no RNG, ascending lexicographic family search) |

## Command

```bash
./target/release/sunflower-core-closure \
  --n 5 --w 3 --target-m 4 \
  --emit-safe-candidates-v1 \
  --out /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/proof_experiment_bridge/EXP-MATH-ERDOS20-SAFE-CANDIDATES-W3N5M4-20260507-01_RESULTS.json
```

## Parameters

| Field | Value |
|---|---|
| n | 5 |
| w | 3 |
| k | 3 |
| ground_set | [1, 2, 3, 4, 5] |
| ground_set_indexing | 1-based |
| target_m | 4 |
| jamming_m_star | 4 |
| window_offset | 0 |

## Family Chosen

The Rust backtracker walks site indices in ascending order. For w=3, n=5 the site index 0 is `{1,2,3}` (1-based), index 1 is `{1,2,4}`, etc. The first sunflower-free family of size 4 it returns is the all-3-subsets-of-{1,2,3,4} family:

| set_id | elements (1-based) | mask_u64 |
|---|---|---|
| S-0001 | [1, 2, 3] | 7 |
| S-0002 | [1, 2, 4] | 11 |
| S-0003 | [1, 3, 4] | 13 |
| S-0004 | [2, 3, 4] | 14 |

`sunflower_free_check.method = exhaustive_triple_scan, checked_triples = 4, forbidden_witness = null`. Note element 5 is unused by F.

## Summary Counts (Aggregated From Candidate Rows)

| core_size_s | channel_count | candidate_count | local_safe_count | local_unsafe_count | I_core_local_bits |
|---|---|---|---|---|---|
| 0 | 1 | 6 | 6 | 0 | -0.000000 |
| 1 | 5 | 18 | 18 | 0 | -0.000000 |
| 2 | 10 | 18 | 12 | 6 | 0.584963 |

Totals across all emitted cores: candidate_count=42, local_safe=36, local_unsafe=6.

For s=0 the empty core admits every w-set not in F as a candidate; for k=3 a single empty-core obstruction would require three pairwise-disjoint w-sets in F, and our family has no triple of pairwise-disjoint sets, so all 6 candidates at s=0 are locally safe. For s=1 every petal pair across F shares ground-set element {1,2,3,4} so there is also no local obstruction. For s=2 the unsafe rows are exactly the cores {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4} — each a 2-subset of {1,2,3,4} that lies in two F-sets, with a single candidate B that completes the 3-sunflower through that core (specifically B = core ∪ {5}).

## Per-Core Breakdown (16 emitted core rows)

| core_id | core_elements | s | candidate_count | local_safe | local_unsafe |
|---|---|---|---|---|---|
| C-0001 | [] | 0 | 6 | 6 | 0 |
| C-0002 | [1] | 1 | 3 | 3 | 0 |
| C-0003 | [2] | 1 | 3 | 3 | 0 |
| C-0004 | [3] | 1 | 3 | 3 | 0 |
| C-0005 | [4] | 1 | 3 | 3 | 0 |
| C-0006 | [5] | 1 | 6 | 6 | 0 |
| C-0007 | [1, 2] | 2 | 1 | 0 | 1 |
| C-0008 | [1, 3] | 2 | 1 | 0 | 1 |
| C-0009 | [1, 4] | 2 | 1 | 0 | 1 |
| C-0010 | [1, 5] | 2 | 3 | 3 | 0 |
| C-0011 | [2, 3] | 2 | 1 | 0 | 1 |
| C-0012 | [2, 4] | 2 | 1 | 0 | 1 |
| C-0013 | [2, 5] | 2 | 3 | 3 | 0 |
| C-0014 | [3, 4] | 2 | 1 | 0 | 1 |
| C-0015 | [3, 5] | 2 | 3 | 3 | 0 |
| C-0016 | [4, 5] | 2 | 3 | 3 | 0 |

Cores at s=2 that lie in two F-sets emit exactly one unsafe candidate; cores that touch element 5 are saturated with safe candidates because element 5 is unused by F.

## Invariant Checks

All 15 invariants from the bridge packet are enforced inline in Rust (panic on violation) and re-verified by an external Python `json.load` pass over the written file.

| # | Invariant | Status |
|---|---|---|
| 1 | Ground-set range `elements ⊆ [1..n]` | PASS |
| 2 | `mask_u64` matches `elements` bit-encoding | PASS (mask built from `elements`) |
| 3 | Uniformity `|set| = w` and `|candidate| = w` | PASS |
| 4 | Family size `m = sets.length` | PASS |
| 5 | Family-set masks pairwise distinct | PASS |
| 6 | Candidate masks pairwise distinct within each core row | PASS |
| 7 | `core_elements ⊆ candidate.elements` for every candidate | PASS |
| 8 | Candidate is w-uniform and not in family | PASS |
| 9 | `petal_elements = candidate \\ core` and `petal_size = w − s` | PASS |
| 10 | `sunflower_free_k3` certified by exhaustive triple scan | PASS (4 triples, no witness) |
| 11 | Unsafe candidate carries nonempty checkable witness pair list | PASS (each witness rechecked: pic + pdj) |
| 12 | Safe candidate has empty witness list (exhaustive search closed) | PASS |
| 13 | `candidate_count = local_safe_count + local_unsafe_count` per core | PASS |
| 14 | `candidate_count = candidate_rows.length` per core | PASS |
| 15 | Summary totals derived from candidate rows, not a parallel implementation | PASS (Rust aggregates rows; Python re-aggregates and matches) |

## JSON Re-Read Sanity Check

The written JSON parses cleanly under `python3 -c "json.load(open(...))"`. Aggregating candidate rows externally yields the same `candidate_count_total = 42`, `local_safe_total = 36`, `local_unsafe_total = 6` as the Rust-emitted `summary_rows`. This confirms invariant 15 from outside the runner.

## Default Mode Continuity

The pre-existing W3N8 stdout summary path was re-run after the source change and still emits the same `core_rows` shape (`m=8, site_count=56, core_count=24, core_rows length=3`). No regression in that path.

## Blockers / Partial Completion

None. Task completed in full. The emitted file is ready for downstream Worker C consumption (Abbott-Hansen-Sauer denominator computation) and for hand-translation of one core-row block into Lean literals per the bridge packet's "Proof-to-Experiment Bridge Procedure".

## Token Embargo Confirmation

The embargoed token (per Ken's 2026-04-24 rule) does not appear anywhere in the modified Rust source, the emitted JSON, this report, or the run command. Verified by case-insensitive grep on all four artifacts.

## Required Phrase

shadow signature, not universal law.

## Claim Ceiling

A0. The status of Erdős #20 is unchanged. No theorem is COMPILED. No lower bound is proved or improved. This run only emits a candidate-level proof object that the Lean `safeCandidates` slot can be bound against in a future translation step.
