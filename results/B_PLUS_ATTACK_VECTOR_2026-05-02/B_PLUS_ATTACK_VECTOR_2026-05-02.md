# B+ Attack Vector Decision Memo - 2026-05-02

## Verdict

The conditional claim is not yet true as a public claim. It is a credible attack vector.

Current tier: B.

Conditional B+ path: clear the blockers below, then claim a systematic deployment instrument, not theorem-resolution parity with autonomous solvers.

## What B+ Requires

Operational B+ means all five conditions hold:

1. A nontrivial Lean portfolio compiles cleanly: lake build PASS, zero sorry, current Mathlib, registry refreshed from D1.
2. The morphism framework has a public-safe definition surface.
3. At least one morphism/collider lane has reproducible evidence plus a demoted failure/control.
4. The public packet distinguishes theorem progress, formalization progress, and platform/instrument progress.
5. A current SOTA comparison shows the platform differentiator without claiming solver superiority.

This packet clears parts of 1, 2, 3, 4, and 5 locally. It does not clear public release, direct D1 refresh, or publication.

## Lean Compile Truth

Local builds run:

| Workspace | Command | Result | Toolchain | Mathlib | Public status |
|---|---|---:|---|---|---|
| `Lean4/transdimensional-painter` | `lake build` | PASS | Lean 4.27.0 | a3a10db0... | COMPILED local, zero sorry, explicit axioms |
| `export-packets/erdos30` | `lake build` | PASS | Lean 4.24.0 | f897ebcf... | COMPILED local, zero sorry, explicit axioms |
| `erdos-experiments/private-attacks/lean4` | `lake build` | PASS | Lean 4.24.0 | f897ebcf... | COMPILED local, zero sorry, explicit axioms |
| `erdos-experiments/Erdos30` | `lake build` | PASS | Lean 4.27.0 | a3a10db0... | MIXED: zero-sorry core targets plus one scratch file with 4 sorrys |

Direct Wrangler D1 queries failed because the session had no `CLOUDFLARE_API_TOKEN`. The live API reported 21 lake-pass, zero-sorry artifacts, but strict external copy should rerun the direct D1 query before quoting a number.

What compiled: multiple local Lean workspaces and the active Erdos30 certificate workspace.

What failed: the direct D1 query; the active Erdos30 scratch file `scratch/Erdos30_SpectralSidon.lean` remains non-public because it has 4 `sorry`s.

Public number: do not quote a final direct-D1 count from this session. The audit observation is that the live API reported 21 lake-pass zero-sorry artifacts on 2026-05-02.

## Morphism Surface

Public-safe definition:

> A putative morphism is a structural correspondence tested by triangulated evidence: literature support, a target result, transfer to a different problem, a collider or failure test, and a formal checkpoint when available.

Safe public lane: map, check, demote, formalize.

Internal/quarantined lane: scoring algorithm, transfer-operator construction details, dimensional decomposition, calibration heuristics, auto-discovery recipes, internal feature rankings.

## Collider Evidence

Evidence packet includes:

| Lane | Evidence role | Result |
|---|---|---|
| Sidon #30 field response | pass support | finite field/frontier response packet exists and sidecar verifies |
| B_2[2] transfer | cross-problem support | 11/11 split rows in n=20..30, exact finite counts |
| Sum-free #166 | negative control | 11/11 formula matches, 0 split rows |
| Hypercube exp8_v2 | demoted failure | `LEG4_FAIL`, all three criteria failed |
| Chiral exp9_v2 | ambiguous/demoted | `CHIRAL_CLASS_REJECTED`, C3 passed but C1/C2 failed |

Evidence ceiling: finite laboratory. This does not prove #30, does not close an Erdos problem, and does not make any morphism discovered fact.

## SOTA Crosswalk

Current SOTA check demotes the comparison:

`Boris Alexeev-style systematic deployment`: PARTIAL.

The external AI-Erdos ecosystem is stronger on theorem resolution and autonomous Lean formalization. The Atlas can compare only on systematic deployment: repeated corpus routing, explicit evidence manifests, demotions, and formalization targets.

The differentiator is not solver strength. It is an auditable map-check-formalize instrument with a falsification layer.

## B+ Blockers

1. Rerun direct D1 Lean count with Wrangler auth and reconcile the API/scorecard/local build surfaces.
2. Quarantine or remove non-zero-sorry scratch targets from any public build claim.
3. Publish only the public-safe morphism definition surface after the IP gate or explicit approval.
4. Build the outward packet that separates A-axis theorem progress, B-axis formalization progress, and C-axis platform/instrument progress.
5. Keep the SOTA comparison dimensional: no claim of autonomous theorem-resolution parity.

## What Can Be Publicly Claimed

ErdosAtlas can be described as a systematic deployment instrument that maps putative structural correspondences, tests them with reproducible finite packets, demotes failures, and routes survivors toward formal proof targets.

It can say it has local and API-observed compiled Lean artifacts, but external count copy must be re-pulled from D1 before use.

It can say the evidence packet includes both pass-support and demotion/control examples.

## What Must Remain Internal

- Scoring algorithm internals.
- Transfer-operator construction details.
- Dimensional decomposition and calibration heuristics.
- Scout-swarm/autodiscovery mechanisms.
- Claims that #30 is solved, proved, or closed.
- Claims that morphisms are discovered facts.
- Claims that the Collider proves theorems.

## Narrow Next Move

Make the B+ claim true by clearing the registry/publication gate, not by strengthening rhetoric:

1. Restore D1 read auth and rerun the exact direct D1 Lean count.
2. Mark the active Erdos30 scratch sorry file as non-public or move it out of public build surfaces.
3. Prepare the public-safe morphism page using `MORPHISM_PUBLIC_SURFACE_2026-05-02.md`.
4. Attach this evidence manifest to the public packet only after IP/publication approval.

## Exact Sentence We Are Allowed To Say

> ErdosAtlas is not yet a theorem-resolution platform at the level of autonomous solvers, but its defensible lane is now visible: a systematic map-check-formalize instrument that turns putative structural correspondences into auditable proof targets, while recording finite passes, controls, and demotions instead of treating analogy as proof.
