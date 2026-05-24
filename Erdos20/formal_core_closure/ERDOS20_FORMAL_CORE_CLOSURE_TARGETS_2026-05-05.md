# Erdős #20 Formal Core-Closure Targets

**Date:** 2026-05-05  
**Lane:** sunflower core/petal/core-closure formal bookkeeping  
**Claim ceiling:** structural encoding / proof target only. No sunflower theorem progress is claimed here.

## Inputs Inspected

- `erdos-experiments/Erdos20/SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md`
- `Math-Problems/proofs/archive/Sunflower_MDL.lean`
- `erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md`
- `erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md`
- `erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md`

The controlling interpretation is the May 5 packet: the April measurements define an internal Leg-4 observable, not a public lower-bound contribution. The phrase to preserve is **shadow signature, not universal law**.

## Lean-Shaped Vocabulary

The minimum finite-set vocabulary is:

```lean
def Petal (C A : Finset α) : Finset α :=
  A \ C

def HasCore (C : Finset α) (F : Finset (Finset α)) : Prop :=
  ∀ A ∈ F, ∀ B ∈ F, A ≠ B → A ∩ B = C

def DisjointAwayFromCore (C : Finset α) (F : Finset (Finset α)) : Prop :=
  ∀ A ∈ F, ∀ B ∈ F, A ≠ B → Disjoint (Petal C A) (Petal C B)
```

For core-closure bookkeeping, the safest first carrier is a finite state object:

```lean
structure CoreClosureState where
  family : Finset (Finset α)
  core : Finset α
  candidateExtensions : Finset (Finset α)
  safeCandidates : Finset (Finset α)
  safe_subset_candidates : safeCandidates ⊆ candidateExtensions

def CoreClosureCost (S : CoreClosureState (α := α)) : ℕ :=
  S.candidateExtensions.card - S.safeCandidates.card
```

This `CoreClosureCost` is intentionally only a Nat-valued structural proxy: candidate extensions minus safe candidate extensions. It is not the analytic closure information `I_core(C,s,m) = -log₂ p_safe(C,s,m)`. The analytic version should wait until per-core instrumentation exists and the numerator is pre-registered.

## Smallest Compiled Lemma Target

The smallest real lemma target is pure bookkeeping:

```lean
lemma hasCore_disjointAwayFromCore {C : Finset α} {F : Finset (Finset α)}
    (hCore : HasCore C F) :
    DisjointAwayFromCore C F
```

Meaning: if distinct sets in the family all share the same pairwise intersection `C`, then after removing `C`, their petals are disjoint. This formalizes the core/petal split; it does not bound sunflower-free family size, does not prove a sunflower exists, and does not move Erdős #20.

## Scratch Lean Status

Optional scratch file created:

- `Erdos20CoreClosureScratch.lean`

Verified command:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/Math-Problems/proofs/archive
lake env lean /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos20/formal_core_closure/Erdos20CoreClosureScratch.lean
```

This command passed cleanly in the current checkout. Artifact status: compiled scratch theorem, local only, not registered in D1, not scorecarded, not public.

## What Not To Claim

- Do not claim progress on the sunflower conjecture.
- Do not claim the April data improves Erdős-Rado lower bounds.
- Do not claim a Leg-4 PASS.
- Do not call this a compiled portfolio theorem unless it is moved into a maintained Lean project and built through that project.

Safe wording: "We have a compiled local scratch lemma for the core/petal bookkeeping target."

## Next Blocker

The next blocker is not Lean syntax. It is semantic: define the real per-core safe-extension predicate against the enumerator output. Until `safeCandidates` is bound to actual per-core instrumentation, the formal lane should stay at structural encoding and bookkeeping lemmas.
