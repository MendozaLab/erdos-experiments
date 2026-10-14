/-
Copyright 2025-2026 The Formal Conjectures Authors.
Released under Apache 2.0 license as described in the file LICENSE.
Definitions copied from FormalConjecturesForMathlib/Geometry/Metric.lean
and FormalConjectures/ErdosProblems/604.lean at
57f2e0cf2b374563dcb42bdd9adec97aa7e27517.
Modified for this package: placed in Erdos604Proof.Target and expanded
the Euclidean plane notation. The definition bodies are preserved.
-/
import Erdos604Proof.Definitions

open scoped Finset

namespace Erdos604Proof.Target

variable {X : Type*} [MetricSpace X]

open Classical in
noncomputable def distinctDistancesFrom (points : Finset X) (pt : X) : ℕ :=
  #((points.erase pt).image fun x => dist x pt)


noncomputable def pinnedDistinctDistances (n : ℕ) : ℕ :=
  sInf {m : ℕ | ∃ A : Finset OAI.WeakPinned.Plane,
    #A = n ∧ A.sup (distinctDistancesFrom A) = m}

end Erdos604Proof.Target
