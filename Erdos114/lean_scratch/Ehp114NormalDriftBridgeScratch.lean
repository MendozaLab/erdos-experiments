import Mathlib

/-!
# EHP114 n=14 Normal-Drift Bridge Scratch

This file records only the finite-sum implication needed after the accepted
Taylor-collar row contract. It does not encode complex analysis, exact
lemniscate length, or a proof of Erdős #114.
-/

namespace Erdos
namespace EHP114
namespace NormalDrift

open scoped BigOperators

/-- One row in a finite bridge certificate. -/
structure BridgeRow where
  normalDriftError : ℚ
  exactLengthError : ℚ
  regularityResolved : Prop

/--
Finite normal-drift contract: every row is regularity-resolved, every exact
length error contribution is controlled by the normal-drift contribution, and
the finite normal-drift sum is within the available budget.
-/
def FiniteNormalDriftContract {ι : Type} [Fintype ι]
    (rows : ι → BridgeRow) (budget : ℚ) : Prop :=
  (∀ i, (rows i).regularityResolved ∧
    (rows i).exactLengthError ≤ (rows i).normalDriftError) ∧
  (∑ i, (rows i).normalDriftError) ≤ budget

/--
The finite certificate implication: once the analytic row-wise comparison is
available, the accepted finite normal-drift sum carries the exact-error budget.
-/
theorem finite_normal_drift_contract_implies_exact_error_budget
    {ι : Type} [Fintype ι] (rows : ι → BridgeRow) (budget : ℚ)
    (h : FiniteNormalDriftContract rows budget) :
    (∑ i, (rows i).exactLengthError) ≤ budget := by
  rcases h with ⟨hrows, hbudget⟩
  exact le_trans (Finset.sum_le_sum (fun i _hi => (hrows i).2)) hbudget

/-- The current accepted row count from the immutable row-contract artifact. -/
def acceptedRowCount : ℕ := 15793

/-- The row-contract artifact hash, kept as metadata rather than proof data. -/
def acceptedRowContractSha256 : String :=
  "8417f5903d3696fe300961359e65e4771c1c94269d067e4dd81db7341c422e47"

end NormalDrift
end EHP114
end Erdos

