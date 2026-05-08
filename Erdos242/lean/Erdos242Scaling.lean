/-
Copyright 2026 The Formal Conjectures Authors.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-/

import FormalConjectures.Util.ProblemImports

/-!
# Erdos 242 Scaling Infrastructure

Review-only Lean infrastructure for the Erdos-Straus residue-coverage lane.
This file proves the cleared-denominator scaling lemma used by the local
coverage operator: once a certificate exists for `n`, multiplying all
denominators by `k` gives a certificate for `k * n`.
-/

namespace Erdos242Scout

/-- Cleared-denominator form of `4 / n = 1 / x + 1 / y + 1 / z`. -/
@[category research solved, AMS 11]
def erdosStrausClearedSolves (n x y z : Nat) : Prop :=
  2 < n /\ 0 < x /\ x < y /\ y < z /\
    4 * x * y * z = n * (y * z + x * z + x * y)

/--
Scaling lemma for the cleared Erdos-Straus equation.

This does not prove Erdos #242. It proves the reusable step behind the
divisor-descent coverage operator: a certificate for `n` induces one for
`k * n`.
-/
@[category research solved, AMS 11]
theorem erdosStrausClearedSolves_mul
    {n x y z k : Nat}
    (hk : 0 < k)
    (h : erdosStrausClearedSolves n x y z) :
    erdosStrausClearedSolves (k * n) (k * x) (k * y) (k * z) := by
  rcases h with ⟨hn, hx, hxy, hyz, hEq⟩
  refine ⟨?_, ?_, ?_, ?_, ?_⟩
  · have hkn : n <= k * n := by
      nth_rewrite 1 [← one_mul n]
      exact Nat.mul_le_mul_right n (Nat.succ_le_of_lt hk)
    exact lt_of_lt_of_le hn hkn
  · exact Nat.mul_pos hk hx
  · exact Nat.mul_lt_mul_of_pos_left hxy hk
  · exact Nat.mul_lt_mul_of_pos_left hyz hk
  · calc
      4 * (k * x) * (k * y) * (k * z)
          = k ^ 3 * (4 * x * y * z) := by ring
      _ = k ^ 3 * (n * (y * z + x * z + x * y)) := by rw [hEq]
      _ = (k * n) * ((k * y) * (k * z) + (k * x) * (k * z) + (k * x) * (k * y)) := by ring

end Erdos242Scout
