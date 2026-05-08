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
# Erdos 242 Salez Filter Certificate Lemmas

Review-only Lean infrastructure for the Erdos-Straus certificate lane.
This file proves the algebraic certificate step behind the seven local
Salez-style filter labels used by the #242 reconstruction scripts.

The lemmas below are sufficiency lemmas for producing cleared-denominator
certificates. They do not prove Erdos-Straus, do not reproduce the `10^18`
Mihnea/Bogdan computation, and do not formalize the optimized upstream sieve.
-/

namespace Erdos242Scout

/--
Order-free cleared-denominator certificate for `4 / n = 1 / x + 1 / y + 1 / z`.

The audit scripts may additionally sort denominators for presentation, but the
algebraic certificate does not need a strict denominator order.
-/
@[category research solved, AMS 11]
def erdosStrausClearedWitness (n x y z : Nat) : Prop :=
  0 < n /\ 0 < x /\ 0 < y /\ 0 < z /\
    4 * x * y * z = n * (y * z + x * z + x * y)

/-- Rosati case 1: `4ABCD = A + B + pC` yields a cleared certificate. -/
@[category research solved, AMS 11]
theorem rosati1_clearedWitness
    {p A B C D : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hRos : 4 * A * B * C * D = A + B + p * C) :
    erdosStrausClearedWitness p (p * B * C * D) (p * A * C * D) (A * B * D) := by
  refine ⟨hp, ?_, ?_, ?_, ?_⟩
  · positivity
  · positivity
  · positivity
  · calc
      4 * (p * B * C * D) * (p * A * C * D) * (A * B * D)
          = p * p * A * B * C * D * D * (A + B + p * C) := by
            rw [← hRos]
            ring
      _ = p * ((p * A * C * D) * (A * B * D)
              + (p * B * C * D) * (A * B * D)
              + (p * B * C * D) * (p * A * C * D)) := by ring

/-- Rosati case 2: `4ABCD = p(A+B) + C` yields a cleared certificate. -/
@[category research solved, AMS 11]
theorem rosati2_clearedWitness
    {p A B C D : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hRos : 4 * A * B * C * D = p * (A + B) + C) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) := by
  refine ⟨hp, ?_, ?_, ?_, ?_⟩
  · positivity
  · positivity
  · positivity
  · calc
      4 * (B * C * D) * (A * C * D) * (p * A * B * D)
          = p * A * B * C * D * D * (p * (A + B) + C) := by
            rw [← hRos]
            ring
      _ = p * ((A * C * D) * (p * A * B * D)
              + (B * C * D) * (p * A * B * D)
              + (B * C * D) * (A * C * D)) := by ring

@[category research solved, AMS 11]
private theorem rosati1_from_AB_E
    {p A B C D E : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hAB : A + B = C * E)
    (hpE : p + E = 4 * A * B * D) :
    erdosStrausClearedWitness p (p * B * C * D) (p * A * C * D) (A * B * D) := by
  apply rosati1_clearedWitness hp hA hB hC hD
  calc
    4 * A * B * C * D = C * (4 * A * B * D) := by ring
    _ = C * (p + E) := by rw [← hpE]
    _ = A + B + p * C := by
      rw [hAB]
      ring

@[category research solved, AMS 11]
private theorem rosati2_from_AB_E
    {p A B C D E : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hAB : A + B = C * E)
    (hpE : p * E + 1 = 4 * A * B * D) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) := by
  apply rosati2_clearedWitness hp hA hB hC hD
  calc
    4 * A * B * C * D = C * (4 * A * B * D) := by ring
    _ = C * (p * E + 1) := by rw [← hpE]
    _ = p * (A + B) + C := by
      rw [hAB]
      ring

/-- `eqmod1a` quotient form used by the local Salez reconstruction. -/
@[category research solved, AMS 11]
theorem salez_eqmod1a_clearedWitness
    {p A B C D : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hQ : B + p * C + A = 4 * A * B * C * D) :
    erdosStrausClearedWitness p (p * B * C * D) (p * A * C * D) (A * B * D) := by
  apply rosati1_clearedWitness hp hA hB hC hD
  rw [← hQ]
  ring

/-- `eqmod1b` quotient form: `A+B = CE` and `p+E = 4ABD`. -/
@[category research solved, AMS 11]
theorem salez_eqmod1b_clearedWitness
    {p A B C D E : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hAB : A + B = C * E)
    (hpE : p + E = 4 * A * B * D) :
    erdosStrausClearedWitness p (p * B * C * D) (p * A * C * D) (A * B * D) :=
  rosati1_from_AB_E hp hA hB hC hD hAB hpE

/-- `eqmod1c` quotient form; the second quotient recovers `A+B = CE`. -/
@[category research solved, AMS 11]
theorem salez_eqmod1c_clearedWitness
    {p A B C D E : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hpE : p + E = 4 * A * B * D)
    (hCq : 4 * B * D * E * C = p + E + 4 * B * B * D) :
    erdosStrausClearedWitness p (p * B * C * D) (p * A * C * D) (A * B * D) := by
  have hEC : E * C = A + B := by
    apply Nat.mul_left_cancel (n := 4 * B * D) (by positivity)
    calc
      (4 * B * D) * (E * C) = 4 * B * D * E * C := by ring
      _ = p + E + 4 * B * B * D := hCq
      _ = 4 * A * B * D + 4 * B * B * D := by rw [hpE]
      _ = (4 * B * D) * (A + B) := by ring
  have hAB : A + B = C * E := by
    rw [← hEC]
    ring
  exact rosati1_from_AB_E hp hA hB hC hD hAB hpE

/-- `eqmod2a` quotient form: `A+B = CE` and `pE+1 = 4ABD`. -/
@[category research solved, AMS 11]
theorem salez_eqmod2a_clearedWitness
    {p A B C D E : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hAB : A + B = C * E)
    (hpE : p * E + 1 = 4 * A * B * D) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) :=
  rosati2_from_AB_E hp hA hB hC hD hAB hpE

/-- `eqmod2b` quotient form: `p+F = 4BCD` and `pB+C = AF`. -/
@[category research solved, AMS 11]
theorem salez_eqmod2b_clearedWitness
    {p A B C D F : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hpF : p + F = 4 * B * C * D)
    (hAF : p * B + C = A * F) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) := by
  apply rosati2_clearedWitness hp hA hB hC hD
  calc
    4 * A * B * C * D = A * (4 * B * C * D) := by ring
    _ = A * (p + F) := by rw [← hpF]
    _ = A * p + A * F := by ring
    _ = A * p + (p * B + C) := by rw [← hAF]
    _ = p * (A + B) + C := by
      ring

/-- `eqmod2c` quotient form; the two quotients recover `pE+1 = 4ABD`. -/
@[category research solved, AMS 11]
theorem salez_eqmod2c_clearedWitness
    {p A B C D E F : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hAB : A + B = C * E)
    (hpF : p + F = 4 * B * D * C)
    (hEF : 4 * B * B * D + 1 = E * F) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) := by
  have hpE : p * E + 1 = 4 * A * B * D := by
    apply Nat.add_right_cancel (m := 4 * B * B * D)
    calc
      (p * E + 1) + 4 * B * B * D = p * E + (4 * B * B * D + 1) := by ring
      _ = p * E + E * F := by rw [hEF]
      _ = (p + F) * E := by ring
      _ = (4 * B * D * C) * E := by rw [hpF]
      _ = 4 * B * D * (C * E) := by ring
      _ = 4 * B * D * (A + B) := by rw [← hAB]
      _ = 4 * A * B * D + 4 * B * B * D := by ring
  exact rosati2_from_AB_E hp hA hB hC hD hAB hpE

/-- `eqmod2d` quotient form: `p+F = 4BCD` and `pB+C = AF`. -/
@[category research solved, AMS 11]
theorem salez_eqmod2d_clearedWitness
    {p A B C D F : Nat}
    (hp : 0 < p) (hA : 0 < A) (hB : 0 < B) (hC : 0 < C) (hD : 0 < D)
    (hpF : p + F = 4 * B * C * D)
    (hAF : p * B + C = A * F) :
    erdosStrausClearedWitness p (B * C * D) (A * C * D) (p * A * B * D) :=
  salez_eqmod2b_clearedWitness hp hA hB hC hD hpF hAF

end Erdos242Scout
