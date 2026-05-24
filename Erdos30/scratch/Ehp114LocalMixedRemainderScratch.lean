/-
  Ehp114LocalMixedRemainderScratch.lean
  =====================================

  Scratch theorem-target scaffold for Erdős #114, n = 14.

  This file does not prove the analytic mixed-remainder theorem. It defines the
  local objects and checks the closure algebra:

    radial Puiseux deficit + positive shape cone - mixed remainder > 0.

  The immediate live missing theorem is now `EpsScaledDeficit14`. The older
  `MixedRemainderAbsorption14` target remains as a stronger later route.

  Build target:

    lake build Ehp114LocalMixedRemainderScratch

  Claim ceiling:
  - build PASS means the theorem target and algebraic splice are well typed;
  - it is not a proof of EHP #114;
  - it is not a public COMPILED theorem claim for the scorecard.
-/

import Mathlib

noncomputable section

namespace Erdos114
namespace LocalN14

open Real

/-! ## Local quotient objects -/

/-- The n = 14 shape quotient has observed rank `2n - 3 = 25`.

This concrete model is a scratch representation of quotient coordinates, not
the final analytic construction of the quotient by scale/rotation symmetries. -/
abbrev ShapeQuotient14 := Fin 25 → ℝ

/-- Squared Euclidean norm in the current quotient-coordinate model. -/
def quotientNormSq (s : ShapeQuotient14) : ℝ :=
  ∑ i : Fin 25, (s i) ^ 2

/-- Euclidean norm in the current quotient-coordinate model. -/
def quotientNorm (s : ShapeQuotient14) : ℝ :=
  Real.sqrt (quotientNormSq s)

theorem quotientNormSq_nonneg (s : ShapeQuotient14) :
    0 ≤ quotientNormSq s := by
  unfold quotientNormSq
  exact Finset.sum_nonneg (fun i _ => sq_nonneg (s i))

/-! ## Constants from the artifact stack -/

/-- Certified radial Puiseux constant used by the n14 interval packet. -/
def radialPuiseuxConstant14 : ℝ := 24

/-- Half-radial reserve reserved after mixed-remainder absorption. -/
def radialReserve14 : ℝ := 12

/-- Working shape-cone lower bound below the interval Gershgorin certificate. -/
def shapeLambda14 : ℝ := 100000

/-- Half-shape reserve reserved after mixed-remainder absorption. -/
def shapeReserve14 : ℝ := 50000

/-- Current local scout cap. This is not a universal constant. -/
def localConeCap14 : ℝ := (8 : ℝ) / 1000

/-- Current epsilon-scaled cone constant supported by finite n14 axis and
spectral-direction searches. This is a theorem-target parameter, not a
universal constant. -/
def eta0_14 : ℝ := (14 : ℝ) / 1000

/-! ## Abstract analytic objects -/

/-- Analytic admissible boundary radius still to be proved. -/
axiom eta14Boundary : ℝ → ℝ

/-- Radial Puiseux deficit `L14(1) - L14(1 - eps)`. -/
axiom radialDeficit14 : ℝ → ℝ

/-- Shape-cone quadratic energy after quotienting the symmetry modes. -/
axiom shapeQuadratic14 : ShapeQuotient14 → ℝ

/-- Mixed radial/shape remainder that must be absorbed. -/
axiom mixedRemainder14 : ℝ → ShapeQuotient14 → ℝ

/-- Unit-disk root admissibility for the radial plus shape perturbation. -/
axiom RootsInClosedUnitDisk14 : ℝ → ShapeQuotient14 → Prop

/-- The local cone where the finite mixed scout passed. -/
def LocalCone14 (eps : ℝ) (s : ShapeQuotient14) : Prop :=
  RootsInClosedUnitDisk14 eps s ∧
    quotientNorm s ≤ min (eta14Boundary eps) localConeCap14

/-- The current epsilon-scaled cone target after the admissible spectral and
axis/spectral-direction scans. -/
def EpsilonScaledCone14 (eps : ℝ) (s : ShapeQuotient14) : Prop :=
  RootsInClosedUnitDisk14 eps s ∧
    quotientNorm s ≤ eta0_14 * Real.rpow eps ((1 : ℝ) / 28)

/-- Local deficit after decomposing into radial, shape, and mixed terms. -/
def totalDeficit14 (eps : ℝ) (s : ShapeQuotient14) : ℝ :=
  radialDeficit14 eps + shapeQuadratic14 s - mixedRemainder14 eps s

/-! ## Component propositions -/

/-- Radial Puiseux lower bound. Artifact-backed, but not proved natively here. -/
def RadialCertificate14 : Prop :=
  ∀ eps : ℝ, 0 < eps → eps ≤ (1 : ℝ) / 10 →
    radialPuiseuxConstant14 * Real.rpow eps ((1 : ℝ) / 14) ≤
      radialDeficit14 eps

/-- Shape-cone positivity. Artifact-backed by the interval matrix certificate. -/
def ShapeConeCertificate14 : Prop :=
  ∀ s : ShapeQuotient14,
    shapeLambda14 * quotientNormSq s ≤ shapeQuadratic14 s

/-- The single live theorem target. This is not proved by the finite scout. -/
def MixedRemainderAbsorption14 : Prop :=
  ∀ eps : ℝ, ∀ s : ShapeQuotient14,
    0 < eps → eps ≤ (1 : ℝ) / 10 → LocalCone14 eps s →
      mixedRemainder14 eps s ≤
        radialReserve14 * Real.rpow eps ((1 : ℝ) / 14) +
          shapeReserve14 * quotientNormSq s

/-- Current scalar theorem target. This bypasses the too-strong transported
positive shape-cone claim and asks directly for the radial Puiseux reserve floor
on the epsilon-scaled admissible cone. -/
def EpsScaledDeficit14 : Prop :=
  ∀ eps : ℝ, ∀ s : ShapeQuotient14,
    0 < eps → eps ≤ (1 : ℝ) / 10 → EpsilonScaledCone14 eps s →
      radialReserve14 * Real.rpow eps ((1 : ℝ) / 14) ≤
        totalDeficit14 eps s

/-! ## Closure algebra -/

/-- If the three component estimates are supplied, the local deficit retains
half of the radial budget and half of the shape budget.

This theorem is just algebra. The analytic work is isolated in
`MixedRemainderAbsorption14`. -/
theorem local_deficit_reserve_from_components
    (hRadial : RadialCertificate14)
    (hShape : ShapeConeCertificate14)
    (hMixed : MixedRemainderAbsorption14)
    (eps : ℝ) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps ≤ (1 : ℝ) / 10)
    (hlocal : LocalCone14 eps s) :
    radialReserve14 * Real.rpow eps ((1 : ℝ) / 14) +
        shapeReserve14 * quotientNormSq s ≤ totalDeficit14 eps s := by
  unfold RadialCertificate14 at hRadial
  unfold ShapeConeCertificate14 at hShape
  unfold MixedRemainderAbsorption14 at hMixed
  unfold totalDeficit14
  unfold radialPuiseuxConstant14 radialReserve14 shapeLambda14 shapeReserve14 at *
  have hR := hRadial eps hpos hsmall
  have hQ := hShape s
  have hM := hMixed eps s hpos hsmall hlocal
  linarith

/-- The stronger mixed-remainder route implies the weaker scalar reserve on
the old local cone. This is only algebra; it does not prove the analytic mixed
remainder theorem. -/
theorem scalar_reserve_from_components
    (hRadial : RadialCertificate14)
    (hShape : ShapeConeCertificate14)
    (hMixed : MixedRemainderAbsorption14)
    (eps : ℝ) (s : ShapeQuotient14)
    (hpos : 0 < eps)
    (hsmall : eps ≤ (1 : ℝ) / 10)
    (hlocal : LocalCone14 eps s) :
    radialReserve14 * Real.rpow eps ((1 : ℝ) / 14) ≤
      totalDeficit14 eps s := by
  have h := local_deficit_reserve_from_components hRadial hShape hMixed eps s hpos hsmall hlocal
  have hShapeReserveNonneg : 0 ≤ shapeReserve14 * quotientNormSq s := by
    unfold shapeReserve14
    nlinarith [quotientNormSq_nonneg s]
  linarith

/-- The current theorem target restated as a named declaration shape for
downstream agents. It is intentionally an axiom in this scratch file: proving or
interval-certifying this proposition is the next real task. -/
axiom ehp114_n14_eps_scaled_cone_deficit :
    EpsScaledDeficit14

/-- The stronger theorem target restated as a named declaration shape for
downstream agents. It is intentionally an axiom in this scratch file: proving or
interval-certifying this proposition is the next real task. -/
axiom ehp114_n14_local_mixed_remainder_absorption :
    MixedRemainderAbsorption14

end LocalN14
end Erdos114
