import Mathlib

/-!
# Face/Field Response Lemmas

This file formalizes the first proof target from the Collider face-field packet.

It is intentionally abstract.  There is no Sidon theorem here and no claim about
Erdos #30.  The point is the small finite-set lemma that justifies the language
"the extremal face is the object": if two observables expose different
zero-temperature witnesses on the same finite family, then the family cannot be
treated as a unique optimizer selected by enumeration order.
-/

namespace Erdos
namespace Collider

variable {α β γ : Type*}

/-- `IsFieldMinOn F φ x` says that `x` is a minimizer of observable `φ` on
the finite family `F`.  The word "field" is only paper vocabulary; formally this
is just minimization of an ordered function on a `Finset`. -/
def IsFieldMinOn [LinearOrder β] (F : Finset α) (φ : α → β) (x : α) : Prop :=
  x ∈ F ∧ ∀ y ∈ F, φ x ≤ φ y

/-- `IsFieldMinimizerSet F S φ` says that `S` is exactly the selected
minimizer surface of observable `φ` on the finite family `F`. This is the
surface-valued version of `IsFieldMinOn`, used once a field selects a tied face
instead of a unique point. -/
def IsFieldMinimizerSet [LinearOrder β] (F S : Finset α) (φ : α → β) : Prop :=
  S ⊆ F ∧ ∀ x : α, x ∈ S ↔ IsFieldMinOn F φ x

/-- Certify a whole minimizer surface by proving every point on the surface
minimizes and every off-surface point in the finite face is strictly beaten by
some surface point. -/
theorem isFieldMinimizerSet_of_le_on_and_lt_off [LinearOrder β]
    {F S : Finset α} {φ : α → β}
    (hS : S ⊆ F)
    (hle : ∀ x ∈ S, ∀ y ∈ F, φ x ≤ φ y)
    (hoff : ∀ y ∈ F, y ∉ S → ∃ x ∈ S, φ x < φ y) :
    IsFieldMinimizerSet F S φ := by
  constructor
  · exact hS
  · intro x
    constructor
    · intro hx
      exact ⟨hS hx, hle x hx⟩
    · intro hxmin
      by_contra hxnot
      rcases hoff x hxmin.1 hxnot with ⟨y, hyS, hylt⟩
      exact (not_lt_of_ge (hxmin.2 y (hS hyS))) hylt

/-- Minimizer transfer across two observables that agree on the finite face.
This is the reusable core behind replacing a packet-probe field by a scalar
field after proving pointwise equality on the exported face. -/
theorem isFieldMinOn_of_eq_on [LinearOrder β]
    {F : Finset α} {φ ψ : α → β} {x : α}
    (hmin : IsFieldMinOn F φ x)
    (heq : ∀ y ∈ F, ψ y = φ y) :
    IsFieldMinOn F ψ x := by
  constructor
  · exact hmin.1
  · intro y hy
    rw [heq x hmin.1, heq y hy]
    exact hmin.2 y hy

/-- Weighted-joint transfer: if the prefix component is replaced by an equal
scalar observable on the finite face, then a minimizer of the probe joint field
is also a minimizer of the scalar joint field. -/
theorem isFieldMinOn_weightedJoint_of_prefix_eq_on
    {ι : Type*} {F : Finset ι} {x : ι}
    {prefixScalar prefixProbe mass : ι → ℝ} {weight : ℝ}
    (hmin :
      IsFieldMinOn F
        (fun i => prefixProbe i * weight + mass i) x)
    (hprefix : ∀ i ∈ F, prefixScalar i = prefixProbe i) :
    IsFieldMinOn F
      (fun i => prefixScalar i * weight + mass i) x := by
  exact isFieldMinOn_of_eq_on hmin (by
    intro i hi
    rw [hprefix i hi])

/-- Minimizer-surface transfer across two observables that agree on the finite
face. -/
theorem isFieldMinimizerSet_of_eq_on [LinearOrder β]
    {F S : Finset α} {φ ψ : α → β}
    (hset : IsFieldMinimizerSet F S φ)
    (heq : ∀ y ∈ F, ψ y = φ y) :
    IsFieldMinimizerSet F S ψ := by
  constructor
  · exact hset.1
  · intro x
    constructor
    · intro hx
      exact isFieldMinOn_of_eq_on ((hset.2 x).mp hx) heq
    · intro hxmin
      exact (hset.2 x).mpr
        (isFieldMinOn_of_eq_on hxmin (by
          intro y hy
          exact (heq y hy).symm))

/-- Weighted-joint transfer for selected surfaces. -/
theorem isFieldMinimizerSet_weightedJoint_of_prefix_eq_on
    {ι : Type*} {F S : Finset ι}
    {prefixScalar prefixProbe mass : ι → ℝ} {weight : ℝ}
    (hset :
      IsFieldMinimizerSet F S
        (fun i => prefixProbe i * weight + mass i))
    (hprefix : ∀ i ∈ F, prefixScalar i = prefixProbe i) :
    IsFieldMinimizerSet F S
      (fun i => prefixScalar i * weight + mass i) := by
  exact isFieldMinimizerSet_of_eq_on hset (by
    intro i hi
    rw [hprefix i hi])

/-- `SurfaceTransition S I N` says that a selected surface `S` splits into an
inherited part `I` and a new part `N`. This is deliberately abstract: the
meaning of inherited is supplied by the finite certificate, e.g. previous-row
or shifted-previous-row membership. -/
def SurfaceTransition [DecidableEq α] (S inherited novel : Finset α) : Prop :=
  inherited ⊆ S ∧ novel ⊆ S ∧ Disjoint inherited novel ∧ S = inherited ∪ novel

/-- A selected-surface transition decomposes cardinality additively. -/
theorem surfaceTransition_card_eq [DecidableEq α]
    {S inherited novel : Finset α}
    (h : SurfaceTransition S inherited novel) :
    S.card = inherited.card + novel.card := by
  rcases h with ⟨_, _, hdisjoint, hunion⟩
  rw [hunion]
  exact Finset.card_union_of_disjoint hdisjoint

/-- Two observables split the face if their selected minimizing witnesses can be
chosen differently. -/
def FieldSplit [LinearOrder β] [LinearOrder γ] (F : Finset α)
    (φ : α → β) (ψ : α → γ) : Prop :=
  ∃ x y : α, IsFieldMinOn F φ x ∧ IsFieldMinOn F ψ y ∧ x ≠ y

/-- The abstract exposed-face lemma: a field split certifies that the finite
family has at least two points. -/
theorem fieldSplit_card_two_le [DecidableEq α] [LinearOrder β] [LinearOrder γ]
    {F : Finset α} {φ : α → β} {ψ : α → γ}
    (h : FieldSplit F φ ψ) :
    2 ≤ F.card := by
  rcases h with ⟨x, y, hx, hy, hxy⟩
  have hxmem : x ∈ F := hx.1
  have hymem : y ∈ F := hy.1
  have hsubset : ({x, y} : Finset α) ⊆ F := by
    intro z hz
    simp only [Finset.mem_insert, Finset.mem_singleton] at hz
    rcases hz with rfl | rfl
    · exact hxmem
    · exact hymem
  have hcard_pair : ({x, y} : Finset α).card = 2 := by
    rw [Finset.card_pair]
    exact hxy
  calc
    2 = ({x, y} : Finset α).card := hcard_pair.symm
    _ ≤ F.card := Finset.card_le_card hsubset

/-- A face with a field split cannot be a singleton face. -/
theorem fieldSplit_not_card_le_one [DecidableEq α] [LinearOrder β] [LinearOrder γ]
    {F : Finset α} {φ : α → β} {ψ : α → γ}
    (h : FieldSplit F φ ψ) :
    ¬ F.card ≤ 1 := by
  intro hle
  have htwo : 2 ≤ F.card := fieldSplit_card_two_le h
  exact Nat.not_succ_le_self 1 (le_trans htwo hle)

end Collider
end Erdos
