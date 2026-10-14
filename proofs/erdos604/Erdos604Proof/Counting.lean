import Erdos604Proof.Definitions

/-! Finite counting adapter. No theorem from the external proof closure is imported. -/

noncomputable section
open scoped BigOperators

namespace PinnedDistanceAdapter

/-- Partition a finite set into fibers. Light fibers contribute at most `T` each;
the remaining points are counted individually by the rich-point term. -/
theorem card_le_image_mul_threshold_add_rich
    {α β : Type*} [DecidableEq β] (P : Finset α) (f : α → β)
    (T : ℝ) (hT : 0 ≤ T) :
    (P.card : ℝ) ≤ ((P.image f).card : ℝ) * T +
      ((P.filter (fun y => T ≤ ((P.filter (fun z => f z = f y)).card : ℝ))).card : ℝ) := by
  classical
  let R := P.filter (fun y => T ≤ ((P.filter (fun z => f z = f y)).card : ℝ))
  have hP : (P.card : ℝ) =
      ∑ b ∈ P.image f, ((P.filter (fun a => f a = b)).card : ℝ) := by
    exact_mod_cast Finset.card_eq_sum_card_image f P
  have hR : (R.card : ℝ) =
      ∑ b ∈ P.image f, ((R.filter (fun a => f a = b)).card : ℝ) := by
    exact_mod_cast (Finset.card_eq_sum_card_fiberwise
      (s := R) (t := P.image f) (f := f)
      (fun a ha => Finset.mem_image_of_mem f (Finset.mem_filter.mp ha).1))
  have hb (b : β) : ((P.filter (fun a => f a = b)).card : ℝ) ≤
      T + ((R.filter (fun a => f a = b)).card : ℝ) := by
    by_cases h : T ≤ ((P.filter (fun a => f a = b)).card : ℝ)
    · have heq : R.filter (fun a => f a = b) = P.filter (fun a => f a = b) := by
        ext a
        simp only [R, Finset.mem_filter]
        constructor
        · exact fun ha => ⟨ha.1.1, ha.2⟩
        · intro ha
          exact ⟨⟨ha.1, by simpa only [ha.2] using h⟩, ha.2⟩
      rw [heq]
      linarith
    · have hnonneg : 0 ≤ ((R.filter (fun a => f a = b)).card : ℝ) := by positivity
      linarith
  calc
    (P.card : ℝ) = ∑ b ∈ P.image f, ((P.filter (fun a => f a = b)).card : ℝ) := hP
    _ ≤ ∑ b ∈ P.image f, (T + ((R.filter (fun a => f a = b)).card : ℝ)) :=
      Finset.sum_le_sum (fun b _ => hb b)
    _ = ((P.image f).card : ℝ) * T + (R.card : ℝ) := by
      rw [Finset.sum_add_distrib, ← hR]
      simp

open OAI.WeakPinned

/-- Distances from `x` to the other points of `P`; zero is excluded even when `x ∉ P`. -/
def pinnedDistances (P : Finset Plane) (x : Plane) : Finset ℝ :=
  (P.erase x).image (fun y => dist y x)

/-- Pins with at most `D` distinct distances to other points, with a real cutoff. -/
def badPins (P : Finset Plane) (D : ℝ) : Finset Plane :=
  P.filter (fun x => ((pinnedDistances P x).card : ℝ) ≤ D)

def richNeighbors (P : Finset Plane) (s : ℝ) (x : Plane) : Finset Plane :=
  (P.erase x).filter (fun y => (P.card : ℝ) ^ s ≤ (k P x y : ℝ))

theorem richNeighbors_eq_filter (P : Finset Plane) (s : ℝ) (x : Plane) :
    richNeighbors P s x =
      P.filter (fun y => x ≠ y ∧ (P.card : ℝ) ^ s ≤ (k P x y : ℝ)) := by
  classical
  ext y
  simp [richNeighbors, ne_comm, and_left_comm, and_assoc]

theorem richPairs_card_eq_sum (P : Finset Plane) (s : ℝ) :
    (richPairs P s).card = ∑ x ∈ P, (richNeighbors P s x).card := by
  classical
  simp only [richPairs, richNeighbors_eq_filter, Finset.card_eq_sum_ones,
    Finset.sum_filter, Finset.sum_product]

/-- The finite adapter to the exact external ordered-pair definition.
No lower cardinality hypothesis is needed, and `D` and `s` are arbitrary reals. -/
theorem badPins_mul_le_richPairs (P : Finset Plane) (D s : ℝ) :
    ((badPins P D).card : ℝ) * ((P.card : ℝ) - 1 - D * (P.card : ℝ) ^ s) ≤
      ((richPairs P s).card : ℝ) := by
  classical
  have hpin (x : Plane) (hx : x ∈ badPins P D) :
      (P.card : ℝ) - 1 - D * (P.card : ℝ) ^ s ≤
        ((richNeighbors P s x).card : ℝ) := by
    obtain ⟨hxP, hxD⟩ := Finset.mem_filter.mp hx
    have hT : 0 ≤ (P.card : ℝ) ^ s := Real.rpow_nonneg (Nat.cast_nonneg _) _
    have h := card_le_image_mul_threshold_add_rich
      (P.erase x) (fun y => dist y x) ((P.card : ℝ) ^ s) hT
    have herase : ((P.erase x).card : ℝ) = (P.card : ℝ) - 1 :=
      Finset.cast_card_erase_of_mem hxP
    change ((pinnedDistances P x).card : ℝ) ≤ D at hxD
    change ((P.erase x).card : ℝ) ≤ ((pinnedDistances P x).card : ℝ) *
      (P.card : ℝ) ^ s + ((richNeighbors P s x).card : ℝ) at h
    rw [herase] at h
    have hmul := mul_le_mul_of_nonneg_right hxD hT
    linarith
  have hsubset : badPins P D ⊆ P := Finset.filter_subset _ _
  calc
    ((badPins P D).card : ℝ) * ((P.card : ℝ) - 1 - D * (P.card : ℝ) ^ s) =
        ∑ x ∈ badPins P D, ((P.card : ℝ) - 1 - D * (P.card : ℝ) ^ s) := by
      simp only [Finset.sum_const, nsmul_eq_mul]
    _ ≤ ∑ x ∈ badPins P D, ((richNeighbors P s x).card : ℝ) :=
      Finset.sum_le_sum hpin
    _ ≤ ∑ x ∈ P, ((richNeighbors P s x).card : ℝ) := by
      apply Finset.sum_le_sum_of_subset_of_nonneg hsubset
      intros
      positivity
    _ = ((richPairs P s).card : ℝ) := by
      exact_mod_cast (richPairs_card_eq_sum P s).symm

/-- A useful normalization when the light-fiber contribution is at most half
the `n - 1` neighbors. The denominator matches the external `pairFraction`. -/
theorem badPins_fraction_le_two_pairFraction (P : Finset Plane) (D s : ℝ)
    (hn : 2 ≤ P.card)
    (hsmall : D * (P.card : ℝ) ^ s ≤ ((P.card : ℝ) - 1) / 2) :
    ((badPins P D).card : ℝ) / (P.card : ℝ) ≤ 2 * pairFraction P s := by
  have hn' : (2 : ℝ) ≤ (P.card : ℝ) := by exact_mod_cast hn
  have hnpos : (0 : ℝ) < (P.card : ℝ) := by linarith
  have hmpos : (0 : ℝ) < (P.card : ℝ) - 1 := by linarith
  have hcoeff : ((P.card : ℝ) - 1) / 2 ≤
      (P.card : ℝ) - 1 - D * (P.card : ℝ) ^ s := by linarith
  have hBnonneg : (0 : ℝ) ≤ ((badPins P D).card : ℝ) := Nat.cast_nonneg _
  have hhalf := (mul_le_mul_of_nonneg_left hcoeff hBnonneg).trans
    (badPins_mul_le_richPairs P D s)
  have hfull : ((badPins P D).card : ℝ) * ((P.card : ℝ) - 1) ≤
      2 * ((richPairs P s).card : ℝ) := by linarith
  calc
    ((badPins P D).card : ℝ) / (P.card : ℝ) =
        (((badPins P D).card : ℝ) * ((P.card : ℝ) - 1)) /
          ((P.card : ℝ) * ((P.card : ℝ) - 1)) := by
      field_simp
    _ ≤ (2 * ((richPairs P s).card : ℝ)) /
        ((P.card : ℝ) * ((P.card : ℝ) - 1)) :=
      div_le_div_of_nonneg_right hfull (le_of_lt (mul_pos hnpos hmpos))
    _ = 2 * pairFraction P s := by unfold pairFraction; ring

end PinnedDistanceAdapter
