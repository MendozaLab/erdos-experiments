import Mathlib
import Erdos30_Sidon_Defs

open Finset Nat Erdos.Sidon

/--
The number of distinct ordered pairwise sums from a Sidon set of size k
is at least k*(k+1)/2. These sums lie in {0, ..., 2*max(A)}.
Therefore k*(k+1)/2 ≤ 2*max(A) + 1.

We prove the slightly weaker: k*(k+1)/2 ≤ 2*max(A) + 1, which suffices
for the upper bound theorem.

PROVIDED SOLUTION
Step 1: Define the set of ordered pairs {(a,b) ∈ A×A : a ≤ b}. Its cardinality is
  |A|*(|A|+1)/2 because there are C(|A|,2) pairs with a<b plus |A| diagonal pairs.

Step 2: The map (a,b) ↦ a+b is INJECTIVE on this set by the Sidon property.
  If a₁+b₁ = a₂+b₂ with a₁≤b₁ and a₂≤b₂, then a₁=a₂ and b₁=b₂.

Step 3: All sums a+b lie in {0, ..., 2*sup(A)} since a,b ≤ sup(A).
  So the IMAGE of the injection has cardinality ≤ 2*sup(A) + 1.

Step 4: By injectivity, |pairs| ≤ |image| ≤ 2*sup(A) + 1.

Key Lean tactics:
- Use Finset.card_image_of_injOn for injectivity
- Use Finset.card_le_card with subset of Finset.range (2 * A.sup id + 1)
- The Sidon property gives the injOn hypothesis directly
- For the pair count: use Finset.card_filter with the (a ≤ b) predicate on A ×ˢ A
-/
theorem sidon_sum_count (A : Finset ℕ) (hS : IsSidonSet A) :
    A.card * (A.card + 1) / 2 ≤ 2 * (A.sup id) + 1 := by
  -- Proof backfilled from Sidon_SumCount_Fix_Solved.lean
  have h_card : (Finset.image (fun (p : ℕ × ℕ) => p.fst + p.snd) (Finset.filter (fun (p : ℕ × ℕ) => p.fst ≤ p.snd) (A ×ˢ A))).card ≥ (A.card * (A.card + 1)) / 2 := by
    rw [ Finset.card_image_of_injOn ];
    · have h_pairs_eq : ((Finset.filter (fun p => p.1 ≤ p.2) (A ×ˢ A)).card = (Finset.filter (fun p => p.1 < p.2) (A ×ˢ A)).card + (Finset.image (fun a => (a, a)) A).card) := by
        rw [ show { p ∈ A ×ˢ A | p.1 ≤ p.2 } = { p ∈ A ×ˢ A | p.1 < p.2 } ∪ Finset.image ( fun a => ( a, a ) ) A from ?_, Finset.card_union_of_disjoint ];
        · rw [ Finset.disjoint_right ] ; aesop;
        · grind +ring;
      have h_pairs_lt : ((Finset.filter (fun p => p.1 < p.2) (A ×ˢ A)).card = A.card * (A.card - 1) / 2) := by
        have h_pairs_lt : ((Finset.filter (fun p => p.1 < p.2) (A ×ˢ A)).card = Finset.card (Finset.powersetCard 2 A)) := by
          refine' Finset.card_bij ( fun p hp => { p.1, p.2 } ) _ _ _;
          · grind;
          · simp +contextual [ Finset.Subset.antisymm_iff, Finset.subset_iff ];
            intros; omega;
          · intro b hb; rw [ Finset.mem_powersetCard ] at hb; rcases Finset.card_eq_two.mp hb.2 with ⟨ x, y, hxy, rfl ⟩ ; cases lt_trichotomy x y <;> aesop;
        rw [ h_pairs_lt, Finset.card_powersetCard, Nat.choose_two_right ];
      cases k : Finset.card A <;> simp_all +decide [ Nat.mul_succ, Finset.card_image_of_injective, Function.Injective ] ; omega;
    · intro p hp q hq; specialize hS p.1 ( by aesop ) p.2 ( by aesop ) q.1 ( by aesop ) q.2 ( by aesop ) ; aesop;
  refine le_trans h_card ?_;
  refine' le_trans ( Finset.card_le_card _ ) _;
  exact Finset.range ( 2 * A.sup id + 1 );
  · exact Finset.image_subset_iff.mpr fun p hp => Finset.mem_range.mpr ( by linarith [ show p.1 ≤ A.sup id from Finset.le_sup ( f := id ) ( Finset.mem_filter.mp hp |>.1 |> Finset.mem_product.mp |>.1 ), show p.2 ≤ A.sup id from Finset.le_sup ( f := id ) ( Finset.mem_filter.mp hp |>.1 |> Finset.mem_product.mp |>.2 ) ] );
  · norm_num

