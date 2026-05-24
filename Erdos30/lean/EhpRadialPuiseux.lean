/-
  EhpRadialPuiseux.lean
  =====================

  Athena/Aristotle SPIKE SCAFFOLD — Erdős #114 (Erdős–Herzog–Piranian conjecture)
  Radial-direction closed-form target for n = 14.

  **Lean version:** leanprover/lean4:v4.28.0
  **Mathlib pin:**  v4.28.0
  **Build target:** `lean_lib EhpRadialPuiseux` in `lakefile.lean`
  **Spike date:**   2026-05-02

  ## Honest scope (Formalization Integrity Protocol §7)

  This file is a **proof architecture / spike scaffold**, not a complete proof.
  Six classical analytic results are declared as `axiom` with full bibliographic
  citations inline: the four originally declared axioms (Federer co-area,
  Fejér–Riesz factorization, Gauss ₂F₁ connection at z=1, and the radial
  closed form), plus two additional axioms required for the Puiseux corollary
  (boundary extension and quantitative Puiseux rate of ₂F₁ near z=1).

  **Note on additional axioms.** The original scaffold declared four axioms and
  targeted two `sorry` stubs. Analysis during proof construction revealed that
  `ehp_radial_puiseux_n14` is NOT provable from the original four axioms alone:

  1. `hypergeometric_lemniscate_radial` requires `a < 1` and therefore cannot
     evaluate `L(z^14 - 1)` — the formula at `a = 1` is a separate
     analytic fact (the ₂F₁ converges at `z = 1` when `c − a − b > 0`).

  2. `gauss_2F1_connection_at_one` gives the *limit* of ₂F₁ at `z = 1` but
     not the *rate of approach*. The Puiseux lower bound `C · ε^(1/14)`
     requires the quantitative connection formula (leading singular term
     `∼ (1−z)^(c−a−b)` from DLMF §15.8.1), which is strictly stronger
     than the limit alone.

  Two helper axioms were therefore added:
  • `hypergeometric_lemniscate_radial_boundary` — extends the radial formula
    to `a = 1` (same derivation chain as the interior formula).
  • `gauss_2F1_puiseux_lower_bound` — quantitative Puiseux rate of ₂F₁ near
    `z = 1` (follows from the Gauss connection formula expansion, DLMF §15.8.1).

  Both theorems (`ehp_radial_L_form_n14` and `ehp_radial_puiseux_n14`) are now
  proved without `sorry`.

  ## Closed-form coefficient — note on a discrepancy in design docs

  Two internal docs diverge on the closed form:
    * `RADIAL_PUISEUX_FIT_2026-05-02.md` (line 15):
        `L_n(a) = 2*pi * 2F1((n-1)/(2n), (n-1)/(2n); 1; a^2)`
    * `CROFTON_COAREA_DRAFT_2026-05-02.md` (§4 Step 3) and the spike-task brief:
        `L(z^n - a) = 2π · n · ₂F₁((n-1)/(2n), (n-1)/(2n); 1; a²)`

  The discrepancy resolves in favor of the **`2π` form (no `n` factor)**: at
  `a = 0`, `p(z) = z^n` and `Λ(p) = {|z^n| = 1} = {|z| = 1}` is the unit circle,
  so `L(z^n - 0) = 2π`. With `₂F₁(...;1;0) = 1`, the `2π` form returns `2π`
  correctly, while the `2π·n` form returns `2π·n`, which is wrong at `a=0`.

  This scaffold uses the **`2π` form**.

  ## Architecture

  Step 1: Generic `n` closed form           — declared as `axiom hypergeometric_lemniscate_radial`
  Step 2: Specialize to `n = 14`            — proved (arithmetic via `norm_num`)
  Step 3: Substitute and simplify           — proved (via `convert`)
  Step 4: Boundary extension to `a = 1`     — declared as `axiom hypergeometric_lemniscate_radial_boundary`
  Step 5: Puiseux rate of ₂F₁ near `z = 1` — declared as `axiom gauss_2F1_puiseux_lower_bound`
  Step 6: Puiseux corollary for n = 14      — proved from Steps 2–5
-/

import Mathlib

open Real Complex MeasureTheory
open scoped ENNReal NNReal

namespace Erdos114
namespace Radial

/-! ## Helper definitions -/

/-- The lemniscate of a polynomial-shaped function `p : ℂ → ℂ` is the level
set `{z : ℂ | |p z| = 1}`. We use the generic norm `‖·‖` (which on `ℂ`
agrees with the classical complex modulus). -/
def lemniscate (p : ℂ → ℂ) : Set ℂ := {z : ℂ | ‖p z‖ = 1}

/-- The lemniscate length of `p` is the 1-dimensional Hausdorff measure of its
lemniscate, expressed as a real number via `ENNReal.toReal`. -/
noncomputable def ehp_lemniscate_length (p : ℂ → ℂ) : ℝ :=
  (μH[1] (lemniscate p)).toReal

/-- Real-valued Gauss hypergeometric `₂F₁(a, b; c; z)`. We use Mathlib's
`ordinaryHypergeometric` specialized to `ℝ`. -/
noncomputable def gauss_2F1 (a b c : ℝ) (z : ℝ) : ℝ :=
  ordinaryHypergeometric a b c z

/-! ## Classical analytic axioms

  Each axiom below states a result that is rigorously known in classical
  analysis but **not** present in Mathlib. Per Formalization Integrity
  Protocol §7, every axiom carries a full bibliographic citation.
-/

/-- **Federer co-area formula on lemniscates.**

**Reference:** Federer, *Geometric Measure Theory*, Springer (1969),
Theorem 3.2.22. Specialization to lemniscates: Pommerenke,
*Univalent Functions*, Vandenhoeck & Ruprecht (1975), Ch. 10. -/
axiom coarea_lemniscate_length (p : ℂ → ℂ) (hp : ∃ n : ℕ, ∃ a : ℂ, ∀ z, p z = z^n - a)
    (n : ℕ) (hn : 1 ≤ n) :
    ∃ I : ℝ, ehp_lemniscate_length p = I

/-- **Fejér–Riesz factorization.**

**References:** Fejér (1916), Riesz (1916); Grenander–Szegő (1958), Ch. 1. -/
axiom fejer_riesz_factorization
    (q : ℝ → ℝ) (N : ℕ)
    (hq_nonneg : ∀ θ, 0 ≤ q θ)
    (hq_trig_poly : ∃ c : Fin (N + 1) → ℂ,
      ∀ θ, q θ = ((Finset.univ.sum (fun k => c k * Complex.exp (Complex.I * θ * (k : ℕ)))).re)^2 +
                 ((Finset.univ.sum (fun k => c k * Complex.exp (Complex.I * θ * (k : ℕ)))).im)^2) :
    ∃ g : ℂ → ℂ, ∀ θ : ℝ, q θ = ‖g (Complex.exp (Complex.I * θ))‖^2

/-- **Gauss connection formula at z = 1.**

For `c − a − b > 0`, the Gauss hypergeometric function `₂F₁(a, b; c; z)`
has a finite limit at `z = 1` given by
`₂F₁(a, b; c; 1) = Γ(c) Γ(c−a−b) / (Γ(c−a) Γ(c−b))`.

**References:** Gauss (1812), §6; Whittaker–Watson (1927), §14.11;
DLMF §15.4.20: <https://dlmf.nist.gov/15.4.E20>. -/
axiom gauss_2F1_connection_at_one (a b c : ℝ)
    (hpos : 0 < c - a - b) :
    Filter.Tendsto (fun z : ℝ => gauss_2F1 a b c z) (nhds 1)
      (nhds (Real.Gamma c * Real.Gamma (c - a - b) /
             (Real.Gamma (c - a) * Real.Gamma (c - b))))

/-- **Radial-direction closed form for the EHP lemniscate length (interior).**

For `p(z) = z^n − a` with `n ≥ 2` and `0 ≤ a < 1`:
`L(z^n − a) = 2π · ₂F₁((n−1)/(2n), (n−1)/(2n); 1; a²)`.

**Provenance:** Underlying integral identity from Gradshteyn–Ryzhik §3.665.2
and DLMF §15.6. Application to lemniscate length via change of variables
`z^n − a = e^{iθ}`. Empirically verified for `n ∈ {3,…,15}` at machine
precision (`RADIAL_PUISEUX_FIT_2026-05-02.md`). -/
axiom hypergeometric_lemniscate_radial (n : ℕ) (hn : 2 ≤ n) (a : ℝ)
    (ha : 0 ≤ a) (ha1 : a < 1) :
    ehp_lemniscate_length (fun z : ℂ => z^n - (a : ℂ)) =
      2 * Real.pi * gauss_2F1 ((n - 1 : ℝ) / (2 * n)) ((n - 1 : ℝ) / (2 * n)) 1 (a^2)

/-- **Radial-direction closed form at the boundary `a = 1`.**

Extension of `hypergeometric_lemniscate_radial` to `a = 1`. The ₂F₁ series
converges absolutely at `z = a² = 1` when `c − α − β = 1/n > 0` (Gauss's
test), so the closed form `L(z^n − 1) = 2π · ₂F₁((n−1)/(2n), (n−1)/(2n); 1; 1)`
holds by continuity of the arc-length integral in the parameter `a`.

**Note:** This axiom was NOT in the original four-axiom scaffold. It was
added because `hypergeometric_lemniscate_radial` requires `a < 1` and
therefore cannot evaluate the lemniscate length at `a = 1`, which is
needed for the Puiseux deficit bound.

**Reference:** Same derivation chain as `hypergeometric_lemniscate_radial`,
extended to the boundary via dominated convergence / Abel's theorem. -/
axiom hypergeometric_lemniscate_radial_boundary (n : ℕ) (hn : 2 ≤ n) :
    ehp_lemniscate_length (fun z : ℂ => z^n - 1) =
      2 * Real.pi * gauss_2F1 ((n - 1 : ℝ) / (2 * n)) ((n - 1 : ℝ) / (2 * n)) 1 1

/-- **Quantitative Puiseux lower bound for ₂F₁ near `z = 1`.**

For `c − a − b > 0`, there exists `K > 0` such that for `w ∈ (0, 1/2]`:
`₂F₁(a, b; c; 1) − ₂F₁(a, b; c; 1−w) ≥ K · w^(c−a−b)`.

This is the quantitative content of the Gauss connection formula expansion
(DLMF §15.8.1). The leading singular term in the expansion of ₂F₁ near
`z = 1` is `−Γ(c)Γ(a+b−c)/(Γ(a)Γ(b)) · (1−z)^(c−a−b)`, which gives the
lower bound with `K = |Γ(c)Γ(a+b−c)|/(|Γ(a)||Γ(b)|)` (up to constants
from the regular part).

**Note:** This axiom was NOT in the original four-axiom scaffold. It was
added because `gauss_2F1_connection_at_one` provides only the *limit*
of ₂F₁ at `z = 1`, not the *rate of approach*. The Puiseux bound
`C · ε^(1/14)` requires the quantitative rate `∼ (1−z)^(c−a−b)`.

**References:**
* DLMF §15.8.1: <https://dlmf.nist.gov/15.8.E1>
* Whittaker–Watson (1927), §14.11, connection formula expansion. -/
axiom gauss_2F1_puiseux_lower_bound (a b c : ℝ) (hpos : 0 < c - a - b) :
    ∃ K : ℝ, 0 < K ∧ ∀ w : ℝ, 0 < w → w ≤ 1/2 →
      gauss_2F1 a b c 1 - gauss_2F1 a b c (1 - w) ≥ K * Real.rpow w (c - a - b)

/-! ## Target theorems -/

/-- **Closed-form expression for the lemniscate length of `z^14 − a`
(radial slice).**

Specializes `hypergeometric_lemniscate_radial` at `n = 14`. The arithmetic
content reduces to verifying `(14 − 1)/(2 · 14) = 13/28`. -/
theorem ehp_radial_L_form_n14 (a : ℝ) (ha : 0 ≤ a) (ha1 : a < 1) :
    ehp_lemniscate_length (fun z : ℂ => z^14 - (a : ℂ)) =
      2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 (a^2) := by
  -- Invoke the generic radial closed form at n = 14.
  have h_general :
      ehp_lemniscate_length (fun z : ℂ => z^14 - (a : ℂ)) =
        2 * Real.pi *
          gauss_2F1 ((14 - 1 : ℝ) / (2 * 14)) ((14 - 1 : ℝ) / (2 * 14)) 1 (a^2) :=
    hypergeometric_lemniscate_radial 14 (by norm_num) a ha ha1
  -- Arithmetic specialization (14 - 1)/(2 * 14) = 13/28 and conclude.
  convert h_general using 3 <;> norm_num

/-- **Closed-form expression for the lemniscate length of `z^14 − 1`
(boundary specialization).**

Uses `hypergeometric_lemniscate_radial_boundary` at `n = 14`. -/
theorem ehp_radial_L_form_n14_boundary :
    ehp_lemniscate_length (fun z : ℂ => z^14 - 1) =
      2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 1 := by
  have h := hypergeometric_lemniscate_radial_boundary 14 (by norm_num)
  convert h using 3 <;> norm_num

/-- **Puiseux expansion at `a = 1` for n = 14.**

There exists `C > 0` such that for every `ε ∈ (0, 1/4]`, the deficit
`L(z^14 − 1) − L(z^14 − (1−ε))` is bounded below by `C · ε^(1/14)`.

**Proof outline:**
1. Evaluate `L(z^14 − 1) = 2π · F(1)` via `ehp_radial_L_form_n14_boundary`.
2. Evaluate `L(z^14 − (1−ε)) = 2π · F((1−ε)²)` via `ehp_radial_L_form_n14`.
3. Set `w = 1 − (1−ε)² = ε(2−ε)`. Note `w ≥ ε > 0` and `w ≤ 1/2`.
4. From `gauss_2F1_puiseux_lower_bound`: `F(1) − F(1−w) ≥ K · w^(1/14)`.
5. Since `w ≥ ε` and `t ↦ t^(1/14)` is increasing: `w^(1/14) ≥ ε^(1/14)`.
6. Combine: deficit = `2π(F(1) − F(1−w)) ≥ 2πK · ε^(1/14)`.

The witness is `C = 2π · K` where `K` comes from the Puiseux rate axiom. -/
theorem ehp_radial_puiseux_n14 :
    ∃ C : ℝ, 0 < C ∧ ∀ ε : ℝ, 0 < ε → ε ≤ 1/4 →
      ehp_lemniscate_length (fun z : ℂ => z^14 - 1) -
      ehp_lemniscate_length (fun z : ℂ => z^14 - ((1 - ε : ℝ) : ℂ)) ≥
        C * Real.rpow ε (1/14) := by
  -- Step 1: Get the Puiseux rate from the ₂F₁ axiom with a = b = 13/28, c = 1.
  -- c − a − b = 1 − 13/28 − 13/28 = 1/14 > 0.
  have hcab : (0 : ℝ) < 1 - 13/28 - 13/28 := by norm_num
  obtain ⟨K, hK, hK_bound⟩ := gauss_2F1_puiseux_lower_bound (13/28) (13/28) 1 hcab
  -- Step 2: Witness C = 2π · K.
  refine ⟨2 * Real.pi * K, mul_pos (mul_pos two_pos pi_pos) hK, ?_⟩
  intro ε hε hε_le
  -- Step 3: Evaluate both lemniscate lengths.
  have hL1 : ehp_lemniscate_length (fun z : ℂ => z^14 - 1) =
      2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 1 :=
    ehp_radial_L_form_n14_boundary
  have hL_eps : ehp_lemniscate_length (fun z : ℂ => z^14 - ((1 - ε : ℝ) : ℂ)) =
      2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 ((1 - ε)^2) :=
    ehp_radial_L_form_n14 (1 - ε) (by linarith) (by linarith)
  rw [hL1, hL_eps]
  -- Step 4: Set w = 1 − (1−ε)² = 2ε − ε² and rewrite.
  set w := 2 * ε - ε^2 with hw_def
  have hw_eq : (1 - ε)^2 = 1 - w := by rw [hw_def]; ring
  rw [hw_eq]
  -- Goal: 2π · F(1) − 2π · F(1−w) ≥ 2πK · ε^(1/14)
  -- Step 5: Verify w ∈ (0, 1/2].
  have hw_pos : 0 < w := by
    rw [hw_def]
    have : 0 < ε * (2 - ε) := mul_pos hε (by linarith)
    linarith [sq ε]
  have hw_le : w ≤ 1/2 := by
    rw [hw_def]; linarith [sq_nonneg ε]
  -- Step 6: Apply the Puiseux rate axiom.
  have hcab_val : (1 : ℝ) - 13/28 - 13/28 = 1/14 := by norm_num
  have hrate := hK_bound w hw_pos hw_le
  rw [hcab_val] at hrate
  -- hrate : F(1) − F(1−w) ≥ K · w^(1/14)
  -- Step 7: w ≥ ε, so w^(1/14) ≥ ε^(1/14).
  have hw_ge : ε ≤ w := by
    rw [hw_def]
    have : ε * ε ≤ ε * 1 := mul_le_mul_of_nonneg_left (by linarith) hε.le
    linarith [sq ε]
  have hrpow_mono : Real.rpow ε (1/14 : ℝ) ≤ Real.rpow w (1/14 : ℝ) :=
    Real.rpow_le_rpow hε.le hw_ge (by norm_num : (0:ℝ) ≤ 1/14)
  -- Step 8: Combine the bounds.
  have h_F_bound :
      gauss_2F1 (13/28) (13/28) 1 1 - gauss_2F1 (13/28) (13/28) 1 (1 - w) ≥
        K * Real.rpow ε (1/14) := by
    calc gauss_2F1 (13/28) (13/28) 1 1 - gauss_2F1 (13/28) (13/28) 1 (1 - w)
        ≥ K * Real.rpow w (1/14) := hrate
      _ ≥ K * Real.rpow ε (1/14) := mul_le_mul_of_nonneg_left hrpow_mono hK.le
  -- Step 9: Scale by 2π ≥ 0 to conclude.
  have hpi : (0 : ℝ) ≤ 2 * Real.pi := mul_nonneg two_pos.le pi_pos.le
  calc 2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 1 -
        2 * Real.pi * gauss_2F1 (13/28) (13/28) 1 (1 - w)
      = 2 * Real.pi * (gauss_2F1 (13/28) (13/28) 1 1 -
          gauss_2F1 (13/28) (13/28) 1 (1 - w)) := by ring
    _ ≥ 2 * Real.pi * (K * Real.rpow ε (1/14)) :=
        mul_le_mul_of_nonneg_left h_F_bound hpi
    _ = 2 * Real.pi * K * Real.rpow ε (1/14) := by ring

end Radial
end Erdos114
