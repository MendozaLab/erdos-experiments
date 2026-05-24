# Closure Attempt: `lindstrom_bound` (T2 axiom)

**Date:** 2026-05-02
**File:** `lean/Erdos30_Lindstrom.lean:521`
**Status:** **DEFERRED — statement-level obstruction (not a tactic-level gap)**
**Outcome:** Axiom left in place; theoretical analysis below shows the stated form does not follow from `lindstrom_quadratic` alone.

## Target

```lean
axiom lindstrom_bound (A : Finset ℕ) (N : ℕ) (hS : IsSidonSet A)
    (hA : A ⊆ Finset.range (N + 1)) (hN : 0 < N) :
    A.card ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 1
```

## Strategy attempted

Set `s = Nat.sqrt N` (so `s² ≤ N ≤ s² + 2s`) and `ℓ = Nat.sqrt s` (so `ℓ² ≤ s ≤ ℓ² + 2ℓ`). Apply `lindstrom_quadratic A N ℓ ...` with `ℓ ≥ 1` to get
```
ℓ · (2k − ℓ − 1)² ≤ 4 · (ℓ + 1) · N
```
By contradiction, assume `k ≥ s + ℓ + 2`. Then `2k − ℓ − 1 ≥ 2s + ℓ + 3`, so the LHS is at least
```
ℓ · (2s + ℓ + 3)²  =  4ℓs² + 4ℓs(ℓ + 3) + ℓ(ℓ + 3)².
```
The RHS, using `N ≤ s² + 2s`, is at most
```
4(ℓ + 1)(s² + 2s)  =  4ℓs² + 4s² + 8ℓs + 8s.
```
Subtracting and simplifying:
```
LHS − RHS  ≥  4s(ℓ − 1)(ℓ + 2) + ℓ(ℓ + 3)² − 4s²        (★)
```
For the contradiction we need (★) > 0, i.e.,
```
ℓ(ℓ + 3)²  >  4s² − 4s(ℓ − 1)(ℓ + 2)  =  4s · (s − ℓ² − ℓ + 2).
```
In the worst case `s = ℓ² + 2ℓ` (largest allowed `s` for a given `ℓ`), the RHS becomes
```
4(ℓ² + 2ℓ)(ℓ + 2)  =  4ℓ(ℓ + 2)².
```
The required inequality is then `ℓ(ℓ + 3)² > 4ℓ(ℓ + 2)²`, i.e., `(ℓ + 3)² > 4(ℓ + 2)²`, equivalently `ℓ + 3 > 2ℓ + 4`, i.e., `ℓ < −1` — **false for every nonnegative ℓ**.

In fact `4ℓ(ℓ + 2)² − ℓ(ℓ + 3)² = ℓ(ℓ + 1)(3ℓ + 7) > 0` for `ℓ ≥ 1`.

## Concrete worked counterexample to the proof strategy

Take `N = 15`. Then `s = Nat.sqrt 15 = 3`, `ℓ = Nat.sqrt 3 = 1`. The Lean statement claims `k ≤ 5`. Plugging into `lindstrom_quadratic` with `ℓ = 1`:
```
1 · (2k − 2)²  ≤  4 · 2 · 15  =  120
⟹ 2k − 2 ≤ 10
⟹ k ≤ 6.
```
With `ℓ = 2`: `2 · (2k − 3)² ≤ 4 · 3 · 15 = 180 ⟹ (2k − 3)² ≤ 90 ⟹ 2k − 3 ≤ 9 ⟹ k ≤ 6`.
With `ℓ = 3`: `3 · (2k − 4)² ≤ 4 · 4 · 15 = 240 ⟹ (2k − 4)² ≤ 80 ⟹ 2k − 4 ≤ 8 ⟹ k ≤ 6`.

So `lindstrom_quadratic` (the only ingredient available below the axiom line) forces only `k ≤ 6`, while the Lean statement claims `k ≤ 5`. The empirical fact that no Sidon set of size 6 fits in `[0,15]` (max is 5, e.g. `{0,1,3,7,12}`) does not follow from `lindstrom_quadratic` — it is a **strictly stronger** combinatorial fact.

## Root cause

The classical real-arithmetic Lindström bound is `k ≤ √N + N^{1/4} + 1` (with the floor giving `⌊√N + N^{1/4}⌋ + 1` when integerized). The Lean statement's
```
Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 1
```
is **strictly tighter** than `⌊√N + N^{1/4}⌋ + 1` whenever `√N` and `⁴√N` both have substantial fractional parts. For `N = 15`, the real bound gives `k ≤ ⌊6.84⌋ = 6`, but the ℕ statement claims `k ≤ 5`. This means closing the axiom in the form written requires combinatorial input **beyond** `lindstrom_quadratic`.

## What WOULD close from `lindstrom_quadratic` (provable in pure ℕ)

A weaker variant is provable by the real-arithmetic chain plus a Real-to-Nat round-trip:
```
k ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 2     (one-extra slack)
```
or, with a small constant absorbed into the rounding,
```
k ≤ Nat.sqrt (4 * N) + Nat.sqrt (Nat.sqrt (4 * N)) + 1
```
These would close from `lindstrom_quadratic` plus `Real.sqrt` arithmetic, but the file already uses `Real.sqrt` in `Erdos30_FaceField_60_FullFace_Certificate.lean`, so the dependency is already present. A `+1` slack version is the honest provable form.

## Mathlib lemmas examined

- `Nat.sqrt_le n : sqrt n * sqrt n ≤ n`
- `Nat.lt_succ_sqrt n : n < (sqrt n + 1) * (sqrt n + 1)` ⟹ `n ≤ s² + 2s`
- `Nat.le_sqrt : m ≤ sqrt n ↔ m * m ≤ n`
- `Nat.sqrt_lt : sqrt m < n ↔ m < n * n`
- `Nat.eq_sqrt`, `Nat.sqrt_le_sqrt`, `Nat.sqrt_le_add`

None of these change the underlying issue: the ℕ statement asserts a tighter inequality than is derivable from `lindstrom_quadratic`'s `ℓq² ≤ 4(ℓ+1)N`.

## Sharper combinatorial form available, but still insufficient

Going one level deeper, `order_diff_counting + sum_distinct_pos_ge` actually deliver the tighter
```
m(m + 1) ≤ ℓ(ℓ + 1)N           where 2m = ℓ(2k − ℓ − 1).
```
Multiplying by 4:
```
ℓ(2k − ℓ − 1)·(ℓ(2k − ℓ − 1) + 2)  ≤  4ℓ(ℓ + 1)N
⟹ (2k − ℓ − 1)·(ℓ(2k − ℓ − 1) + 2)  ≤  4(ℓ + 1)N.
```
For `N = 15`, `ℓ = 1`: `(2k − 2)·((2k − 2) + 2) = (2k − 2) · 2k ≤ 120`, giving `2k(2k − 2) ≤ 120 ⟹ 4k² − 4k ≤ 120 ⟹ k² − k ≤ 30 ⟹ k ≤ 6`. Still 6, not 5.

So even the **tightest** consequence of `order_diff_counting` does not deliver the stated bound for boundary `N`.

## Recommendation

**Three closure paths, in order of cost:**

1. **(Cheapest, recommended.)** Restate the theorem with one extra unit of slack:
   ```lean
   theorem lindstrom_bound_relaxed ... :
       A.card ≤ Nat.sqrt N + Nat.sqrt (Nat.sqrt N) + 2
   ```
   This follows from `lindstrom_quadratic` plus the `Real.sqrt` arithmetic chain (already imported in the package). Perhaps 60–80 lines.

2. **(Medium.)** Keep the current statement but close it via a **sharper combinatorial argument** that goes beyond `order_diff_counting` — specifically, an extension that exploits the position of the smallest and largest elements of `A` to refine the telescoping sum. This is what Lindström's 1969 paper does in its full form (as opposed to the simplified BFR §2 version axiomatized here). Estimated ~150 lines plus a strengthened axiom statement.

3. **(Expensive.)** Brute-force enumeration for small `N` (up to whatever boundary makes the real bound and the ℕ bound coincide), combined with the real-arithmetic chain for large `N`. Probably 200+ lines and depends on `decide`-style automation working at that scale.

## Verdict

**Axiom NOT closed in this session.** The strategy outlined in the task description — apply `lindstrom_quadratic` with `ℓ = Nat.sqrt(Nat.sqrt N)` and round via `Nat.sqrt` — provably **cannot** close the axiom as stated, because in worst-case `N` the quadratic alone leaves a unit gap that the stated bound does not afford.

**Build state:** Unchanged. `lake build Erdos30_Lindstrom` PASS.
**Axiom count:** Still 4 (bfr_core_bound, order_diff_counting, lindstrom_bound, singer_sidon_exists).
**Recommended next action for Ken:** Decide between (1) accepting the +1 slack and closing the relaxed form from `lindstrom_quadratic` + Real.sqrt, or (2) committing to the sharper combinatorial argument that closes the original statement.
