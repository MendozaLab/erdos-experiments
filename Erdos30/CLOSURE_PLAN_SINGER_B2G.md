# Closure Plan — `singer_b2g_exists` (Bose–Chowla 1962)

**Scope:** Research-tier (T3) scoping document. **Not** a closure attempt.
**Date:** 2026-05-02
**File:** `lean/Erdos755_Singer_BhG.lean:69`
**Sidecar status:** Out of #30 publication scope. The sidecar package (Erdős #755 / B_h[g] generalization) is companion infrastructure to the main Sidon (#30) bundle, not part of any active filing or preprint.

---

## 1. Axiom statement (verbatim)

```lean
/-- **AXIOM — Singer/Bose–Chowla construction for B_2[g].**

For every g ≥ 1 and every sufficiently large N, there exists a B_2[g]
set A ⊆ [0, N] with |A|² ≥ g · (2N + 1) / C for some absolute constant C.

This matches the Lindström upper bound |A|² ≤ g(2N+1) + O((gN)^{3/4})
up to a constant factor.

Proof (in the literature, NOT formalized here):
  Take q = prime power ≈ √(N/g). In GF(q²), pick primitive element α.
  Set A = {i : α^i ∈ GF(q) ⊆ GF(q²)} under a suitable embedding.
  Verify |A| = q and the B_2[g] property via multiplicative structure.

Status: Classical result. Formalization in Lean requires:
  - Prime-power existence (`Nat.exists_prime_pow_near`)
  - GF(q^2) primitive element (Mathlib has `IsPrimitiveRoot`)
  - Embedding ℤ/q²ℤ ↪ ℕ preserving Sidon property
  - Quadratic reciprocity for the g > 1 extension -/
axiom singer_b2g_exists (N g : ℕ) (hg : g ≥ 1) (hN : N ≥ 16) :
    ∃ (A : Finset ℕ), IsB2GSet A g ∧ A ⊆ Finset.range (N + 1) ∧
      g * (2 * N + 1) ≤ 4 * A.card * A.card
```

(`IsB2GSet A g` is the ordered-form definition: every sum has at most `2g`
ordered representations as `(a, b) ∈ A × A`. See `Erdos755_BhG.lean:36`.)

---

## 2. Where it is consumed

Two downstream sites, both inside the #755 sidecar package — never imported
by the #30 main bundle:

| Theorem | File:Line | Role |
|---|---|---|
| `b2g_tight_bound` | `Erdos755_Singer_BhG.lean:83` | Pairs the axiom (lower) with `b2g_card_sq_bound` (upper) to expose Θ(gN) tightness |
| `sidon_specialization` | `Erdos755_Singer_BhG.lean:97` | g=1 instantiation: 4·\|A\|² ≥ 2N+1 |
| `b2g_optimal_density` | `Erdos755_Complete.lean:59` | Headline #755 result combining T1 (sum count, proven) + T4 (this axiom) |
| `sidon_optimal_density` | `Erdos755_Complete.lean:73` | g=1 of the headline |

If the axiom is removed, all four theorems collapse — but nothing in the #30
publication artifact depends on any of them.

---

## 3. Bose–Chowla proof outline (in the literature)

**Setup.** Pick a prime power `q`. Let `θ` be a primitive element of the
multiplicative group `GF(q^h)*`, which is cyclic of order `q^h − 1`. The
Bose–Chowla 1962 construction sets

```
A = { i ∈ {0, 1, …, q^h − 2} : θ^i − 1 ∈ GF(q) }   (or an isomorphic image)
```

For the `h = 2` case relevant here, one takes `GF(q^2)` and considers the
image under the canonical map `ψ : GF(q^2)* → ℤ/(q^2 − 1)ℤ` of the line
`{1, θ, θ², …, θ^k}` for `k = q + 1`, where `θ` is a primitive root.

**Sidon property (g = 1).** The proof rests on the fact that distinct sums
`a + b` correspond to distinct multiplicative pairs because `GF(q^h)` is an
integral domain — any collision in the additive structure would force a
nontrivial polynomial relation among `q + 1` distinct field elements, which
the primitive root rules out by minimality of its order.

**B_2[g] extension.** Bose–Chowla observed that allowing `g` field-level
collisions yields exactly `2g` ordered pair preimages, recovering the
ordered-form `IsB2GSet A g` predicate. Quadratic reciprocity enters only
when one wants explicit constants `C` better than `4`; the existence claim
only needs the underlying primitive-root cardinality count.

**Density.** The construction yields `|A| = q + 1`, and a counting argument
gives `q² + q ≤ N`, so `|A|² ≥ q² ≥ N − q ≈ N`. With the `4 ·` slack in the
axiom statement, the constant works out cleanly for any `g ≥ 1` and
`N ≥ 16`.

**Reference.** R. C. Bose and S. Chowla, *Theorems in the additive theory
of numbers*, Comment. Math. Helv. 37 (1962), 141–147.

---

## 4. Dependency on `singer_sidon_exists`

`singer_b2g_exists` is the natural `B_2[g]` generalization of
`singer_sidon_exists` (`Erdos30_Singer.lean:177`):

| Axiom | Range | Construction |
|---|---|---|
| `singer_sidon_exists` | g = 1, every prime q | Singer 1938 cycle in PG(2, q) → perfect difference set in ℤ/(q²+q+1) → Sidon set of size q+1 in `[0, q²+q]` |
| `singer_b2g_exists` | g ≥ 1, sufficiently large N | Bose–Chowla 1962 line-image in GF(q²) → B_2[g] set of size ≈ q in `[0, N]` |

The `g = 1` Bose–Chowla construction is essentially the **same object** as
the Singer perfect difference set, viewed through the multiplicative group
of `GF(q²)` rather than the projective plane PG(2, q). Concretely,
`singer_sidon_exists` already gives a Sidon set of size `q + 1` in
`[0, q² + q]`; `singer_b2g_exists` at `g = 1` only needs the
density-rescaling step `q² + q ≤ N → q ≥ √N − O(1)`.

**Implication:** Once `singer_sidon_exists` is closed in Lean, the `g = 1`
specialization `sidon_specialization` becomes a derived theorem rather than
an axiom corollary. The `g > 1` case still requires the Bose–Chowla
extension argument. So closing Singer first **partially** unblocks this
axiom — it would let `Erdos755_Singer_BhG.sidon_specialization` discharge
without referencing `singer_b2g_exists` at all, but `b2g_optimal_density`
for `g ≥ 2` would still need the full construction.

---

## 5. Mathlib gaps

The full closure needs four pieces, in dependency order:

1. **Prime-power density:** a constructive form of `Nat.exists_prime_pow_near`
   yielding `q` with `q² ≤ N/g < (q+1)²` (or similar). Mathlib has
   `Nat.exists_prime_lt_and_le_two_mul` and Bertrand's postulate but not
   the prime-power refinement; manageable to derive.

2. **`GaloisField q 2` primitive element:** Mathlib has `GaloisField` and
   `IsPrimitiveRoot` but lacks the explicit cyclic-generator extraction
   for the multiplicative group of a finite field. (Same gap blocks
   `singer_sidon_exists`.)

3. **Trace/embedding `GF(q²) → ℤ/(q²−1)ℤ`:** Mathlib `Algebra.trace` is
   partial; the Bose–Chowla construction needs a discrete-log-style index
   map preserving multiplicative structure, which is **not formalized**
   anywhere in Lean / Coq / Isabelle / Mizar (Perplexity gate 2026-03-29
   logged in `AXIOM_INVENTORY.md` for `singer_sidon_exists`).

4. **Sum-collision counting in GF(q²):** the `B_2[g]` upgrade over Sidon
   needs a polynomial-degree argument bounding additive collisions by
   multiplicative degree — straightforward once (1)–(3) are in place, but
   requires a clean `Polynomial.card_roots`-style bridge.

The dominant blocker is (3): the Singer/Bose–Chowla cycle construction is
**not formalized in any prover**. This is the same chokepoint flagged for
`singer_sidon_exists` in the canonical axiom inventory.

---

## 6. Effort estimate

Assuming `singer_sidon_exists` ships first (so (1)–(3) are infrastructure
that already exists):

| Stage | Agent-cycles |
|---|---|
| Prime-power density refinement (gap 1) | shared with Singer closure |
| Bose–Chowla line-image construction (gap 4 + index map specialization) | 8–12 |
| `IsB2GSet` ordered-form proof from multiplicative collision count | 4–6 |
| `b2g_tight_bound` and `sidon_specialization` cleanup once axiom drops | 1–2 |
| **Total (post-Singer)** | **~13–20 agent-cycles** |

If `singer_sidon_exists` is **not** closed first, multiply by ~3× because
the GaloisField + Singer-cycle infrastructure has to be built from scratch
— and that is the multi-week piece flagged as research-tier in the
`AXIOM_INVENTORY` (entry #4, T3).

Standalone closure (no Singer prereq): **~40–60 agent-cycles** spread over
multiple sessions, with the bulk going to Mathlib infrastructure rather
than B_h[g]-specific reasoning.

---

## 7. Recommendation

**Defer indefinitely.** Three reinforcing reasons:

1. **Out of publication scope.** The #755 sidecar is not part of the #30
   publication artifact; no Zenodo / arXiv / Cooley filing depends on
   closing this axiom. The `AXIOM_INVENTORY.md` already records it as
   stable T3 axiomatization of a 1962 published result with explicit
   citation — meeting the H² Formalization Integrity Protocol §7 honest
   scope bar.

2. **Strict dependency on Singer.** The `g = 1` case is genuinely the
   same object as `singer_sidon_exists`; the `g ≥ 2` extension only adds
   value once Singer is in. There is no path where closing
   `singer_b2g_exists` first makes sense — it would force a redundant
   re-derivation of the Singer/Bose–Chowla cycle.

3. **Shared chokepoint upstream.** The Singer cycle construction is
   nowhere formalized in any prover; the gating effort is a Mathlib-PR
   contribution to `GaloisField` infrastructure, not a B_h[g]-specific
   project. That work belongs in the research-grant cycle, not in
   sidecar maintenance.

**Conditional priority bump:** If and only if Ken closes
`singer_sidon_exists` (or upstream Mathlib lands the GaloisField primitive
element + Singer cycle infrastructure), revisit this plan and execute
Section 6 stages 2–4. Until then, leave the axiom in place with its
current citation block. Expected priority: **lowest in the T3 queue**,
behind `bfr_core_bound`, `singer_sidon_exists`, and `lindstrom_sieve`.

---

## 8. Provenance

- Axiom file: `lean/Erdos755_Singer_BhG.lean` (axiom at line 69, body lines 51–71)
- Consumers: `lean/Erdos755_Singer_BhG.lean` (lines 83, 97), `lean/Erdos755_Complete.lean` (lines 59, 73)
- Companion axiom: `singer_sidon_exists`, `lean/Erdos30_Singer.lean:177`
- Inventory: `Erdos30/AXIOM_INVENTORY.md` (entry #5)
- This plan: scoping only. No `.lean` files modified. No `sorry` / `admit` introduced.
