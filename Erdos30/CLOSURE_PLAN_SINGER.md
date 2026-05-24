# CLOSURE_PLAN_SINGER — Erdős #30 Singer-axiom retirement plan

**Date:** 2026-05-02
**Owner:** Math/erdos-experiments/Erdos30 working group
**Tier:** T3 (research-tier formalization)
**Status of axiom:** OPEN (still in `lean/Erdos30_Singer.lean:177`)
**Status of file:** `lake build Erdos30_Singer` PASS @ Mathlib v4.27.0, lean v4.27.0
**Novelty:** No formal proof of Singer's theorem in Lean / Coq / Isabelle / Mizar (Perplexity 2026-03-29).

---

## 1. Verbatim axiom statement

```lean
axiom singer_sidon_exists (q : ℕ) (hq : Nat.Prime q) :
    ∃ A : Finset ℕ, IsSidonSet A ∧ A.card = q + 1 ∧ ∀ a ∈ A, a ≤ q * q + q
```

`IsSidonSet` (from `Erdos30_Sidon_Defs.lean`):

```lean
abbrev IsSidonSet (A : Finset ℕ) : Prop :=
  ∀ a₁ ∈ A, ∀ b₁ ∈ A, ∀ a₂ ∈ A, ∀ b₂ ∈ A,
    a₁ ≤ b₁ → a₂ ≤ b₂ → a₁ + b₁ = a₂ + b₂ → (a₁ = a₂ ∧ b₁ = b₂)
```

The axiom is consumed by `singer_lower_bound` and indirectly by the h(N) ≥ √N
direction of Erdős #30. It is currently the only non-explicit axiom in the
Erdős #30 formalization; closing it is the last step toward a fully axiom-free
Singer lower-bound chain (modulo Lindström, which is a separate workstream).

A note on scope: Singer's theorem is stated for prime *powers* q. Our axiom
restricts to primes for the specific use in Erdős #30; the closure plan below
naturally targets the prime-power version, then specializes.

---

## 2. Mathematical proof outline (Singer 1938)

For prime q (or prime power), the projective plane PG(2, q) is the set of
1-dimensional subspaces of (GF(q))³. It has |PG(2, q)| = q² + q + 1 points.
A line in PG(2, q) is the projectivization of a 2-dim subspace; each line has
q + 1 points. Any two distinct lines meet in exactly one point.

The multiplicative group GF(q³)* has order q³ − 1 = (q − 1)(q² + q + 1) and is
cyclic. Let α generate GF(q³)*; let β = α^{q−1}. Then β has order q² + q + 1.

GF(q³) is a 3-dim vector space over GF(q). Multiplication by α (or β) is a
GF(q)-linear automorphism, so it descends to a permutation σ on PG(2, q). The
order of σ on PG(2, q) is q² + q + 1: σ acts transitively on points and is
cyclic of full order. This σ is the Singer cycle.

Pick any line L₀ ⊂ PG(2, q). Number the q² + q + 1 points by powers of σ
applied to a chosen base point. Let

  D = { i ∈ Z_{q²+q+1} : σ^i(p₀) ∈ L₀ }

be the set of indices of points of L₀. Then |D| = q + 1, and D is a *perfect
difference set* mod q² + q + 1: every nonzero residue d mod q² + q + 1 is
realized as exactly one ordered difference i − j with i, j ∈ D.

Reason: σ acts simply-transitively on point–line incidence pairs after some
identification (via the orbit of incidence), so the multiset of nonzero
differences from D is exactly Z_{q²+q+1} \ {0}, each with multiplicity 1.
Equivalently: any two points of PG(2, q) lie on a unique line, so every nonzero
"shift" σ^d takes L₀ to a line meeting L₀ in exactly one point, which forces
exactly one (i, j) pair with i − j = d.

Perfect difference set ⇒ Sidon set in ℕ: if a₁ + b₁ = a₂ + b₂ with all four in
D ⊂ {0, …, q²+q}, then a₁ − a₂ = b₂ − b₁ as integers; reducing mod q²+q+1, this
is the same nonzero residue (or zero). Zero forces a₁ = a₂, b₁ = b₂. Nonzero
contradicts the perfect-difference-set property (which says each nonzero
residue is hit exactly once as an ordered difference). Done.

---

## 3. Four-step Lean formalization plan

| Step | Statement | Mathlib has | Mathlib missing | Effort |
|---|---|---|---|---|
| (i) | GF(q³) exists, has cardinality q³, and its unit group is cyclic of order q³−1 | `GaloisField p n` (Mathlib/FieldTheory/Finite/GaloisField.lean), `instance [Finite Rˣ] : IsCyclic Rˣ` (Mathlib/RingTheory/IntegralDomain.lean), `GaloisField.finrank` | Smooth specialization to q prime: need `∃ α : (GaloisField q 3)ˣ, orderOf α = q^3 − 1` as a usable form. Mathlib gives `IsCyclic.exists_generator` directly. | LOW (1 agent-cycle: rename + glue) |
| (ii) | The Singer cycle σ: ℙ_{GF(q)} (GF(q)³) → ℙ_{GF(q)} (GF(q)³) given by [v] ↦ [α · v] is a bijection of order q² + q + 1 | `Projectivization` (Mathlib/LinearAlgebra/Projectivization/Basic.lean), `Projectivization.smul_mk` (Action.lean), cardinality `q²+q+1` (Cardinality.lean: `card_of_finrank` gives ∑ q^i over range 3) | The order-of-σ argument: needs (a) GF(q)³ ≃ₐ GF(q³) as GF(q)-modules — provable via `Module.finrank_eq_iff_iso` plus the fact that finite-dim spaces of same dim over a field are iso; (b) the action of α descends modulo Kˣ to ℙ; (c) the order of the descended action is exactly (q³−1)/(q−1) = q²+q+1. The order computation uses orbit-stabilizer + the explicit kernel of the descent map (Kˣ ↪ (GF(q³))ˣ). | MEDIUM-HIGH (3–5 agent-cycles: this is the central step; the descent-and-order argument is roughly 100–200 lines of careful module/group theory; needs a small auxiliary file) |
| (iii) | Pick a line L₀ ⊂ PG(2, q) and number its q+1 points as σ-powers of a base point: D = {i₀, i₁, …, i_q} ⊂ Z_{q²+q+1} | `Projectivization.Subspace` (Mathlib/LinearAlgebra/Projectivization/Subspace.lean) defines projective subspaces; cardinality of a line in ℙ over a finite field follows from `card_of_finrank` applied to a 2-dim subquotient | "Line" as the projectivization of a 2-dim subspace is in Mathlib but the specific formula |line| = q+1 needs to be wrapped. The σ-orbit numbering of points on L₀ is mostly bookkeeping (orbit of σ on PG(2, q), restriction to L₀). | LOW-MEDIUM (1–2 agent-cycles: pure definitional packaging once (ii) is done) |
| (iv) | D is a perfect difference set mod q²+q+1, hence (when realized in ℕ via canonical embedding Z_{q²+q+1} ↪ {0,...,q²+q}) is a Sidon set in ℕ with the bound | None — this is the heart of Singer's argument. Counting: 56 ordered differences (q² + q in general) hit q²+q nonzero residues exactly once (pigeonhole + uniqueness of line through two points). Then the PDS → Sidon-in-ℕ bridge is a 5-line modular-arithmetic lemma. | All missing in formal form: (a) the simple-transitivity argument that the multiset of differences in D equals Z\\{0} as a multiset (this is the "two points lie on a unique line, hence every nonzero shift produces exactly one collision" lemma — fundamentally an incidence-geometry argument that has no analog in Mathlib); (b) the PDS → Sidon-in-ℕ bridge (very short, but needs `ZMod n → ℕ` lift). | HIGH (4–6 agent-cycles: (a) is the deepest piece; (b) is trivial) |
| **Total** | | | | **9–14 agent-cycles** ≈ 4–7 weeks of focused human-equivalent effort |

The dominant subgraph of the work is steps (ii) + (iv). They are mostly
independent: (ii) is finite-field/group-theory; (iv) is finite-incidence-geometry
plus modular arithmetic. Run them as parallel tracks.

---

## 4. Effort summary

- **Agent-cycle estimate:** 9–14 cycles
- **Critical-path length:** ~5 cycles (step (ii) is the longest single chain)
- **Parallelizable:** steps (i)+(ii) and (iv) can run as two independent tracks; (iii) joins them
- **Human-equivalent:** 4–7 weeks of dedicated formalization effort, comparable to a small Mathlib PR (e.g. one new file ~500–800 lines plus enabling lemmas in Projectivization)
- **Confidence:** medium-high. Proof structure is classical and well-understood; risk is mostly in step (iv) part (a), which requires a cleanly stated incidence lemma that Mathlib has not yet abstracted

---

## 5. Direct port vs. Mathlib contribution — recommendation: **Mathlib contribution (with a private staging file)**

| Option | Pros | Cons |
|---|---|---|
| Direct port (private file in `Erdos30/lean/`) | Fastest path to closing the axiom; no PR review cycle; full control of API choices | Reinvents Projectivization-incidence machinery in private; not reusable by Mathlib downstreams; no community visibility for the formalization (which is the *novel* part — Singer is not formalized anywhere) |
| Mathlib PR ("Mathlib/Combinatorics/SingerDifferenceSet.lean" + a small extension to `Projectivization/Subspace.lean`) | First formal Singer's theorem in any prover — publishable as a standalone PR-paper ("Formalizing Singer's 1938 theorem in Lean"); reusable for design-theory work in Mathlib (BIBDs, projective planes, finite geometries); keeps the Erdős-30 module pure | Slower (review + back-and-forth ~2x effort); may surface API design questions in `Projectivization.Subspace` that need broader Mathlib-team buy-in |

**Recommendation:** Mathlib contribution path, with a private staging file
inside `Erdos30/lean/` that keeps the build green during the PR cycle. Once the
upstream PR lands, retire the private file and import from Mathlib. The
community visibility, the standalone publication ("first formalization of
Singer 1938 in any prover"), and the reusability for downstream design-theory
work all favor this path. The Erdős #30 axiom-closure is the *consumer* of the
work, not its main artifact.

A useful intermediate publication is the *partial* result already in
`Erdos30_Singer.lean`: explicit verification for q ∈ {2, 3, 5, 7, 11, 13} via
`native_decide`, together with the four-step plan above. That is publishable
today as "Extended Singer construction verified for q ≤ 13."

---

## 6. Intermediate target — `native_decide` extension to more primes

The pattern in `lean/Erdos30_Singer.lean` for q ∈ {2, 3, 5} works directly for
any q for which a published Singer (or planar) difference set list is
available. The state space for the `IsSidonSet` predicate is q⁴ unordered
quadruples, well within `native_decide` reach for q ≲ 20. Beyond that, kernel
runtime starts to matter and we may need a more efficient predicate (e.g.
`Finset.image (·.1 + ·.2) (A ×ˢ A)` size check) — but for q ≤ 19 the present
form is fine.

**Discharged this session (2026-05-02):**

| q | n = q²+q+1 | Set | Verified |
|---|---|---|---|
| 2 | 7   | {0, 1, 3} | already in file |
| 3 | 13  | {0, 1, 3, 9} | already in file |
| 5 | 31  | {0, 1, 3, 8, 12, 18} | already in file |
| **7**  | **57**  | **{0, 1, 3, 13, 32, 36, 43, 52}** | **NEW — added & built** |
| **11** | **133** | **{0, 1, 3, 12, 20, 34, 38, 81, 88, 94, 104, 109}** | **NEW — added & built** |
| **13** | **183** | **{0, 1, 3, 16, 23, 28, 42, 76, 82, 86, 119, 137, 154, 175}** | **NEW — added & built** |

**Build state:** `lake build Erdos30_Singer` PASS in ~13s on baseline +
13/11/7 layers; all six theorems `singer_sidon_q{2,3,5,7,11,13}` produced
witnesses to the axiom statement at those specific primes.

**Next-q candidates (for future sessions, all from published planar-DS lists):**

- q = 17, n = 307 — published Singer DS exists (length 18)
- q = 19, n = 381 — published Singer DS exists (length 20)
- q = 23, n = 553 — published Singer DS exists (length 24)
- (Note that q = 25 = 5² is not prime; the axiom as stated does *not* cover it.
  Singer's theorem proper does, and would land in the Mathlib-PR path of §5.)

Each new q adds two `native_decide` calls (one for `_sidon`, one is `decide`
for `_card`/`_range`); per-q cost grows roughly as q⁴ but with low constants
in the LEAN kernel; q = 19 should still fit comfortably in a single build.

---

## 7. Alignment with Math/CLAUDE.md priority stack

- **Protect** — closing the Singer axiom does NOT disclose ErdosAtlas-novel
  methods (transfer operators, MDL classification, morphism index). Singer
  1938 is fully classical mathematics. Cooley filter clears.
- **Articulate** — the *partial* result (q ∈ {2, 3, 5, 7, 11, 13} discharged)
  is publishable: "Extended Singer construction verified for q ≤ 13," DOI via
  Zenodo, social via Mathstodon + Lean Zulip. The full Mathlib PR, when
  eventually landed, would be a separate publication.
- **Verify** — `lake build Erdos30_Singer` PASS confirms verification.
- **Connect** — Singer construction in Mathlib would connect to
  `Mathlib.Combinatorics.PerfectDifferenceSet` (NEW), `Mathlib.Combinatorics.Sidon`
  (NEW or extended), `Mathlib.LinearAlgebra.Projectivization.Subspace` (extended
  with line-cardinality), and downstream BIBD / finite-geometry work.
- **Solve** — the axiom stays open *as a Mathlib target*, not as a private
  workaround. We have 6 explicit primes discharged; the general case is now
  scoped at 9–14 agent-cycles with a clear Mathlib path.

---

## 8. References

- Singer, J. (1938). "A theorem in finite projective geometry and some
  applications to number theory." *Trans. AMS* 43(3), 377–385.
- Erdős, P. & Turán, P. (1941). "On a problem of Sidon in additive number theory."
  *J. London Math. Soc.* 16, 212–215.
- O'Bryant, K. (2004). "A complete annotated bibliography of work related to
  Sidon sets." *Electron. J. Combin.* DS11.
- La Jolla Difference Set Repository — published planar difference sets for
  q = 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, …
  (Used to source the q = 7, 11, 13 candidate sets above.)
- Mathlib v4.27.0:
  `Mathlib/FieldTheory/Finite/GaloisField.lean`,
  `Mathlib/RingTheory/IntegralDomain.lean` (`[Finite Rˣ] : IsCyclic Rˣ`),
  `Mathlib/LinearAlgebra/Projectivization/Basic.lean`,
  `Mathlib/LinearAlgebra/Projectivization/Cardinality.lean`,
  `Mathlib/LinearAlgebra/Projectivization/Subspace.lean`,
  `Mathlib/LinearAlgebra/Projectivization/Action.lean`.

---

## 9. Decision queue for Ken

1. Approve Mathlib-PR path (§5) for the general-q closure? Or stay on direct
   port to keep velocity?
2. Authorize Zenodo DOI for the partial result ("Extended Singer construction
   verified for q ≤ 13" — 6 explicit primes, axiom remains for general case)?
   Companion to existing Erdős #30 preprint, not a standalone paper.
3. Schedule the q ∈ {17, 19, 23} extension for a future session (low-hanging,
   ~30 minutes per q once the published DS is sourced)?
