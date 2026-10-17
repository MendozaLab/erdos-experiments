# Erdős 123, snug variant for 2, 3, 5: Lean bridge

This package proves the Formal Conjectures statement
`Erdos123.erdos_123.variants.powers_2_3_5_snug` exactly as written in
[google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/ErdosProblems/123.lean).

**The mathematics is Principia Math's.** Their band theorem says that for pairwise coprime
a, b, c > 1 and every real 1 < ρ < min(a, b, c), every large integer is a sum of distinct
elements of {a^k b^l c^m} ∩ [x, ρx) for some x, none dividing another. Their note
*Erdős 123⁺* records the snug corollary for (2, 3, 5). Both are in
[antoshashakov/Principia-Math-Solutions](https://github.com/antoshashakov/Principia-Math-Solutions/tree/30d7e5327cdf5cda55137d70abec8f4d3b7ff2b3/erdos123)
(Apache-2.0), with a complete Lean 4 formalization.

This package adds only `SnugBridge.lean` (about 70 lines). It specialises their theorem
`Erdos123.Statement.erdos123_dcomplete_real` to (a, b, c) = (2, 3, 5) with
ρ = min(1 + ε, 3/2), and rewrites the conclusion in the Formal Conjectures vocabulary.
`IsSnug` and the theorem statement are verbatim copies of the upstream file; the copies were
checked character for character against upstream `main`.

## Build

```sh
lake exe cache get
lake build
```

Lean `v4.31.0`, Mathlib `db127794c79fdeb86f6b0cf6ff2c804026fbaff1`, Principia pinned at
`30d7e5327cdf5cda55137d70abec8f4d3b7ff2b3`. See `VERIFICATION.md` for the recorded build, axiom
audit and hashes.

## Axioms

`SnugBridge.erdos_123_powers_2_3_5_snug` depends on `[propext, Classical.choice, Quot.sound]` only.
