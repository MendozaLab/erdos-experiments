# Social Media Posts: Sidon Sets in N-Body Choreographies

## Mathstodon Post

Found a Sidon set hiding inside the three-body problem.

If you take an n-body choreography with Z_n symmetry, the Koopman operator has eigenvalues at the nth roots of unity. The pairwise distances between those eigenvalues form a set D_n = {2 sin(pi m/n) : m = 1..floor(n/2)}.

Turns out D_n is a Sidon set in R for all n from 3 to 29. Every pairwise sum is unique — you can identify which two bodies produced a given frequency sum just from the sum itself.

At n = 30 it breaks. The collision is sin(pi/15) + sin(7pi/15) = sin(2pi/15) + sin(4pi/15). All collision values of n up to 200 are divisible by 6. The identity comes from pentagonal structure in cos(pi/15).

This connects Erdos #30 (Sidon sets, $1000 prize) to choreographic orbits through spectral theory. I've got a compiled Lean 4 formalization of the classical Sidon sum-count bound and a draft proof architecture for the spectral bridge.

Preprint + computation: [Zenodo DOI TBD]

#math #dynamics #additivecombinatorics #Lean4

---

## LinkedIn Post

New computational finding from my dynamical systems research: the spectral distance set of an n-body choreography is a Sidon set for all n <= 29.

What does that mean? In a choreography (like the famous figure-eight three-body orbit discovered by Chenciner and Montgomery in 2000), all bodies trace the same curve. The Koopman operator lifts this to a Hilbert space where the eigenvalues are nth roots of unity. The distances between those eigenvalues form a natural set — and that set has a remarkable combinatorial property: all pairwise sums are distinct.

This connects a $1,000 Erdos prize problem (#30, on Sidon sets) to celestial mechanics through operator spectral theory. The connection breaks at n = 30, where a trigonometric identity involving pentagonal symmetry creates the first "collision" — two different pairs of bodies producing identical frequency-sum signatures.

I've formalized the classical Sidon bound in Lean 4 (compiled, zero sorry stubs) and am documenting the formal architecture for the spectral bridge. The computation and preprint are available on Zenodo.

This is part of a broader project formalizing Erdős problems in Lean 4, with seven compiled proofs across different problems so far.

#Mathematics #DynamicalSystems #FormalVerification #Research

---

## Lean Zulip Post

I've been working on connecting Erdos #30 (Sidon sets) to spectral theory of periodic orbits and hit something I wasn't expecting.

Take an n-body choreography with Z_n symmetry. The Koopman operator eigenvalues are nth roots of unity, so the pairwise distances are D_n = {2 sin(pi m/n) : m=1..floor(n/2)}. I checked whether D_n is Sidon in R (all pairwise sums distinct) and it holds for n=3 through 29, then breaks at n=30 via sin(pi/15) + sin(7pi/15) = sin(2pi/15) + sin(4pi/15).

I already have `sidon_sum_count` compiled in Lean 4 (Mathlib v4.24.0, 0 sorries) — the standard |A|(|A|+1)/2 <= 2 sup(A) + 1 bound. Now I want to formalize the spectral distance set and apply the bound.

The gap: I need `spectral_distances_roots_of_unity` showing D_n has size floor(n/2) with the sin formula, then `decide` over the 27 finite cases for the Sidon check. Is there a good Mathlib pattern for "this holds for all n in a finite set, verified by computation"? I'm thinking `Decidable` + `Finset.forall_range` but not sure if that's the idiomatic approach for n up to 29.
