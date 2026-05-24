# Finite-Degree Verification of the Erdos-Herzog-Piranian Lemniscate Conjecture Below Degree 15

## Theorem Statement

`FiniteEHPBelow15`: For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the length for p(z)=z^n-1.

## Public-Facing Abstract

We give a finite-degree verification of the Erdos-Herzog-Piranian lemniscate length conjecture for all degrees 1 <= n < 15, using a combination of elementary reductions, literature input in degree two, and reproducible interval certificates for degrees three through fourteen. This result is independent of, and complementary to, Tao's recent large-degree theorem.

## Relation To Tao's Large-Degree Theorem

Tao proves the conjecture for all sufficiently large degrees; this packet verifies the opposite finite frontier through degree fourteen. A full proof would follow from an effective extraction of Tao's threshold showing that no degrees beyond fourteen remain uncovered, or from additional finite certificates below that extracted threshold.

## Proof Packet Contents

- `n=1`: direct analytic lemma; a monic linear polynomial has a translated unit-circle lemniscate.
- `n=2`: literature row as recorded for the Erdos-Herzog-Piranian problem.
- `3 <= n <= 14`: DOI-backed Rust/inari interval certificate rows with local SHA verification.
- `n=14`: atlas audit included as additional stress-case support, not as a replacement for the degree row.

## Toward a Small-to-Large Degree Completion

Is Tao's sufficiently-large threshold effective in a range that begins at or before n=15? If so, the present finite packet and Tao's theorem would combine into a full proof of Erdos #114.

## Claim Ceiling

Presentation packet only. It presents the finite n<15 theorem and a conditional bridge question; it does not prove the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.

## Local Audit Summary

- Finite packet status: `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`
- Covered degrees: `14`
- Missing degrees: `0`
- Source SHA failures: `0`
- Tao threshold status: `TAO_THRESHOLD_REMAINS_OPAQUE`
- Candidate `N0`: `null`
- `n=15+` computation authorized: `false`
