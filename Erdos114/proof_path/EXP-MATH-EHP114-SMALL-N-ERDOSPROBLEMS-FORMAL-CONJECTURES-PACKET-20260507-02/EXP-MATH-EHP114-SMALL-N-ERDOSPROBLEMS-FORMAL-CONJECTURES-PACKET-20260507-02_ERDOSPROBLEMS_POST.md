# Draft erdosproblems.com post for Problem 114

I would like to record a finite-degree verification related to the Erdos-Herzog-Piranian lemniscate conjecture.

**Finite statement.** For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the length for p(z)=z^n-1.

The dependency split is deliberately finite:

1. `n=1`: elementary, since a monic linear polynomial has a translated unit circle as its unit lemniscate.
2. `n=2`: the Eremenko-Hayman case.
3. `3 <= n <= 14`: reproducible interval-certificate rows, archived under DOI `10.5281/zenodo.19480329`, with local SHA checks in the accompanying packet.

The current packet status is `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`: it records `14` covered positive degrees below `15`, `0` missing degree rows, and `0` source SHA failures.

This is not meant to assert the full conjecture. Rather, it is the small-degree frontier complementary to Tao's sufficiently-large-`n` theorem. The remaining bridge question is whether Tao's threshold can be made effective at or below `n=15`; if not, the intervening finite degrees would still need separate certification.

I would appreciate guidance on whether this is the right form in which to record the finite side on the problem page, and whether the formal-conjectures entry should name the theorem as a finite variant rather than changing the status of the main conjecture.

Certificate packet: `EXP-MATH-EHP114-FINITE-N-LESS-15-PROOF-PACKET-20260507-01`
Claim ceiling: finite `n < 15` result only; conditional Tao-threshold bridge only.

Disclosure: the local packet was assembled with AI-assisted tooling, but the submitted claim is only the SHA-checked finite certificate surface described above.
