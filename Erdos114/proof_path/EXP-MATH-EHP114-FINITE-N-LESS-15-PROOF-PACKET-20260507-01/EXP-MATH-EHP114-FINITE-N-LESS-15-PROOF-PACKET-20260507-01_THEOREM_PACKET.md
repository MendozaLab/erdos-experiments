# FiniteEHPBelow15 Theorem Packet

## Theorem Statement

`FiniteEHPBelow15`: For every monic complex polynomial `p` of degree `n` with `1 <= n < 15`, the length of the lemniscate `{ z in C : |p(z)| = 1 }` is at most the corresponding length for `p(z)=z^n-1`, under the dependency model listed in this packet.

## Proof Decomposition

1. `n=1`: direct analytic lemma. A monic linear polynomial has a translated unit-circle lemniscate of length `2*pi`, matching `z-1`.
2. `n=2`: literature lemma. The Eremenko-Hayman row is cited as recorded by Erdos Problems #114.
3. `3 <= n <= 14`: machine-certificate lemma. Each degree is covered by a Rust/inari IEEE-1788 interval certificate row from DOI `10.5281/zenodo.19480329` with local SHA verification.
4. Finite union lemma. The rows above exhaust every integer degree satisfying `1 <= n < 15`.

## Audit Summary

- Status: `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`
- Covered degree count: `14`
- Missing degree count: `0`
- Source SHA fail count: `0`
- n=14 atlas cross-check: `PASS`
- n=13 route exception retained: `true`

## Claim Ceiling

This packet proves only the finite range `n < 15` under the stated literature and interval-certificate dependencies. It does not prove Erdos #114 for all degrees and does not bridge Tao's opaque high-degree threshold.
