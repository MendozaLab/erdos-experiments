# EXP-MM-030-PMF-SINGER57-CERTIFICATE-2026-04-30

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_SINGER57_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Derived from `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`.

External finite-geometry comparison uses the three cited mod-57 perfect difference set representatives from `https://arxiv.org/html/2502.09536v1`.

## C1-C5 Checks

- C1 all six witnesses are interval Sidon: `True`
- C2 all six share the same positive-difference skeleton: `True`
- C3 all six cover every nonzero residue mod 57: `True`
- C4 mass winner index 3 translates by +1 to joint winner index 5: `True`
- C5 cited Singer PDS representatives verify and overlap diagnostics computed: `True`

All C checks pass: `True`

## Witness Table

| index | interval Sidon | mod-57 residues covered | mod-57 multiplicity histogram | best affine PDS point overlap |
|---:|---|---:|---|---:|
| 0 | True | 56 | {1: 22, 2: 34} | 5 |
| 1 | True | 56 | {1: 22, 2: 34} | 5 |
| 2 | True | 56 | {1: 22, 2: 34} | 5 |
| 3 | True | 56 | {1: 22, 2: 34} | 5 |
| 4 | True | 56 | {1: 22, 2: 34} | 5 |
| 5 | True | 56 | {1: 22, 2: 34} | 5 |

## Singer PDS Comparison

Maximum affine point overlap between any size-10 witness and any cited size-8 PDS representative: `5`.

Contains any affine image of a cited PDS representative: `False`.

## Interpretation

The n=57 exact interval Sidon embeddings sit on the Singer modulus 57 and cover all nonzero residues mod 57 with controlled redundancy, but they do not literally contain an affine image of the cited size-8 Singer PDS representatives.

## Claim Boundary

Finite certificate only: no asymptotic claim, no theorem proof, and no claim that the size-10 interval embeddings are Singer PDS objects.
