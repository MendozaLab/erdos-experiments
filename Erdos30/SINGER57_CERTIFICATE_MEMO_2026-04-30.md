# Singer-Mod-57 Certificate Memo

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_SINGER57_CERTIFICATE / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Exact source packet:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Derived certificate:

```text
EXP-MM-030-PMF-SINGER57-CERTIFICATE-V2-2026-04-30
```

Verification script:

```text
erdos-experiments/Erdos30/scripts/singer57_certificate.py
```

External finite-geometry reference:

```text
https://arxiv.org/html/2502.09536v1
```

The three cited mod-57 perfect difference set representatives checked here are:

```text
PDS_A = {0, 11, 19, 20, 24, 26, 36, 54}
PDS_B = {0, 1, 5, 7, 17, 35, 38, 49}
PDS_C = {16, 19, 30, 38, 39, 43, 45, 55}
```

Each verifies as a size-8 `(57,8,1)` perfect difference set: its ordered
differences cover every nonzero residue mod `57` exactly once.

The V2 certificate also verifies that each cited representative is stabilized
without translation by the multiplier subgroup:

```text
{1, 7, 49}
```

## C1-C5 Result

| check | claim | result |
|---|---|---|
| C1 | all six exported `n=57` witnesses are interval Sidon sets | PASS |
| C2 | all six witnesses share the same positive-difference skeleton | PASS |
| C3 | all six witnesses cover every nonzero residue mod `57` | PASS |
| C4 | mass winner index `3` translates by `+1` to joint winner index `5` | PASS |
| C5 | cited Singer PDS representatives verify, have multiplier stabilizer `{1,7,49}`, and overlap diagnostics are computed | PASS |

The external-facing finite claim is therefore:

> The exported `n=57` exact Sidon maximizers sit on the Singer modulus `57`.
> They are interval Sidon embeddings of size `10`; all six share one
> positive-difference skeleton; each covers every nonzero residue mod `57`;
> and the exposed mass-to-joint handoff is a `+1` translation.

## Boundary Correction

The size-10 witnesses are not literally the cited size-8 Singer PDS objects.

The overlap diagnostic says:

```text
maximum affine point overlap with any cited size-8 PDS representative = 5
contains an affine image of any cited PDS representative = false
```

So the right language is:

```text
Singer-modulus lock-on
```

not:

```text
Singer PDS containment
```

## Core Witnesses

Mass winner:

```text
index 3
[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
```

Joint winner:

```text
index 5
[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
```

Certified relation:

```text
index 5 = index 3 + 1
```

For every exported witness:

```text
mod-57 nonzero residues covered = 56 / 56
multiplicity histogram = {1: 22, 2: 34}
```

## Meaning

This is a stronger and narrower result than the earlier narrative.

The instrument did not find an arbitrary finite wrinkle. It found a
field-sensitive translated embedding at the Singer modulus. But the certificate
also blocks the tempting overclaim: the size-10 witnesses are not Singer
perfect difference sets and do not contain an affine copy of the cited size-8
PDS representatives.

## Claim Boundary

Safe:

- `Singer-modulus finite certificate`
- `size-10 interval Sidon embeddings at modulus 57`
- `full nonzero residue coverage mod 57 with controlled redundancy`
- `translated embedding handoff: index 5 = index 3 + 1`

Unsafe:

- `our witness is a Singer PDS`
- `the witness contains a Singer PDS`
- `PMF proves Sidon`
- `physics solves Erdos #30`
- `SOTA theorem result`

## Next Gate

Translate this certificate into formal proof scaffolding:

```text
1. define the two exposed witnesses;
2. prove interval Sidon by finite difference uniqueness;
3. prove mod-57 nonzero coverage and multiplicity histogram;
4. prove index5 = index3 + 1;
5. prove no affine image of the three cited PDS representatives is contained in
   any of the six exported witnesses.
```

The Lean skeleton is recorded separately as a non-compiled scaffold, not a proof
artifact.

The adjacent `n = 58` branch gate has now been settled:

```text
N58_BRANCH_CERTIFICATE_2026-04-30.md
EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30
```

It confirms the next local row has a second difference skeleton, but that new
skeleton enters on the prefix side. The mass/joint/Pareto candidates remain on
the original translated chain.
