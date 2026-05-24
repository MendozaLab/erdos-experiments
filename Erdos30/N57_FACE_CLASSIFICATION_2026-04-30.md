# n = 57 Exact-Face Classification

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_FROM_EXACT_PACKET / INTERPRETIVE

## Source Packet

Primary packet:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.json
```

Experiment ID:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

Required packet files:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.json
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_REPORT.md
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.sha256
```

SHA256 verification was rerun from the packet directory and returned `OK`.
The sidecar uses a relative filename, so it must be checked from:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30
```

## Question

What is the exact structure of the six-point maximizer face at `n = 57`?

The Maxwell note already shows that `n = 57` is the first tested point where
the exposed face splits. This note classifies the six exact `h = 10`
maximizers into inherited remnants, endpoint-shift states, and exposed
field-selected candidates.

## Answer

At `n = 57`, the exact maximizer face has:

- four inherited remnant states from `n = 56`;
- two new endpoint-shift states containing the newly available endpoint `57`;
- exactly two Pareto-minimal exposed candidates;
- one remnant selected by the mass field;
- one endpoint-shift state selected by the joint field.

That is the finite face-to-face handoff:

```text
mass field  -> inherited remnant face
joint field -> endpoint-shift face
```

This is not a theorem about the Sidon asymptotic. It is a packet-backed finite
classification of the exact `n = 57` maximizer face.

## Classification Table

| index | class | inherited from n=56 | persists in n=58 | contains 0 | contains 57 | prefix winner | mass winner | joint winner | Pareto minimal | joint score |
|---:|---|---:|---:|---|---|---|---|---|---|---:|
| 0 | remnant | 0 | 0 | yes | no | no | no | no | no | 0.306607895725 |
| 1 | remnant | 1 | 1 | yes | no | no | no | no | no | 0.106310447206 |
| 2 | remnant | 2 | 4 | no | no | no | no | no | no | 0.239008144584 |
| 3 | remnant | 3 | 5 | no | no | no | yes | no | yes | 0.038710696066 |
| 4 | endpoint-shift | none | 6 | no | yes | yes | no | no | no | 0.171408393444 |
| 5 | endpoint-shift | none | 7 | no | yes | yes | no | yes | yes | 0.028889055075 |

The two Pareto-minimal exposed candidates are indices `3` and `5`.

Index `3` is the inherited remnant selected by the mass field:

```text
[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
prefix_ratio = 0.029081011040818568
mass_ratio   = 0.009629685024932394
joint_score  = 0.03871069606575096
```

Index `5` is the endpoint-shift state selected by the joint field:

```text
[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
prefix_ratio = 3.0994951715419028e-16
mass_ratio   = 0.02888905507479718
joint_score  = 0.028889055074797488
```

The symmetric-difference distance between these two exposed candidates is `18`.
That is a large jump on a six-point face, not a one-coordinate local edit.

## Finite Lemma Candidate

For the exported exact maximizer face in the verified packet:

```text
At n = 56, all tracked fields select the same Pareto-minimal maximizer.
At n = 57, the exact h = 10 face decomposes into four inherited remnant
maximizers and two endpoint-shift maximizers. The exposed Pareto face has two
candidates: a remnant mass winner and an endpoint-shift joint winner.
At n = 58, both n = 57 exposed candidates persist, and the split remains.
```

This is the cleanest current statement of the Maxwell Face-Handoff Litmus.
The missing information in the first-hit description is exactly the exposed
normal-fan structure of the full maximizer face.

## Claim Boundary

Safe language:

- `finite exact-face handoff`
- `field-selected exposed candidates`
- `remnant-to-endpoint face switch`
- `first-hit accounting loses face geometry`

Unsafe language:

- `PMF proves Sidon`
- `physics solves Erdos #30`
- `Mendoza Limit proves the handoff`
- `Ramanujan shadow proves the defect`
- `SOTA theorem result`

## Next Gate

A certificate-style derivation from the exported face now exists:

```text
EXP-MM-030-PMF-N57-FACE-CERTIFICATE-2026-04-30
```

It verifies:

- all six exported `n = 57` witnesses are Sidon;
- the recomputed prefix, mass, joint, and Pareto roles match the source packet;
- the face decomposes into four remnants and two endpoint-shift states;
- exposed indices `3` and `5` have symmetric-difference distance `18`.

The only piece not reproved inside the derived certificate is exhaustiveness of
the six `h = 10` witnesses; that remains inherited from the source exact packet.
