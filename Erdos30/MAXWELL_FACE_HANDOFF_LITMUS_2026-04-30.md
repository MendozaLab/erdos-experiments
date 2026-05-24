# Maxwell Face-Handoff Litmus

Date: 2026-04-30
Problem: Erdos #30
Status: INTERNAL_ONLY / EXACT + INTERPRETIVE
Primary packet: `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`

## Question

Does `n = 57` show a normal-fan / exposed-face split that is absent at
`n = 56`?

## Answer

Yes for the tested finite packet.

The exact ground face at `n = 56` has four maximizers and a single exposed
winner for prefix, mass, joint, and Pareto-minimal selection. The exact ground
face at `n = 57` has six maximizers and splits into two Pareto-minimal exposed
candidates: one mass-selected remnant and one prefix/joint-selected endpoint
candidate.

This supports the internal interpretation:

> `n = 57` is a finite exact-face field handoff.

It does not prove Sidon, does not prove a Mendoza Limit theorem, and does not
make Ramanujan-shadow language public-safe.

## Internal Morphism

The internal analogy is:

```text
first-hit optimizer
-> Ampere-only accounting

full exact maximizer face + field response
-> Maxwell-completed accounting

extra exported face geometry
-> Landauer price paid to preserve continuity across the handoff

residual defect seen by first-hit accounting
-> Ramanujan-shadow-style signal that the visible description is incomplete
```

This is internal summit-climb language. Public language should say:

```text
exact finite computations expose a field-sensitive handoff on the maximizer face
```

## Packet

The first packet with this experiment ID was written before the scanner date
field was corrected. It is superseded for citation by the versioned packet:

- `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`

Required files exist:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.json
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_REPORT.md
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30_RESULTS.sha256
```

SHA256 sidecar verified `OK`.

## Ground-Face Export Summary

| n | h(n) | exact maximizers | export status | exported | field split | Pareto minima |
|---:|---:|---:|---|---:|---|---|
| 56 | 10 | 4 | EXPORTED_ALL | 4 | false | `[3]` |
| 57 | 10 | 6 | EXPORTED_ALL | 6 | true | `[3, 5]` |
| 58 | 10 | 10 | EXPORTED_ALL | 10 | true | `[7, 9]` |

No export cap was hit. The cap was `20`.

## Exposed-Face Winners

| n | prefix winners | mass winners | joint winners | reading |
|---:|---|---|---|---|
| 56 | `[3]` | `[3]` | `[3]` | one exposed point handles all fields |
| 57 | `[4, 5]` | `[3]` | `[5]` | field handoff: mass and joint select different exposed candidates |
| 58 | `[2, 8, 9]` | `[7]` | `[7]` | split persists, with mass and joint aligned |

At `n = 57`, the exact-face geometry is no longer single-faced. The mass field
selects witness `3`; the joint field selects witness `5`; prefix is tied on
`4` and `5`.

## Key Witnesses at n = 57

Mass-selected remnant:

```text
index 3
[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]
prefix_ratio = 0.029081011040818568
mass_ratio   = 0.009629685024932394
joint_score  = 0.03871069606575096
```

Prefix/joint-selected endpoint candidate:

```text
index 5
[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]
prefix_ratio = 3.0994951715419028e-16
mass_ratio   = 0.02888905507479718
joint_score  = 0.028889055074797488
```

The symmetric-difference distance between these two exposed candidates is `18`.
That is a genuine face-to-face jump, not a minor local perturbation.

## Interpretation

The Maxwell Face-Handoff Litmus passes at finite evidence level:

```text
n=56:
  one exposed point
  first-hit accounting and face-aware accounting agree

n=57:
  two Pareto-minimal exposed candidates
  different fields select different ground states
  first-hit accounting is incomplete

n=58:
  split persists after the handoff
```

This is the strongest current internal formulation:

> the missing information is the full exact maximizer face and its normal-fan
> response to fields.

The physics analogy is not decoration. It tells us what to measure: exposed
ground-state faces under infinitesimal fields.

## Claim Boundary

Safe internal language:

- `Maxwell Face-Handoff Litmus`
- `exact-face field handoff`
- `zero-temperature field selection`
- `missing-information completion`
- `first-hit accounting is incomplete`

Unsafe public language:

- `Mendoza Limit proves Sidon`
- `Ramanujan shadow proves the handoff`
- `physics solves Erdos #30`
- `phase transition proved`
- `SOTA theorem result`

## Next Gate

The symbolic classification of the `n = 57` six-point face now lives in:

```text
N57_FACE_CLASSIFICATION_2026-04-30.md
```

The finite lemma-style statement is:

> At `n = 57`, the exact `h = 10` face decomposes into four inherited remnant
> maximizers and two endpoint-shift maximizers. The exposed Pareto face has two
> candidates: a remnant mass winner and an endpoint-shift joint winner. At
> `n = 56`, the exported exact face has a single exposed candidate selected by
> all three tracked fields.

The certificate-style verification now exists:

```text
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-30/EXP-MM-030-PMF-N57-FACE-CERTIFICATE-2026-04-30_RESULTS.json
```

It verifies the six exported witnesses are Sidon, recomputes the field-winner
and Pareto roles, and confirms the `4` remnant / `2` endpoint-shift split. The
certificate inherits exhaustiveness of the six `h = 10` maximizers from the
source exact packet rather than reproving it.

The wider `56 -> 57 -> 58` window certificate now exists:

```text
FACE_HANDOFF_WINDOW_CERTIFICATE_56_58_2026-04-30.md
EXP-MM-030-PMF-FACE-HANDOFF-WINDOW-CERTIFICATE-56-58-2026-04-30
```

It verifies the whole handoff pattern: one exposed face at `56`, split appears
at `57`, and split persists at `58`.
