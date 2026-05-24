# Morphism / Physics Analogy Scan for Erdős #30

Date: 2026-05-01  
Status: INTERNAL / FINITE-EVIDENCE / NOT A THEOREM CLAIM

## Bottom Line

The strongest current morphism is not "Sidon sets are physics." It is sharper:

> Dense additive extremizers can be treated as finite zero-temperature ground
> faces of a constrained lattice gas, and external observables expose different
> points of that face.

That morphism has operational content because it told us what to measure:
ground-state degeneracy, field response, exact-face handoffs, reset rows,
near-ground shells, and cross-problem controls. It is currently strongest for
Sidon and `B_2[g]`, and it is cleanly absent in the sum-free negative control.

The representation-engineering DNA is the important part. We are not asking
whether physics words can make Sidon sets sound deeper. We are asking whether a
chosen finite state space, observable family, and operator geometry make the
problem compress enough to expose proof targets. Hilbert-space and tensor
language are useful here only as filters on the raw configuration space: they
keep the sectors where invariants, field response, and demotions can be
measured.

## Ranked Candidate Morphisms

| Rank | Candidate morphism | Mathematical object | Physics language | Current verdict | Why it matters |
|---:|---|---|---|---|---|
| 1 | Exact maximizer face as ground-state manifold | Exact Sidon maximizers at fixed `n` | Zero-temperature degenerate ground face | STRONG / FINITE | Corrects the false single-witness view. The object is the whole face, not the first optimizer returned by a scanner. |
| 2 | Observable choice as external field | Prefix, mass, joint observables | Field response / exposed face / normal fan | STRONG / FINITE | Sidon splits `11/11`, `B_2[2]` splits `11/11`, sum-free splits `0/11` in the proof-target packet. |
| 3 | Sidon transfer operator | Legal extension of partial Sidon sets | Hard-core lattice gas / finite-state transfer matrix | STRONG / COMPUTATIONAL | The operator reproduces exact `h(n)` and maximizer counts in checked windows and gives a real partition-function direction. |
| 4 | Handoff at `n = 57` | Birth of multiple exposed candidates in the `h=10` face | Domain wall / field handoff | STRONG / LOCAL | `n=56` is single-exposed, `n=57` splits, `n=58` persists. This is the cleanest local proof target. |
| 5 | Reset at cardinality jumps | `h(n)` increases and ground face collapses | First-order-like reset / nucleation | STRONG / FINITE PATTERN | Seen in Sidon around `25`, `72`, `85`, and in `B_2[g]` jump rows. It suggests proof effort should target construction scarcity at jump rows. |
| 6 | `B_2[g]` as capacity deformation of Sidon | Bounded repeated sum representations | Finite-capacity exclusion law | STRONG / TRANSFER | The field-response behavior survives when "forbidden collision" is relaxed to "bounded collision." That is the best cross-problem support. |
| 7 | Sum-free as rigid control | Maximum sum-free sets | Quiet phase / rigid crystal | STRONG NEGATIVE CONTROL | The transfer API ports, but field response disappears. This demotes the claim that PMF simply narrates every additive problem the same way. |
| 8 | Near-ground shell as thermal bath | `h-1`, `h-2` layers | Low-temperature excitations / entropy bath | USEFUL / CAUTIONARY | Near-ground layers can be huge and smooth, so exact ground-face claims should not be casually extended to capped lower-bound shells. |
| 9 | Inherited-plus-translated persistence | Faces contain previous face and `+1` shift | Propagating branch fan / replication rule | PROMISING / FINITE | Strong in `58..71` and restarts after reset; phi/Fibonacci labels are explicitly demoted. |
| 10 | Maxwell-completed accounting | Full face fixes first-hit artifacts | Missing-information completion | INTERNAL METAPHOR | Good for intuition, but public-safe language should be exact-face field response, not Maxwell/Ramanujan-shadow claims. |

## What Is Load-Bearing

The load-bearing structure is the finite operator:

```text
state = occupied sites + used additive memories + cardinality
move  = skip x, or occupy x if the local exclusion law remains valid
Z_n(beta, fields) = sum_A exp(beta |A| - field_prefix P(A) - field_mass M(A))
```

The zero-temperature limit recovers exact maximizers. Small fields select
exposed points of the ground face. That is the concrete bridge from physics
language to mathematics.

The key point is that this predicts measurements a plain first-hit scan can
miss. In the current packets, the field language changed the question from:

```text
Which maximizer did the scanner print first?
```

to:

```text
How does the whole exact maximizer face respond when different observables are
used as infinitesimal fields?
```

That is why `n=70` was demoted from a real defect to a first-hit artifact, while
`n=57` remained a real local handoff candidate.

## Best Physics Analogies by Usefulness

### 1. Ground-State Face / Normal Fan

This is the best analogy because it is not decorative. It maps directly to
finite data:

- states: exact maximizer witnesses;
- Hamiltonian: `-|A|` plus observable penalties;
- fields: prefix, mass, joint;
- measurable response: which witness minimizes each field;
- entropy: log number of exact maximizers or near-ground states.

This should be the public-safe physics language if any physics language is used.

### 2. Hard-Core Lattice Gas

Sidon is a hard-core exclusion law on additive/difference collisions. The
transfer operator is the real object here. It is not just an analogy; it is an
explicit finite state machine.

### 3. Capacity Deformation

`B_2[g]` is the best cross-problem bridge. Sidon forbids repeated additive
collisions; `B_2[g]` allows bounded collision multiplicity. The persistence of
field response across this deformation is stronger evidence than another Sidon
row.

### 4. Reset / Nucleation

At an `h(n)` jump, the old plateau face can collapse and a tiny new ground face
appears. This behaves like a finite first-order reset, but the safe mathematical
phrase is:

> new-cardinality rows have sparse construction supply before plateau expansion.

### 5. Thermal Bath

Near-ground shells are useful because they tell us where the analogy stops. The
exact face can have sharp field structure while `h-1` and `h-2` layers are huge.
So proofs should target ground-face structure first, not vague thermal
smoothness.

## Demoted or Unsafe Analogies

| Analogy | Status | Reason |
|---|---|---|
| Physics proves Sidon | REJECT | No theorem claim; all current evidence is finite. |
| Phase transition proved | REJECT | Reset/handoff behavior is finite and suggestive, not asymptotic critical behavior. |
| Phi/Fibonacci mediation | DEMOTED | Replication passes, but phi/Fibonacci residual tests do not. |
| Maxwell/Ramanujan shadow public language | INTERNAL ONLY | Useful story language, not a claim boundary. |
| Universal additive-combinatorics field response | REJECT | Sum-free is a clean negative control. |

## Strategic Read

The Collider is now finding a real neighborhood: hard additive exclusion systems
with large exact maximizer faces. The strongest proof lane is not "find one
better witness." It is:

1. formalize exact ground-face field response in Sidon;
2. prove strict non-winner lower bounds for the `56..58` exact-prefix field;
3. convert the `n=57` handoff into a small symbolic face classification;
4. transfer the same lemma shape to `B_2[g]`;
5. keep sum-free as the standing negative control.

The strategic bet is that the right theorem target is a stability/face theorem:

> under a Sidon-like local exclusion law, dense exact extremal faces admit
> multiple field-exposed structures, and jump rows reset the available
> construction supply.

That is still below an Erdős #30 proof. But it is above a pretty analogy: it is
a falsifiable, reproducible map of where proof effort should go next.

## Evidence Pointers

- `PMF_TRANSFER_OPERATOR_ATTACK_2026-04-29.md`
- `MAXWELL_FACE_HANDOFF_LITMUS_2026-04-30.md`
- `B2G_TRANSFER_CROSSPROBLEM_2026-04-29.md`
- `SUMFREE_TRANSFER_CROSSPROBLEM_2026-04-29.md`
- `GROUNDFACE_BRANCH_WINDOW_58_64_2026-05-01.md`
- `GROUNDFACE_RESET_72_2026-05-01.md`
- `proof-targets/FACE_FIELD_RESPONSE_PROOF_TARGET_2026-05-01/`
