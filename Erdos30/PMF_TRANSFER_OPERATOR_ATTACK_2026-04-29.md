# PMF Transfer-Operator Attack for Erdos #30

Date: 2026-04-29
Lane: Prima Materia Framework / Erdős Collider
Problem: Erdős #30, finite Sidon set rigidity
Status: attack design, evidence-bound to exact Rust packets

## The Plain Story

We started with what looked like a finite compatibility ridge.

On the exact Sidon maximizer surface, the mass observable and the prefix
observable kept agreeing more often than they had any obvious right to agree.
The direct first-hit scan over `51 <= n <= 70` held in `18/20` values across
`298,968` exact maximizers. The apparent failures were `n = 57` and `n = 70`.

At first that looked like a punctuated ridge: mostly stable, but with two
pinches.

Then the top-k frontier scan changed the story. It showed that `n = 70` was not
a true mathematical obstruction. The exact maximizer face at `n = 70` contains
joint witnesses with essentially zero prefix cost and zero density-adjusted
mass cost. The old failure was a scanner-selection artifact: one tied witness
was reported first, but the face itself had the compatible witness.

That left `n = 57` as the real local event. It is not a collapse; it is a
handoff. Across `56, 57, 58`, neighboring witnesses exchange roles on the
frontier. The face still has structure, but the optimizer label changes.

Then `n = 71` rebounded. It scanned `203,840` exact maximizers, visited
`2,185,808,767` search nodes, and passed the first-hit comparison. Its best
top-k joint witness has near-zero prefix ratio and mass ratio
`0.0014239337207513653`.

That is the discovery arc:

1. Exact enumeration found a ridge.
2. First-hit reporting made the ridge look punctuated at `57` and `70`.
3. Top-k frontier witnesses corrected `70`.
4. The remaining signal is a real local handoff around `57`.
5. The right object is not one optimizer; it is the whole exact maximizer face.

## Atlas-DNA Claim Boundary

The physics language is not the proof. It is the microscope.

The working thesis is MDL first, physics second, proof last-mile. Physics is
useful here because it is a library of low-description-length structures that
nature has already stress-tested: locality, conservation, packing, stability,
resonance, entropy, phase transition, degeneracy, and field response. When an
Erdős problem shows the same compression signature, the Atlas move is not to
declare a physical explanation. The move is to build an explicit finite
operator and ask whether exact enumeration, parity checks, frontier witnesses,
and near-ground diagnostics confirm a real mathematical structure.

So the safe sentence is:

> The analogy generated the operator; exact computation tests it.

The unsafe sentence remains:

> Physics proves Sidon.

## Collider Translation

The collider move is to stop treating this as a bag of Sidon sets and start
treating it as a finite lattice gas.

| Collider object | Sidon / PMF meaning |
|---|---|
| Lattice sites | integers `0, 1, ..., n` |
| Particle | selected element of a Sidon set |
| Hard-core exclusion | no repeated positive pairwise difference |
| Configuration | subset `A subset [0,n]` |
| Energy | deficit from maximal cardinality, plus observable penalties |
| Ground states | exact maximizers `|A| = h(n)` |
| Low-energy face | top-k frontier witnesses on the exact maximizer surface |
| Defect / domain wall | real handoff pinch such as `n = 57` |
| Degeneracy / tied face | `n = 70`, where first-hit selection lied |
| Perfect crystal | Singer/projective-plane configurations |
| Transfer operator | finite-state propagation of admissible configurations |

In this language, the exact Rust scanner is already doing ground-state
enumeration. The PMF transfer-operator attack asks for the next layer:

> What does the operator spectrum of the Sidon lattice gas say about ground
> state degeneracy, handoff pinches, and the correction term above `sqrt n`?

## Transfer-Operator Form

A transfer state should remember only the information needed to extend a
partial Sidon configuration.

For a prefix scan ending at position `x`, a state can be represented as:

- `occupied_suffix`: recent selected sites, or enough selected sites to update
  differences;
- `used_differences`: bitset of positive differences already realized;
- `cardinality`: number of occupied sites;
- optional observable accumulators: prefix residual, mass, density-adjusted
  mass, and joint score.

Adding a new site `x` is legal iff all new differences `x - a` for `a in A`
are unused. The transfer operator has two moves:

```text
skip x:   state -> state
occupy x: state -> state' if all new differences are unused
```

The partition function is:

```text
Z_n(beta, fields) = sum_A exp(beta * |A| - field_prefix * P(A) - field_mass * M(A))
```

The zero-temperature limit recovers exact maximizers. Small fields tilt the
ground-state face toward the prefix, mass, or joint frontier.

This is why the top-k witnesses matter. They are not just examples; they are
finite samples of how the ground-state face responds to external fields.

## Evidence Bound to Current Packets

Primary packets:

- `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`
- `EXP-MM-030-RUST-TOPK-FRONTIER-71-2026-04-29`
- `EXP-MM-030-RUST-TOPK-FRONTIER-57-2026-04-29`
- `EXP-MM-030-RUST-TOPK-FRONTIER-70-2026-04-29`
- `EXP-MM-030-EXACT-REFERENCE-10-30-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PARITY-10-30-V2-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PARITY-10-30-V3-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-FIELD-FRONTIER-56-58-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-69-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-70-2026-04-29`
- `EXP-MM-030-PMF-TRANSFER-PRUNED-FRONTIER-71-2026-04-29`

The `51 <= n <= 70` packet:

- exact maximizers scanned: `298,968`;
- runtime: `284.640270375` seconds;
- first-hit direct compatibility: `18/20`;
- first-hit failures: `57, 70`;
- top-k correction: `70` is recovered at face level.

The joint frontier over `51 <= n <= 70` has exact or near-exact zero-joint
witnesses at:

```text
51, 52, 53, 54, 60, 62, 64, 66, 68, 70
```

Nonzero joint-best scores occur at:

```text
55, 56, 57, 58, 59, 61, 63, 65, 67, 69
```

The largest joint score in that window occurs at `n = 58`:

- joint score: `0.03616352164892123`;
- prefix ratio: `0.028641813045612422`;
- mass ratio: `0.007521708603308812`;
- witness: `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]`.

This is why we should not overclaim "perfect joint witnesses always exist." The
honest claim is softer and more useful:

> exact maximizer faces carry a small-joint frontier, with isolated handoff
> events that look like domain walls in the lattice gas.

## Transfer-Operator Parity Result

The first Rust transfer-operator crate now exists at:

```text
erdos-experiments/Erdos30/rust-transfer-operator/
```

Its state is explicit and bitset-backed:

```text
State = {
  occupied_mask: u128,
  used_differences_mask: u128,
  cardinality: u8
}
```

Each reported row also serializes a representative ground-state encoding:
occupied sites, the exposed occupied suffix, used positive-difference memory,
and cardinality. The transition rule is the lattice-gas rule:

```text
skip x:   state -> state
occupy x: state -> state' iff every new difference |x-a| is unused
```

The first parity packet,
`EXP-MM-030-PMF-TRANSFER-PARITY-10-30-2026-04-29`, deliberately remains in the
artifact trail as a useful failed turn. It used shifted positive-domain
reference counts against a `0..=n` operator. That was not a math failure; it was
a convention mismatch.

The corrected packet,
`EXP-MM-030-PMF-TRANSFER-PARITY-10-30-V2-2026-04-29`, binds against the current
exact Rust reference packet `EXP-MM-030-EXACT-REFERENCE-10-30-2026-04-29`.
The current citation target is
`EXP-MM-030-PMF-TRANSFER-PARITY-10-30-V3-2026-04-29`, which keeps the same
ground-state parity and also uses the exact scanner's raw-prefix frontier
observable rather than the earlier recentered residual. Result:

```text
h(n) parity:              21 / 21
maximizer-count parity:   21 / 21
mismatches:               []
checksum:                 OK
runtime:                  0.10552425 sec
```

This establishes a narrow but important fact: the PMF state representation
reproduces the exact ground-state surface for the trusted small window. It does
not yet establish asymptotic theorem language, and it does not yet prove that
the top-k observable definitions are identical to the older maximizer scanner.

## Field-Frontier Parity Result

The next packet,
`EXP-MM-030-PMF-TRANSFER-FIELD-FRONTIER-56-58-2026-04-29`, answers the
observable question directly. After cardinality selects the exact `h(n)`
ground-state layer, small prefix and mass fields select the same top witnesses
as the Rust exact maximizer scanner for the suspected defect window.

Machine comparison against `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`
for the top-1 prefix, mass, and joint witnesses over `56 <= n <= 58` produced
an empty diff.

| n | prefix-field winner | mass-field winner | joint-field winner | joint score |
|---|---|---|---|---:|
| 56 | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | 0.011840 |
| 57 | `[2, 3, 8, 12, 25, 28, 36, 43, 55, 57]` | `[1, 3, 15, 22, 30, 33, 46, 50, 55, 56]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | 0.028889 |
| 58 | `[0, 2, 15, 21, 22, 32, 46, 50, 55, 58]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | `[2, 4, 16, 23, 31, 34, 47, 51, 56, 57]` | 0.036164 |

This is the current best evidence that `57` is a field-sensitive handoff. The
cardinality spectrum is smooth, but the zero-temperature field response splits
the face into different prefix, mass, and joint winners.

## Pruned 69-71 Scale Gate

The scale report is
`PRUNED_TRANSFER_FRONTIER_69_71_2026-04-29.md`.

Question:

> Can the transfer operator reproduce the 69-71 top-k frontier behavior without
> enumerating every low-cardinality Sidon state?

Answer: yes, for the tested top-1 prefix, mass, and joint frontier witnesses.

The full HashMap transfer operator was too blunt at this scale. The working
version uses the same transfer state but switches to a layer-pruned DFS:

```text
keep only branches that can still reach h(n)-d
```

For the 69-71 scale gate, `d = 0`, so only branches capable of reaching the
exact ground cardinality `h(n) = 10` survive.

| n | h(n) | retained terminal states | pruned states | runtime sec | joint score |
|---|---:|---:|---:|---:|---:|
| 69 | 10 | 66,412 | 218,433,851 | 28.298015 | 0.001481 |
| 70 | 10 | 117,202 | 260,934,836 | 43.711637 | ~0 |
| 71 | 10 | 203,840 | 311,297,246 | 53.424933 | 0.001424 |

Machine comparison against the exact Rust maximizer packets produced an empty
diff for the top-1 prefix, mass, and joint witnesses over `69 <= n <= 71`.

This closes the immediate scale blocker for the PMF lane. It does not prove the
Sidon asymptotic; it says the transfer-operator state machine can reproduce the
known exact maximizer face behavior at the current frontier windows without
enumerating low-cardinality terminal states.

## First Spectral Defect Probe

The near-ground D1 report is
`PRUNED_TRANSFER_NEARGROUND_D1_69_71_2026-04-29.md`.

It extends the 69-71 pruned scale gate from `d = 0` to `d = 1`, retaining both
the exact ground face and the first near-ground layer. Frontier parity still
holds: the top-1 prefix, mass, and joint winners remain identical to the exact
Rust maximizer packets.

| n | ground states | first-excited states | terminal retained | runtime sec |
|---|---:|---:|---:|---:|
| 69 | 66,412 | 16,623,922 | 16,690,334 | 35.101241 |
| 70 | 117,202 | 22,801,688 | 22,918,890 | 38.715517 |
| 71 | 203,840 | 31,058,406 | 31,262,246 | 46.397141 |

The first-excited layer is huge and smoothly increasing. That does not weaken
the handoff story; it narrows it. The visible PMF signal is not a raw
cardinality-layer spike. It is the zero-temperature field response on the exact
maximizer face.

The D2 single-point follow-up is
`PRUNED_TRANSFER_NEARGROUND_D2_70_2026-04-30.md`.

It extends the `n = 70` probe from `d = 1` to `d = 2`, retaining the exact
ground face plus the first two near-ground layers:

| n | ground states | h-1 states | h-2 states | terminal retained | runtime sec | joint score |
|---|---:|---:|---:|---:|---:|---:|
| 70 | 117,202 | 22,801,688 | 150,565,250 | 173,484,140 | 95.237988 | ~0 |

The D2 packet preserves exact top-1 frontier parity with the Rust top-k
reference for prefix, mass, and joint fields. The second near-ground layer is
large, but it does not turn `70` into a visible spectral/degeneracy defect. The
best reading remains that `70` was a first-hit / field-selection artifact on a
degenerate exact maximizer face.

The D2 window follow-up is
`PRUNED_TRANSFER_NEARGROUND_D2_69_71_2026-04-30.md`.

It reruns the same `d = 2` retention rule uniformly over `69 <= n <= 71`.
Parity passes for all three rows: `h(n)` and maximizer counts match the exact
reference. The second near-ground layer remains smooth across the window:

| n | ground states | h-1 states | h-2 states | `(h-1)/h` | `(h-2)/(h-1)` | joint score |
|---|---:|---:|---:|---:|---:|---:|
| 69 | 66,412 | 16,623,922 | 122,996,090 | 250.315033 | 7.398741 | 0.001481 |
| 70 | 117,202 | 22,801,688 | 150,565,250 | 194.550332 | 6.603250 | ~0 |
| 71 | 203,840 | 31,058,406 | 183,623,966 | 152.366591 | 5.912215 | 0.001424 |

This sharpens the interpretation: `70` is special in the zero-temperature
field-selected joint score, not in the raw near-ground state bath. The
near-ground ratios move smoothly downward across the window, so the D2 evidence
does not support a raw spectral singularity at `70`.

The face-level witness-distance report is
`FIELD_HANDOFF_WITNESS_TABLE_51_71_2026-04-30.md`.

It compares the top-1 prefix, mass, and joint witnesses across `51 <= n <= 71`
by symmetric-difference distance. This separates the two first-hit failures:

| n | p-m d | p-j d | m-j d | joint score | reading |
|---:|---:|---:|---:|---:|---|
| 57 | 16 | 16 | 18 | 0.02889 | wide three-way handoff, no zero-joint rescue |
| 70 | 14 | 12 | 14 | ~0 | first-hit artifact; zero-joint witness exists |

That makes `57` the stronger field/optimizer handoff candidate. The next #30
gate should focus on the `h = 10` birth window, especially `56 <= n <= 59`,
rather than continuing to count deeper layers at `69 <= n <= 71`.

The D2 birth-window follow-up is
`PRUNED_TRANSFER_NEARGROUND_D2_56_59_2026-04-30.md`.

It required adding exact reference rows for `56 <= n <= 59` to the Rust
transfer parity table from `EXP-MM-030-RUST-TOPK-FRONTIER-51-70-2026-04-29`.
After rebuild, the D2 packet passed parity for all four rows.

| n | ground states | h-1 states | h-2 states | `(h-1)/h` | joint score | reading |
|---:|---:|---:|---:|---:|---:|---|
| 56 | 4 | 69,564 | 4,813,066 | 17,391.000000 | 0.011840 | rigid top-1 field response |
| 57 | 6 | 120,704 | 6,511,012 | 20,117.333333 | 0.028889 | field-handoff onset |
| 58 | 10 | 200,946 | 8,684,372 | 20,094.600000 | 0.036164 | split persists |
| 59 | 18 | 330,056 | 11,500,722 | 18,336.444444 | 0.016531 | three-way split returns |

This makes the handoff story sharper: the near-ground bath is already enormous
at `56`, but the top-1 field response is still rigid there. The qualitative
change at `57` is on the exact maximizer face, not in the mere existence of a
large near-ground layer.

The internal Maxwell-completion follow-up is
`MAXWELL_FACE_HANDOFF_LITMUS_2026-04-30.md`.

It adds a small exact-face export mode to the Rust transfer scanner and runs:

```text
EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30
```

The first non-versioned packet was superseded because the scanner date field
was still hardcoded to `2026-04-29`; the `V2` packet has `date = 2026-04-30`
and is the citation target. SHA verification passed.

| n | exact maximizers | export status | field split | prefix winners | mass winners | joint winners | Pareto minima |
|---:|---:|---|---|---|---|---|---|
| 56 | 4 | EXPORTED_ALL | false | `[3]` | `[3]` | `[3]` | `[3]` |
| 57 | 6 | EXPORTED_ALL | true | `[4,5]` | `[3]` | `[5]` | `[3,5]` |
| 58 | 10 | EXPORTED_ALL | true | `[2,8,9]` | `[7]` | `[7]` | `[7,9]` |

This answers the litmus question: `n = 57` shows a normal-fan / exposed-face
split absent at `n = 56`. Internally, this supports the Maxwell Face-Handoff
language: first-hit accounting is the incomplete Ampere-only view, while the
exported exact face plus field response is the completed accounting. Public
language remains blocked at finite exact-face field handoff.

The derived classification note is
`N57_FACE_CLASSIFICATION_2026-04-30.md`. It sharpens the handoff into a finite
lemma candidate: the `n = 57` exact `h = 10` face consists of four inherited
remnant maximizers and two endpoint-shift maximizers. The exposed Pareto face
has exactly two candidates in the exported packet: a remnant mass winner and an
endpoint-shift joint winner. Their symmetric-difference distance is `18`, so
the handoff is a genuine face-to-face switch in the finite geometry, not a
minor local edit.

The derived certificate packet is
`EXP-MM-030-PMF-N57-FACE-CERTIFICATE-2026-04-30`. It directly verifies the six
exported witnesses are Sidon, recomputes the field-winner and Pareto roles from
the packet coordinates, and confirms the remnant/endpoint split. It inherits
the no-other-`h=10` exhaustiveness claim from the source exact packet.

The window certificate is
`EXP-MM-030-PMF-FACE-HANDOFF-WINDOW-CERTIFICATE-56-58-2026-04-30`, with the
readable note `FACE_HANDOFF_WINDOW_CERTIFICATE_56_58_2026-04-30.md`. It verifies
the full local pattern: `n = 56` has a single exposed face, `n = 57` introduces
a remnant-vs-endpoint split, and `n = 58` preserves that split. All exported
witnesses in the three-row window pass direct Sidon checks, and all
field/Pareto roles recompute from packet coordinates.

The difference-skeleton scout is
`EXP-MM-030-PMF-DIFFSET-SKELETON-WINDOW-56-58-2026-04-30`, with the readable
note `DIFFSET_SKELETON_WINDOW_56_58_2026-04-30.md`. It is an important
correction: `n = 57` is not a difference-set skeleton split. All six exported
maximizers at `57` share one positive-difference set. The handoff is therefore
an embedding/field-selection switch inside one skeleton. A second local
difference skeleton first appears at `n = 58`.

The affine embedding certificate is
`EXP-MM-030-PMF-AFFINE-EMBEDDING-N57-2026-04-30`, with the readable note
`AFFINE_EMBEDDING_HANDOFF_N57_2026-04-30.md`. It gives the cleanest current
finite mechanism: the `n = 57` mass winner is index `3`, the joint winner is
index `5`, and `index 5 = index 3 + 1`. Thus the handoff is a translated
embedding selected by the field response, not a replacement of the underlying
difference skeleton.

The Singer-mod-57 certificate is
`EXP-MM-030-PMF-SINGER57-CERTIFICATE-V2-2026-04-30`, with the readable note
`SINGER57_CERTIFICATE_MEMO_2026-04-30.md`. It verifies five elementary claims:
all six `n=57` exported witnesses are interval Sidon; all six share one
positive-difference skeleton; all six cover every nonzero residue mod `57`;
the mass-to-joint exposed handoff is the `+1` translation `index 5 = index 3 +
1`; and the three cited size-8 Singer PDS representatives verify as
`(57,8,1)` perfect difference sets with multiplier stabilizer `{1,7,49}`. The
overlap diagnostic is boundary-setting: the maximum affine point overlap with
any cited PDS representative is `5/8`, so the witnesses lock onto the Singer
modulus but do not contain an affine copy of the cited Singer PDS
representatives.

The `n = 58` branch certificate is
`EXP-MM-030-PMF-N58-BRANCH-CERTIFICATE-2026-04-30`, with the readable note
`N58_BRANCH_CERTIFICATE_2026-04-30.md`. It settles the adjacent branch question:
`n = 58` is the first local row with two difference skeletons, but the new
skeleton is prefix-side only. The Pareto candidates `7` and `9` remain on the
original translated skeleton, with `index 9 = index 7 + 1`; `index 7` is the
persisted `n = 57` joint winner. The local mechanism is now:

```text
56: single exposed embedding
57: translated-embedding handoff inside one Singer-modulus skeleton
58: first local skeleton branch appears, but Pareto remains on translated chain
```

## Sum-Free Cross-Problem Gate

The first changed-rule report is
`SUMFREE_TRANSFER_CROSSPROBLEM_2026-04-29.md`.

Question:

> Is the field-sensitive handoff a Sidon-only artifact, or does the same PMF
> state-machine signature appear under a changed local exclusion rule?

Answer:

> The transfer-state API ports cleanly, but the field-sensitive handoff
> signature does not appear in the first sum-free ground-face test.

The new binary is:

```text
erdos-experiments/Erdos30/rust-transfer-operator/src/bin/sumfree_transfer.rs
```

The local rule changes from Sidon difference memory to:

```text
occupy x iff no occupied a,b satisfy a + b = x
```

Packets:

- `EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-20-40-2026-04-29`
- `EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-41-60-2026-04-29`

Combined result:

```text
20 <= n <= 60
h(n)=ceil(n/2) formula parity: 41 / 41
frontier split count: 0
```

The sum-free face is rigid in this finite model: prefix, mass, and joint fields
all select the same odd-set witness throughout the tested window. This is a good
negative control. The API generalizes, but the Sidon handoff signature is not a
generic artifact of the machinery.

The first full-state spectral probe is
`EXP-MM-030-PMF-TRANSFER-DEFECT-56-58-2026-04-29`. It enumerates every exact
Sidon transfer state, not just maximizers, and extracts cardinality-layer
observables around the suspected `57` pinch.

| n | h(n) | ground degeneracy | entropy ln | h-1 states | h-2 states | gap | terminal states | best joint score |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 56 | 10 | 4 | 1.386294 | 69,564 | 4,813,066 | 1 | 33,911,490 | 0.011840 |
| 57 | 10 | 6 | 1.791759 | 120,704 | 6,511,012 | 1 | 40,821,696 | 0.009630 |
| 58 | 10 | 10 | 2.302585 | 200,946 | 8,684,372 | 1 | 48,993,329 | 0.007522 |

The Feynman interpretation changed again. `57` does not currently look like a
raw cardinality-spectrum singularity: the gap stays `1`, ground degeneracy
rises smoothly from `4 -> 6 -> 10`, and near-ground counts rise smoothly too.
So the sharper reading is:

> `57` is more likely a field/optimizer handoff on the exact maximizer face
> than a standalone spectral-gap defect.

That is still PMF-useful. It says the transfer operator should next be tilted
by prefix and mass fields, because the interesting event is not just the
zero-temperature cardinality layer. It is how the ground-state face responds
when the observables pull in different directions.

## The Discovery We Can Tell

The story is not "we solved Sidon with physics."

The stronger and more credible story is:

> We built an exact finite collider for the Sidon problem. It enumerates every
> maximizer, measures competing observables, and then asks whether the apparent
> failures are real defects or artifacts of optimizer selection. The first
> surprise was that the ridge held across hundreds of thousands of exact
> maximizers. The second surprise was more important: one apparent failure
> vanished when we looked at the whole face instead of the first witness. That
> forced the right mathematical object into view: the exact maximizer face,
> not a single selected optimizer.

This story is defensible because it includes the twist. We did not preserve our
first interpretation. We let the top-k frontier falsify it.

That is the Feynman part: the machine found an effect, then the next machine
found our mistake, and the remaining signal got sharper.

## Attack Plan

### Step 1: Transfer-Operator Parity

Build a Rust transfer-operator scanner that reproduces trusted exact outputs:

- `h(n)`;
- maximizer counts;
- top-k joint frontier witnesses;
- first-hit versus face-aware diagnostics.

Parity windows:

- `10 <= n <= 30` against the older Python exact scan;
- `51 <= n <= 60` against the Rust exact maximizer packets;
- single-point parity at `70` and `71`.

The first acceptance gate is not speed or theory. It is cold parity with the
existing exact enumerator.

### Step 2: Field-Tilted Ground States

Add three external fields:

- cardinality field, favoring large `|A|`;
- prefix field, favoring low prefix residual;
- mass field, favoring low density-adjusted mass deviation.

At large cardinality field, the operator should recover exact maximizer faces.
With small prefix/mass tilts, it should reproduce the top-k frontier behavior.

This turns the frontier scan into a thermodynamic probe.

### Step 3: Defect Diagnostics

Measure whether `n = 57` has a spectral signature:

- degeneracy of the ground-state face;
- gap between ground and first excited layer;
- sensitivity of the selected witness to prefix/mass field tilts;
- number of near-ground states sharing the same difference skeleton;
- persistence of the `56 -> 57 -> 58` witness handoff.

If `57` is a real lattice-gas domain wall, these diagnostics should spike or
change slope locally. If it is just another tied face, they should flatten like
`70`.

### Step 4: Cross-Problem PMF Test

Run the same transfer-operator shape on:

- #755 `B_h[g]`, where the exclusion rule changes from no repeated differences
  to bounded repeated sums;
- #166 sum-free sets, where the hard-core rule is additive closure exclusion.

This is the morphism test. If the PMF lane is real, the same operator
vocabulary should transfer across these problems with only the local exclusion
rule changed.

## Claim Boundary

Safe internal claim:

> The PMF lane gives a concrete transfer-operator attack API for #30. Existing
> exact packets already show ground-state face structure and at least one
> field-sensitive handoff candidate. The first transfer-operator parity packet
> reproduces exact `h(n)` and maximizer counts on `10 <= n <= 30`; the
> field-frontier packet reproduces top-1 prefix/mass/joint winners on
> `56 <= n <= 58`; the pruned scale packet reproduces top-1
> prefix/mass/joint winners on `69 <= n <= 71`; the D1 near-ground probe shows
> the first-excited layer is large and smooth, so the event is not a raw
> cardinality-layer singularity; the D2 `69 <= n <= 71` window shows the same
> smooth near-ground bath one layer deeper; the D2 `56 <= n <= 59` birth-window
> probe points to `57` as the stronger exact-face handoff onset.

Unsafe public claim:

> PMF proves the Sidon asymptotic or solves #30.

Current status:

`ACTIVE`, not `VALIDATED`. The transfer-operator implementation has passed the
small exact ground-state parity gate, the first field-frontier parity gate, and
the 69-71 pruned scale gate. The `d = 1` near-ground layer has also been run on
69-71 and supports the narrower field-response interpretation. The `d = 2`
near-ground window has now been run on 69-71 as well; it preserves parity and
shows smooth near-ground ratios, with `70` special only as a field-selected
joint-score dip. The `56 <= n <= 59` D2 birth-window probe has now passed as
well and points to `57` as the stronger exact-face handoff onset. The #166
sum-free port passes as an API generalization but behaves as a rigid negative
control, not a second handoff example. The #755 `B_2[2]` port now restores
persistent field-sensitive face behavior in a Sidon-adjacent bounded-sum
deformation: `12 <= n <= 40` split in `28/29` exact finite rows. The `B_2[3]`
capacity-lane deformation split in `18/19` rows over `12 <= n <= 30`. The
#755 near-ground D2 passes show large smooth `h-1/h-2` layers for both
`B_2[2]` and `B_2[3]`, so the exact-face reading remains the safer one. The
near-ground ratio gate adds one finite diagnostic: increasing capacity from
`g=2` to `g=3` preserves exact-face field sensitivity. The relative
near-ground bath compresses inside fixed-`h` plateaus, then resets upward when
`h(n)` jumps and the new ground face is sparse. This reset now appears in both
lanes and repeats: `B_2[3]` at `n=31`, `B_2[2]` at `n=36`, and `B_2[3]`
again at `n=35`. Both lanes then compress again inside the new plateau,
through `B_2[2] n=40` and `B_2[3] n=38`. The plateau-slope extraction gives
the cleanest finite regularity in `B_2[3]`: the observed `h=14,15,16`
plateaus compress `ln((h-1)/ground)` at about `-1.1` to `-1.2` per added site.
The `B_2[3] n=40` ground-only packet then found the next sparse reset:
`h=17` with only `8` exact maximizers. A corrected imported-ground D2 packet
shows the first post-reset shell is already large: `810,894` exact `h-1`
states and at least `1,000,000` `h-2` states. Because the `h-2` layer hit its
cap, this is lower-bound evidence, not an exact D2 ratio row.
The `B_2[2] n=43` layer-capped packet then produced the paired reset in the
lower-capacity lane: `h=14`, `18` exact maximizers, `352,990` exact `h-1`
states, and a capped `h-2` shell. The paired resets are now the main finite
object to explain.
The `B_2[3] n=41` ground-only packet keeps `h=17` and expands the ground face
to `246` maximizers. Its joint-best score is about `0.1105`, so the face is
field-sensitive but not near-zero joint-compatible. That variation is now part
of the honest story: resets and plateaus persist, but the field frontier can be
costly inside the plateau.
The `B_2[2] n=44` layer-capped row gives the same next-step behavior in the
lower-capacity lane: `h=14`, `122` exact maximizers, `997,280` exact `h-1`
states, and a capped `h-2` layer. Its normalized first-excited ratio falls
relative to the reset row, so the reset-to-plateau-compression mechanism is now
visible one step past the reset.
The imported-ground D2 row for `B_2[3] n=41` is censored in both retained
near-ground layers: `h-1 >= 1,000,000` and `h-2 >= 1,000,000`. This supports
rapid plateau expansion after the `n=40` reset, but it does not supply an exact
compression slope for the `h=17` plateau.
The `B_2[2] n=45` row reaches the same censored regime in the lower-capacity
lane: `h=14`, `724` exact maximizers, and both retained near-ground layers at
the million cap. Exact D2 ratio tracking is now capped out unless we raise the
per-layer budget.
The `B_2[3] n=44` ground-only row extends the exact-face plateau sequence:
`n=40,41,42,43,44` have `8 -> 246 -> 3,134 -> 32,002 -> 212,586` exact
maximizers at `h=17`. That is now the cleanest finite example of reset followed
by plateau expansion in the higher-capacity lane. The joint-best score also
compresses across the plateau, from about `0.1105` at `n=41` to `0.0748` at
`n=44`.
The `B_2[3] n=45` ground-only row then jumps to `h=18` and resets the exact
ground face back to `8` maximizers. This closes one full finite cycle in the
higher-capacity lane: sparse reset, plateau expansion, next sparse reset.
The `B_2[3] n=46` row begins the new plateau: `h=18`, `142` exact maximizers.
It is also a non-split row for the tracked top-k fields, so the current
language should emphasize reset/plateau mechanics with field-sensitive
exceptions, not universal field splitting.
The `B_2[2] n=46` layer-capped row extends the lower-capacity plateau:
`n=43,44,45,46,47` have `18 -> 122 -> 724 -> 3,504 -> 16,036` exact
maximizers at `h=14`. Both lanes now show the same reset-to-expansion motif on
their exact maximizer faces.

## Next Engineering Move

Generalize the layer-pruned transfer operator:

```text
port the same exclusion-rule API to #755 B_h[g]
```

The #755 report answered the first version of this question:

> Does bounded repeated-sum capacity preserve Sidon-like field-sensitive face
> behavior better than the rigid sum-free rule?

The first exact scouts say yes for `B_2[2]` and `B_2[3]`, with isolated
non-split rows. The next report should now answer one scaling question only:

> Does the plateau/reset pattern repeat at the next handoff, especially
> do plateau-compression slopes stabilize across the observed `h` plateaus, or
> do we need a capped/tilted sampler before pushing D2 enumeration farther?

The first capped/layer-capped scouts show the next bottleneck: balanced D2
lower bounds are now cheap enough, but exact ground enumeration for `g=3`
dominates beyond `n=39`. The `B_2[3] n=40` ground-only row found the next
sparse reset: `h=17` with only `8` exact maximizers. The imported-ground
follow-up removes that repeated-ground bottleneck for already verified rows,
but the next new `g=3` site still needs a ground pass first. The `B_2[2] n=43`
row now gives the matched lower-capacity reset; the natural next scale moves
are `B_2[3] n=47` ground-only and `B_2[2] n=48` ground-only.

If yes, then the honest theorem-language candidate remains:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> isolated field-sensitive handoffs.

Still unsafe: any claim that PMF supplies a Sidon proof.

## 2026-05-01 Update: Ground-Face Branch Window, n = 58..64

The next exact gate is now packet-backed:

```text
EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01
EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-64-2026-05-01
```

The transfer-operator ground-face export matched exact `h(n)` and exact
maximizer counts for every row `59 <= n <= 64`. The derived branch certificate
then shows that every post-58 row contains the previous exact face and the
previous exact face shifted by `+1`. Difference-skeleton counts grow:

```text
n:         58   59   60   61    62    63    64
skeletons:  2    4   18   49   123   312   669
```

So `58` is no longer best described as an isolated split. It is the first
visible branch point of a growing exact-face fan: inherited persistence,
translated persistence, and rapidly accumulating new skeleton branches.

The claim ceiling remains finite and interpretive. This is exact face geometry,
not a Sidon theorem.
