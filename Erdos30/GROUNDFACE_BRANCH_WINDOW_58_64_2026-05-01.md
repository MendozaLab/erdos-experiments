# Ground-Face Branch Window, n = 58..71

Date: 2026-05-01
Problem: Erdos #30
Status: EXACT_PACKET_BACKED / FINITE / INTERPRETIVE

## Source Packets

- `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-64-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-65-68-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-68-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-68-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-69-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-69-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-69-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-70-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-70-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-70-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-BRANCH-CERTIFICATE-58-71-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-71-2026-05-01`

All cited branch-window claims below are derived from packet-backed exported
ground faces. The `59 <= n <= 64` transfer packet matched exact `h(n)` and
exact maximizer counts for all six rows. The `65 <= n <= 68` transfer packet
matched exact `h(n)` and exact maximizer counts for all four rows.
The `n = 69` transfer packet matched exact `h(69) = 10` and exact maximizer
count `66,412`. The `n = 70` transfer packet matched exact `h(70) = 10` and
exact maximizer count `117,202`. The `n = 71` transfer packet matched exact
`h(71) = 10` and exact maximizer count `203,840`.

## Result

The `n = 58` branch is not a one-off artifact. From `n = 59` through `n = 71`,
each exact ground face contains:

- every witness from the previous row unchanged;
- every witness from the previous row shifted by `+1`;
- additional new witnesses that rapidly increase the number of difference
  skeletons.

This is the clearest current finite formulation:

> after the Singer-modulus handoff, the exact ground face grows by inherited
> persistence plus translated persistence, while new skeleton branches
> accumulate around that inherited core.

## Branch Table

| n | exact maximizers | difference skeletons | previous face | previous face +1 | Pareto count |
|---:|---:|---:|---:|---:|---:|
| 58 | 10 | 2 | 0 | 0 | 2 |
| 59 | 18 | 4 | 10 | 10 | 2 |
| 60 | 54 | 18 | 18 | 18 | 1 |
| 61 | 152 | 49 | 54 | 54 | 4 |
| 62 | 398 | 123 | 152 | 152 | 6 |
| 63 | 1022 | 312 | 398 | 398 | 10 |
| 64 | 2360 | 669 | 1022 | 1022 | 11 |
| 65 | 5018 | 1329 | 2360 | 2360 | 47 |
| 66 | 9994 | 2488 | 5018 | 5018 | 40 |
| 67 | 19418 | 4712 | 9994 | 9994 | 120 |
| 68 | 36234 | 8408 | 19418 | 19418 | 109 |
| 69 | 66412 | 15089 | 36234 | 36234 | 267 |
| 70 | 117202 | 25395 | 66412 | 66412 | 174 |
| 71 | 203840 | 43319 | 117202 | 117202 | 452 |

The skeleton count strictly increases across the window. Field-selection split
is present in every row. At `n = 70`, the Pareto face has `174` witnesses
across `173` difference skeletons, so the field-selected surface remains tiny
relative to the full exact face while no longer being skeleton-injective. At
`n = 71`, the Pareto face rebounds to `452` witnesses across `452` skeletons.

For `65 <= n <= 71`, full witness faces were exported but the all-pairs
distance graph was intentionally skipped by cap. At `n = 71`, the distance graph
would require `20,775,270,880` edges. This is a performance boundary on geometry
serialization, not a truncation of the exact face.

## Interpretation

This strengthens the finite face-handoff story. At `57`, the exposed optimizer
handoff happens inside one Singer-modulus skeleton. At `58`, the first local
second skeleton appears. From `59` onward, the face behaves like a branching
ground-state fan: old states persist, translated old states persist, and new
skeletons grow around them.

The physics analogy remains zero-temperature field selection on a degenerate
ground face. The data now supports "branching exact-face fan" better than
"single defect row."

## Reset Boundary at n = 72

The next exact reference row is not a continuation of the same `h = 10` face.
Packet `EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01` shows `h(72) = 11` with
only `4` exact maximizers. Packet `EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01`
reproduces that row in the transfer lane with full face export and parity match.

The separate reset note is:

```text
GROUNDFACE_RESET_72_2026-05-01.md
```

## Phi / Self-Replication Litmus

Follow-up packet:

```text
EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-64-2026-05-01
EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-68-2026-05-01
EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-69-2026-05-01
EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-70-2026-05-01
EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-58-71-2026-05-01
```

These packets separate two claims:

- self-replication: supported in the finite exact-face sense;
- phi/Fibonacci mediation: not supported on the current finite window.

The replication rule passes through `n = 71` because every post-58 face contains
the inherited union of the previous face and the previous face shifted by `+1`.
The phi label still does not pass: face-count ratios, skeleton-count ratios, and
Fibonacci recurrence residuals are not within the stated tolerance.

The honest language is therefore:

> inherited-plus-translated self-replication with new branch production.

Not:

> phi-mediated Fibonacci growth.

## Claim Boundary

Safe:

- `finite exact ground faces exhibit inherited and translated persistence`
- `difference-skeleton branches accumulate after n = 58`
- `field-selected Pareto faces remain small relative to the full exact face`
- `this is a proof-candidate structure`
- `self-replication in the finite exact-face sense`

Unsafe:

- `PMF proves Sidon`
- `physics solves Erdos #30`
- `phase transition proved`
- `asymptotic theorem`
- `phi/Fibonacci mediation`
