# Ground-Face Reset and Plateau Restart at n = 72..91

Date: 2026-05-01
Problem: Erdos #30
Status: EXACT_PACKET_BACKED / FINITE / RESET_EVENT / INTERPRETIVE

## Source Packets

- `EXP-MM-030-RUST-TOPK-FRONTIER-72-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-72-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-73-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-73-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-73-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-74-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-74-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-74-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-75-LB11-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-75-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-75-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-76-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-76-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-76-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-77-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-77-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-77-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-78-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-78-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-78-V2-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-78-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-79-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-79-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-79-V2-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-79-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-80-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-80-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-80-V2-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-80-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-81-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-81-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-81-V2-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-81-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-82-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-82-PAR8-2026-05-01`
- `EXP-MM-030-PMF-GROUNDFACE-RESET-CERTIFICATE-72-82-V2-2026-05-01`
- `EXP-MM-030-PMF-PHI-REPLICATION-LITMUS-72-82-2026-05-01`
- `EXP-MM-030-RUST-TOPK-FRONTIER-83-LB11-PAR8-2026-05-01`
- `EXP-MM-030-PMF-PLATEAU-EDGE-LIGHT-83-FROM82-2026-05-01`
- `EXP-MM-030-PMF-PLATEAU-EDGE-MASKSRC-83-FROM82-2026-05-01`
- `EXP-MM-030-PMF-PLATEAU-EDGE-MASKSRC-84-FROM83-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-85-FROM84-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-86-FROM85-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-87-FROM86-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-88-FROM87-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-89-FROM88-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-90-FROM89-2026-05-01`
- `EXP-MM-030-RUST-EXACT-PLATEAU-EDGE-MASKSRC-91-FROM90-2026-05-01`

All cited result files have verified SHA sidecars. The exact-maximizer packets
supply the independent reference rows used by the transfer operator. Where a
transfer packet is cited, it reproduces those rows with parity checks:

- `checked_count = 1` for each transfer packet
- `h_n_match_count = 1` for each transfer packet
- `maximizer_count_match_count = 1` for each transfer packet
- `mismatch_ns = []`

## Result

The `58 <= n <= 71` branching exact-face fan does not simply continue at `n =
72`. Instead, the extremal cardinality jumps:

| n | h(n) | exact maximizers | difference skeletons | Pareto count | distance edges |
|---:|---:|---:|---:|---:|---:|
| 72 | 11 | 4 | 2 | 1 | 6 |
| 73 | 11 | 8 | 2 | 1 | 28 |
| 74 | 11 | 34 | 13 | 1 | 561 |
| 75 | 11 | 84 | 25 | 1 | 3486 |
| 76 | 11 | 214 | 65 | 1 | 22791 |
| 77 | 11 | 482 | 134 | 2 | 115921 |
| 78 | 11 | 970 | 244 | 3 | 469965 |
| 79 | 11 | 1974 | 502 | 8 | 1947351 |
| 80 | 11 | 4030 | 1028 | 7 | 8118435 |
| 81 | 11 | 8214 | 2092 | 35 | 33730791 |
| 82 | 11 | 15958 | 3872 | 46 | 127320903 |
| 83 | 11 | 30510 | NA_MASK_SOURCE | NA_MASK_SOURCE | NOT_EXPORTED |
| 84 | 11 | 56110 | NA_MASK_SOURCE | NA_MASK_SOURCE | NOT_EXPORTED |
| 85 | 12 | 2 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 86 | 12 | 4 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 87 | 12 | 6 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 88 | 12 | 8 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 89 | 12 | 10 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 90 | 12 | 14 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |
| 91 | 12 | 28 | NA_EXACT_MASK_SOURCE | NA_EXACT_MASK_SOURCE | NOT_EXPORTED |

The full ground face is exported for both rows. Unlike the large `h = 10` faces
at `69`, `70`, and `71`, no distance-edge cap is needed: the complete graph has
six edges at `72`, twenty-eight edges at `73`, `561` edges at `74`, and `3,486`
edges at `75`. The `n = 76` graph is still fully exported with `22,791` edges.
The `n = 77` graph is also fully exported with `115,921` edges. At `n = 78`,
the full witness face is exported, but pairwise distance edges are intentionally
skipped by cap; the implied complete graph would contain `469,965` pairwise
edges, which is not needed for the replication check.
At `n = 79`, the full witness face is again exported and pairwise distance edges
are skipped by cap; the implied complete graph would contain `1,947,351`
pairwise edges.
At `n = 80`, the full witness face is again exported and pairwise distance edges
are skipped by cap; the implied complete graph would contain `8,118,435`
pairwise edges.
At `n = 81`, the full witness face is again exported and pairwise distance edges
are skipped by cap; the implied complete graph would contain `33,730,791`
pairwise edges.
At `n = 82`, the full witness face is again exported and pairwise distance edges
are skipped by cap; the implied complete graph would contain `127,320,903`
pairwise edges.
At `n = 83`, the exact maximizer count is packet-backed, but the transfer run
uses the plateau-edge mask-source mode: it verifies parity and inheritance,
exports compact ground-face masks for downstream inheritance checks, and does
not export the full witness face, skeleton classes, Pareto face, or pairwise
distance graph.

At `n = 84`, the plateau-edge mask-source mode again verifies parity and
inheritance and exports all compact ground-face masks. It does not claim full
skeleton classes, Pareto face, or pairwise distance graph.

At `n = 85`, the upgraded exact-maximizer runner performs the exact scan,
inheritance probe, and compact mask export in one pass. This row is not a
continuation of the `h = 11` plateau: the extremal cardinality jumps to `h =
12`, the exact ground face collapses to two maximizers, and the inherited
`n = 84` face is no longer present in the new ground face.

At `n = 86` through `n = 89`, the same upgraded exact runner shows a new `h = 12`
plateau beginning. All four rows are pure inherited-union rows: they contain the
previous exact face and its `+1` shift, with no additional ground-face
witnesses.

At `n = 90`, branch production re-enters the `h = 12` plateau. The inherited
union from `n = 89` contributes `12` witnesses, and two additional witnesses
appear.

At `n = 91`, branch production accelerates: the inherited union from `n = 90`
contributes `18` witnesses, and ten additional witnesses appear.

The derived reset-window certificate verifies that `n = 73` contains every
`n = 72` witness unchanged and every `n = 72` witness shifted by `+1`, with no
additional witnesses. In other words, the inherited-plus-translated rule restarts
immediately after the reset, but initially without new branch production.

The `n = 74` row keeps the same memory rule but restarts branch production: it
contains all `8` witnesses from `n = 73`, all `8` shifted witnesses from `n =
73`, and `22` additional witnesses. The difference-skeleton count jumps from
`2` to `13`.

The `n = 75` row preserves that rule again: it contains all `34` witnesses from
`n = 74`, all `34` shifted witnesses from `n = 74`, and `24` additional
witnesses. The difference-skeleton count rises from `13` to `25`.

The `n = 76` row continues the same post-reset branch fan: it contains all `84`
witnesses from `n = 75`, all `84` shifted witnesses from `n = 75`, and `80`
additional witnesses. The difference-skeleton count rises from `25` to `65`.

The `n = 77` row preserves the same rule: it contains all `214` witnesses from
`n = 76`, all `214` shifted witnesses from `n = 76`, and `138` additional
witnesses. The difference-skeleton count rises from `65` to `134`, and the
Pareto face has two witnesses.

The `n = 78` row preserves the same rule again: it contains all `482` witnesses
from `n = 77`, all `482` shifted witnesses from `n = 77`, and `220` additional
witnesses. The difference-skeleton count rises from `134` to `244`, and the
Pareto face has three witnesses.

The `n = 79` row preserves the same rule again: it contains all `970` witnesses
from `n = 78`, all `970` shifted witnesses from `n = 78`, and `516` additional
witnesses. The overlap between the old and shifted faces is `482`, so the
inherited union has `1,458` witnesses before the `516` new ones. The
difference-skeleton count rises from `244` to `502`, and the Pareto face widens
to eight witnesses.

The `n = 80` row preserves the same rule again: it contains all `1,974`
witnesses from `n = 79`, all `1,974` shifted witnesses from `n = 79`, and
`1,052` additional witnesses. The overlap between the old and shifted faces is
`970`, so the inherited union has `2,978` witnesses before the `1,052` new
ones. The difference-skeleton count rises from `502` to `1,028`, and the Pareto
face has seven witnesses.

The `n = 81` row preserves the same rule again: it contains all `4,030`
witnesses from `n = 80`, all `4,030` shifted witnesses from `n = 80`, and
`2,128` additional witnesses. The overlap between the old and shifted faces is
`1,974`, so the inherited union has `6,086` witnesses before the `2,128` new
ones. The difference-skeleton count rises from `1,028` to `2,092`, and the
Pareto face widens sharply to `35` witnesses.

The `n = 82` row preserves the same rule again: it contains all `8,214`
witnesses from `n = 81`, all `8,214` shifted witnesses from `n = 81`, and
`3,560` additional witnesses. The overlap between the old and shifted faces is
`4,030`, so the inherited union has `12,398` witnesses before the `3,560` new
ones. The difference-skeleton count rises from `2,092` to `3,872`, and the
Pareto face widens to `46` witnesses.

The `n = 83` row is mask-source certified rather than full-face exported. It
preserves the same rule again: it contains all `15,958` witnesses from `n = 82`,
all `15,958` shifted witnesses from `n = 82`, and `6,808` additional witnesses.
The overlap between the old and shifted faces is `8,214`, so the inherited union
has `23,702` witnesses before the `6,808` new ones. Full skeleton count and full
Pareto count are intentionally not claimed for `n = 83`.

The `n = 84` row is also mask-source certified. It preserves the same rule: it
contains all `30,510` witnesses from `n = 83`, all `30,510` shifted witnesses
from `n = 83`, and `11,048` additional witnesses. The inherited union has
`45,062` witnesses before the `11,048` new ones. Full skeleton count and full
Pareto count are intentionally not claimed for `n = 84`.

The `n = 85` row is the next reset edge. The exact runner finds `h(85) = 12`
with `2` exact maximizers. Against the `n = 84` mask source, the probe reports
`56,110` previous witnesses, `56,110` shifted previous witnesses, an inherited
union of `81,710`, and `0` inherited-union witnesses present in the new `h =
12` ground face. The new face count is therefore `2`, not inherited growth.

The `n = 86` row restarts the plateau exactly: it contains both `n = 85`
witnesses and both `+1` shifted witnesses, with inherited union count `4` and
new face count `0`.

The `n = 87` row remains a pure inherited-union row: it contains all `4`
witnesses from `n = 86`, all `4` shifted witnesses, inherited union count `6`,
and new face count `0`. This differs from the `h = 11` restart, where new branch
production had already resumed by the second step after reset.

The `n = 88` row remains pure inherited-union again: it contains all `6`
witnesses from `n = 87`, all `6` shifted witnesses, inherited union count `8`,
and new face count `0`.

The `n = 89` row remains pure inherited-union again: it contains all `8`
witnesses from `n = 88`, all `8` shifted witnesses, inherited union count `10`,
and new face count `0`.

The `n = 90` row is the first branch-production row in the `h = 12` plateau: it
contains all `10` witnesses from `n = 89`, all `10` shifted witnesses, inherited
union count `12`, and new face count `2`.

The `n = 91` row is the first acceleration row after branch onset: it contains
all `14` witnesses from `n = 90`, all `14` shifted witnesses, inherited union
count `18`, and new face count `10`.

The V2 reset-window certificate uses the correct reset semantics: the first
post-reset step `72 -> 73` is allowed to be a skeleton plateau (`2 -> 2`), while
later rows must grow strictly. Under that rule, all required checks pass:
complete exports, field-selection split, exact previous-face persistence,
`+1` translated persistence, nondecreasing skeleton count, and strict skeleton
growth after the initial plateau.

The phi sidecar remains negative. The `72..82` litmus supports finite
inherited-plus-translated self-replication, but does not support phi/Fibonacci
mediation at the current tolerance.

The `n = 83` and `n = 84` mask-source packets extend the
inherited-plus-translated memory rule, but they are not included in the phi
sidecar because the full skeleton/Pareto faces were not exported. The `n = 85`
packet ends the `h = 11` plateau and is a reset-edge datum, not positive phi
evidence. The `n = 86` through `n = 89` packets show the new `h = 12` plateau
restarts by inheritance only, `n = 90` shows branch production re-entering with
two new witnesses, and `n = 91` shows branch production accelerating to ten new
witnesses. This still does not establish phi/Fibonacci mediation.

## Engineering Update

The transfer/export scanner now has a parallel prefix lane:
`--parallel-prefix-depth`. A control packet on `n = 76` verified parity against
the existing sequential packet: same `h(n)`, same maximizer count, same exported
214-witness face, same `22,791` distance edges, and same top joint frontier.
That makes the parallel lane acceptable for later exact-face exports, provided
each new packet still passes SHA and parity checks.

The exact-maximizer scanner now also has a plateau-edge lane:
`--inheritance-source-results` plus `--export-ground-face-masks`. A validation
packet on `n = 73` reproduced the known `h = 11`, `8`-maximizer row and verified
that the `n = 72` face plus its `+1` shift accounts for all of `n = 73`. The
`n = 85` packet then used the same one-pass lane to avoid a separate
transfer/export traversal.

## Ground-Face Witnesses

| index | witness | prefix residual | density-adjusted mass dev | role |
|---:|---|---:|---:|---|
| 0 | `[0,1,4,13,28,33,47,54,64,70,72]` | `0.1177490061` | `46` | non-Pareto |
| 1 | `[0,1,9,19,24,31,52,56,58,69,72]` | `0` | `41` | prefix winner |
| 2 | `[0,2,8,18,25,39,44,59,68,71,72]` | `0` | `26` | prefix / mass / joint / Pareto winner |
| 3 | `[0,3,14,16,20,41,48,53,63,71,72]` | `1.0883117546` | `31` | non-Pareto |

The four witnesses form two positive-difference skeleton classes: `{0,2}` and
`{1,3}`. Field selection still splits the face: prefix has a two-way zero
residual tie, while mass and joint select witness `2`.

## Interpretation

This is a reset followed by a plateau restart, not a continuation of the `h =
10` branching fan. The prior window showed memory-like face growth: the previous
exact face and its `+1` translate persisted while new branches accumulated. At
`n = 72`, the system pays the cardinality price to add an eleventh particle and
the degeneracy collapses from `203,840` maximizers at `n = 71` to `4` maximizers
at `n = 72`. At `n = 73`, the same two-copy memory rule restarts exactly: old
face plus shifted old face, but no extra new witnesses yet. At `n = 74`, the
two-copy memory rule survives and new branches re-enter. At `n = 75`, the same
post-reset branching behavior continues. At `n = 76`, the same rule still holds
while new branch production accelerates. At `n = 77`, the rule survives again
and the field-selected Pareto face becomes two-point rather than single-point.
At `n = 78`, the same memory rule survives and the Pareto face becomes
three-point, while the full face remains small enough to export completely. At
`n = 79`, the memory rule still survives, but the Pareto face expands more
noticeably to eight witnesses. At `n = 80`, the memory rule still survives and
the exact face roughly doubles again, while the Pareto face remains small at
seven witnesses. At `n = 81`, the memory rule still survives, and the Pareto
face widens sharply to 35 witnesses without a cardinality reset. At `n = 82`,
the memory rule still survives, but face and skeleton growth ratios slow below
the near-doubling pattern from `79..81`. At `n = 83`, the memory rule still
survives and face growth slows again; at `n = 84`, it survives one more row. At
`n = 85`, the cardinality jumps to `12` and the ground face collapses to two new
maximizers, so the inherited `h = 11` face is no longer the ground-state object.
At `n = 86` through `n = 89`, the new plateau restarts by inheritance only,
with exact face counts growing linearly `2, 4, 6, 8, 10`. At `n = 90`, that
linear front breaks: the inherited union would give `12`, but the exact face has
`14`, so two new branches enter. At `n = 91`, the inherited union would give
`18`, but the exact face has `28`, so ten new branches enter.

The physics analogy is a jammed-packing reset or internal-variable relaxation
event: history builds an enormous degenerate face, then a new admissible particle
count opens and the ground state snaps to a small low-degeneracy face. This is
finite evidence for a reset/branching mechanism, not an asymptotic theorem.

## Claim Boundary

Safe:

- `n = 72 is an exact finite reset event`
- `n = 73 starts a reset plateau by exact two-copy replication`
- `n = 74 restarts branch production while preserving exact-plus-shift memory`
- `n = 75 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 76 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 77 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 78 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 79 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 80 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 81 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 82 continues post-reset branch production while preserving exact-plus-shift memory`
- `n = 83 continues post-reset exact-plus-shift memory under the mask-source plateau-edge probe`
- `n = 84 continues post-reset exact-plus-shift memory under the mask-source plateau-edge probe`
- `n = 85 is an exact reset edge from h = 11 to h = 12`
- `n = 86 restarts the h = 12 plateau by exact inherited-plus-shift memory`
- `n = 87 continues the h = 12 plateau by exact inherited-plus-shift memory`
- `n = 88 continues the h = 12 plateau by exact inherited-plus-shift memory`
- `n = 89 continues the h = 12 plateau by exact inherited-plus-shift memory`
- `n = 90 is the first branch-production row in the h = 12 plateau`
- `n = 91 is the first acceleration row after h = 12 branch onset`
- `the extremal cardinality jumps from 10 to 11`
- `the extremal cardinality later jumps from 11 to 12 at n = 85`
- `the exact ground face collapses from 203,840 to 4 maximizers`
- `the h = 11 exact ground face then grows from 4 to 8 to 34 to 84 to 214 to 482 to 970 to 1,974 to 4,030 to 8,214 to 15,958 to 30,510 to 56,110 maximizers`
- `the n = 85 h = 12 reset face has 2 exact maximizers`
- `the h = 12 exact ground face grows linearly from 2 to 4 to 6 to 8 to 10 maximizers through n = 89, then jumps to 14 at n = 90 and 28 at n = 91`
- `the reset face has two difference-skeleton classes`
- `field selection remains visible on the reset face`
- `phi/Fibonacci mediation is not supported by the 72..82 sidecar`
- `n = 83 and n = 84 full skeleton/Pareto counts are not claimed`
- `n = 85 is a finite reset-edge datum, not an asymptotic theorem`
- `n = 86 through n = 89 are finite plateau-restart data, not an asymptotic theorem`
- `n = 90 is finite branch-onset data, not an asymptotic theorem`
- `n = 91 is finite branch-acceleration data, not an asymptotic theorem`

Unsafe:

- `PMF proves Sidon`
- `physics solves Erdos #30`
- `phase transition proved`
- `asymptotic theorem`
- `phi/Fibonacci mediation`
