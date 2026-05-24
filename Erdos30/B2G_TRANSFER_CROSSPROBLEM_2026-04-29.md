# B_2[2] Transfer-State Cross-Problem Gate

Date: 2026-04-29
Problem: Erdős #755 candidate lane
Reference problem: Erdős #30 PMF transfer operator
Status: EXACT_SCOUT_PASS, FIELD_SENSITIVITY_PERSISTENT_WITH_EXCEPTIONS

## Question

Does bounded repeated-sum capacity preserve Sidon-like field-sensitive face
behavior better than the rigid sum-free rule?

## Answer

Yes, in the first exact finite scout.

The #166 sum-free port told us the transfer API was not just hallucinating the
same story everywhere: its ground face stayed rigid. The #755 `B_2[2]` port
does the opposite. It keeps the same local-state style as Sidon, but relaxes
the hard Sidon constraint from "no repeated sums" to "bounded repeated sums."
Across `12 <= n <= 40`, the split persisted in `28/29` exact finite rows. The
single non-split row was `n = 37`. That matters: the signal is persistent, not
universal row-by-row.

That is the useful shape: #166 is a negative control, #755 is the closer
Sidon-adjacent deformation, and #30 remains the target.

The epistemic rule is the same as the broader Atlas rule:

> The analogy is not the proof; it is the microscope.

Here the physics/MDL language suggested what to measure: ground-state
degeneracy, field response, reset rows, plateau expansion, and near-ground
shells. The evidence is not the language. The evidence is the Rust operator
matching exact finite configuration counts, then separating #166, #755, and #30
by their actual local exclusion laws.

## Operator Port

The new binary lives in the same Rust crate:

```text
erdos-experiments/Erdos30/rust-transfer-operator/src/bin/b2g_transfer.rs
```

State:

```text
State = {
  occupied prefix in [0,n],
  sum_counts[s] = ordered representation count for a + b = s
}
```

Transition rule:

```text
occupy x iff every updated ordered representation count remains <= 2g
```

For `g = 2`, the ordered cap is `4`. The self-pair `(x,x)` adds `1`; each
existing occupied site `a` adds the two ordered pairs `(a,x)` and `(x,a)`.

## Packets

- `EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-12-20-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-21-30-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-31-35-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-SCOUT-36-40-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-SCOUT-12-24-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-SCOUT-25-30-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-30-32-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-27-30-2026-04-29`

All SHA256 sidecars verified `OK`.

## Result Summary

Combined exact finite `B_2[2]` ground-face window:

```text
12 <= n <= 40
checked rows: 29
h(n) range: 7..13
frontier split count: 28 / 29
non-split row: n = 37
total exact maximizers across rows: 81,856
combined runtime: 109.899904666 sec
```

Selected rows:

| n | h(n) | exact maximizers | nodes visited | branch-bound pruned | field split? |
|---|---:|---:|---:|---:|---|
| 12 | 7 | 328 | 4,439 | 1,364 | yes |
| 20 | 9 | 3,998 | 270,434 | 80,546 | yes |
| 25 | 10 | 16,350 | 2,416,053 | 694,159 | yes |
| 30 | 12 | 6 | 17,170,555 | 4,769,079 | yes |
| 35 | 12 | 24,728 | 111,385,011 | 29,581,466 | yes |
| 37 | 13 | 30 | 214,410,083 | 56,210,818 | no |
| 40 | 13 | 7,400 | 626,135,408 | 159,867,698 | yes |

At `n = 35`, the three zero-temperature fields select different witnesses:

```text
prefix: [1, 3, 7, 12, 13, 15, 17, 24, 25, 28, 32, 35]
mass:   [0, 1, 2, 3, 12, 20, 22, 25, 28, 29, 33, 35]
joint:  [1, 3, 7, 12, 13, 15, 17, 24, 25, 28, 32, 35]
```

The `B_2[3]` deformation keeps the same qualitative behavior in the next
capacity lane:

```text
12 <= n <= 30
checked rows: 19
h(n) range: 9..14
frontier split count: 18 / 19
non-split row: n = 13
total exact maximizers across rows: 108,064
combined runtime: 13.685734708 sec
```

The extension block `25 <= n <= 30` split in every row:

| n | h(n) | exact maximizers | nodes visited | branch-bound pruned | field split? |
|---|---:|---:|---:|---:|---|
| 25 | 13 | 5,048 | 6,819,600 | 2,303,863 | yes |
| 26 | 13 | 28,544 | 12,598,080 | 4,193,047 | yes |
| 27 | 14 | 48 | 20,512,275 | 6,785,688 | yes |
| 28 | 14 | 752 | 31,647,045 | 10,408,049 | yes |
| 29 | 14 | 8,652 | 52,680,018 | 17,164,381 | yes |
| 30 | 14 | 52,958 | 91,765,283 | 29,546,540 | yes |
```

The first `B_2[2]` near-ground pass at `30 <= n <= 32` shows the lower layers
are large and smooth:

| n | h(n) | ground | h-1 | h-2 | retained near-ground states |
|---|---:|---:|---:|---:|---:|
| 30 | 12 | 6 | 27,856 | 1,195,604 | 1,223,466 |
| 31 | 12 | 40 | 89,612 | 2,311,480 | 2,401,132 |
| 32 | 12 | 218 | 240,288 | 4,195,380 | 4,435,886 |

The `B_2[3]` near-ground D2 pass shows the same pattern:

| n | h(n) | ground | h-1 | h-2 | retained near-ground states |
|---|---:|---:|---:|---:|---:|
| 27 | 14 | 48 | 139,090 | 2,846,606 | 2,985,744 |
| 28 | 14 | 752 | 495,394 | 6,316,830 | 6,812,976 |
| 29 | 14 | 8,652 | 1,614,908 | 13,413,866 | 15,037,426 |
| 30 | 14 | 52,958 | 4,417,694 | 26,461,952 | 30,932,604 |

The normalized ratio gate is recorded in
`B2G_NEARGROUND_RATIO_GATE_2026-04-29.md`. Its finite reading is:

```text
fixed-h plateau -> ground face expands -> relative near-ground bath compresses
h-jump -> sparse new ground face -> relative near-ground bath resets
```

At the shared point `n = 30`, `(h-1)/ground` drops from `4,642.667` for
`B_2[2]` to `83.419` for `B_2[3]`. The absolute `h-1` bath is larger for
`B_2[3]`, so this is not a collapse of near-ground states. It is a compression
relative to a much larger ground face. But the direct-overlap extension shows
that the compression is plateau-local: when `B_2[3]` jumps from `h=14` to
`h=15` at `n=31`, the ground face resets small and `(h-1)/ground` spikes to
`7,728.944`. The same reset repeats for `B_2[2]` at `n=36`, where `h` jumps
from `12` to `13`, the ground face drops from `24,728` states to `2`, and
`(h-1)/ground` spikes to `40,881.000`.
The next extension repeats the same mechanism for `B_2[3]`: after the `h=15`
plateau compresses through `n=34`, `h` jumps to `16` at `n=35`, the ground
face resets to `14`, and `(h-1)/ground` spikes to `21,458.857`.
Both lanes then compress again inside the new plateau: `B_2[2]` drops to
`553.536` by `n=40`, and `B_2[3]` drops to `824.038` by `n=38`.
The derived slope report `B2G_PLATEAU_COMPRESSION_SLOPES_2026-04-29.md`
shows the cleanest finite regularity in `B_2[3]`: across the observed
`h=14,15,16` plateaus, `ln((h-1)/ground)` compresses at roughly `-1.1` to
`-1.2` per added site.

## Interpretation

This is the first clean "closer mountain" after #30.

The sum-free operator was useful because it *did not* show Sidon-like
field-sensitivity. That protected us from the sloppy conclusion that every
additive exclusion system has the same PMF face. The `B_2[2]` operator is more
informative because it is a relaxed Sidon law: the same sum-collision memory is
still present, but capacity is finite instead of forbidden.

So the discovery trail now has three rungs:

```text
#166 sum-free: API ports, ground face rigid
#755 B_2[2]/B_2[3]: API ports, field-sensitive face persists across capacity lanes
#30 Sidon: exact target, field-sensitive handoff candidate
```

That is exactly the cross-problem behavior we wanted from a collider: not
universal sameness, but related systems separating by local exclusion law.
The near-ground runs now add an important caution: the signal is clearest on
the exact maximizer face. One layer down, the state space becomes a broad
thermal bath.

The latest `B_2[3] n=40` turn adds the useful engineering twist. A ground-only
packet first found the next reset, `h=17` with only `8` exact maximizers. Then
the scanner imported that verified ground result and ran the D2 shell without
repeating the ground pass. The corrected packet gives:

```text
h=17
ground states=8
h-1 states=810,894 exact
h-2 states>=1,000,000 capped lower bound
```

That is the strongest post-reset evidence in the `g=3` lane so far. It also
forced a scanner fix: per-layer caps must mark the packet censored whenever any
capped layer skips terminal states, even if another layer remains below cap.

The `g=2` lane then produced the corresponding next reset at `n=43`:

```text
h=14
ground states=18
h-1 states=352,990 exact
h-2 states>=1,000,000 capped lower bound
```

So the current shape is not a one-lane accident. Both bounded-sum lanes now
show the same reset motif: a tiny new ground face at the `h(n)` jump, with a
large first-excited shell immediately below it.

The next `g=3` ground-only row, `n=41`, keeps `h=17` and expands the exact
ground face from `8` to `246` maximizers. It remains field-sensitive, but the
joint-best score is about `0.1105`, so this is plateau expansion with a costly
joint frontier, not a near-perfect joint witness.

The next `g=2` row, `n=44`, does the same thing in the lower-capacity lane:
`h=14`, ground face `122`, exact `h-1=997,280`, and capped `h-2`. The
normalized first-excited ratio drops from `19,610.556` at the reset row to
`8,174.426`, so the plateau-compression reading survives the next step.

The imported-ground D2 row for `g=3 n=41` is more censored: both `h-1` and
`h-2` hit the `1,000,000` per-layer cap. That still supports plateau expansion
after the `n=40` reset, but it does not give an exact ratio. The safe statement
is `h-1 / ground >= 4,065.041`, not an equality.

The `g=2 n=45` row reaches the same censored regime: `h=14`, ground face
`724`, and both retained near-ground layers hit the million cap. The lower
bound `h-1 / ground >= 1,381.215` keeps the compression direction compatible
with the plateau story, but the exact ratio is now beyond the current scout
budget.

The `g=3 n=44` ground-only row continues the exact-face expansion from the
reset: `8 -> 246 -> 3,134 -> 32,002 -> 212,586` maximizers across
`n=40,41,42,43,44`, all at `h=17`. That is the cleanest finite
reset-to-plateau sequence in the higher-capacity lane. Its joint-best score
also falls from about `0.1105` at `n=41` to `0.0748` at `n=44`, so the
field-tension observable is compressing as the face thickens.

The next row closes the loop: `g=3 n=45` jumps to `h=18` and resets the exact
ground face to `8` maximizers. So the high-capacity lane now has a complete
finite cycle:

```text
h=17 reset at n=40 -> plateau expansion through n=44 -> h=18 reset at n=45
```

This is the strongest evidence so far that the reset/plateau mechanism is a
real finite operator phenomenon rather than a one-off visual pattern.
The `g=3 n=46` row begins the next plateau: it stays at `h=18` and expands to
`142` exact maximizers. Unlike most prior rows, its tracked top-k fields do not
split, so the safe finite claim is reset/plateau structure with field-sensitive
exceptions, not universal field splitting.

The `g=2` lane now gives the matching lower-capacity sequence:
`18 -> 122 -> 724 -> 3,504 -> 16,036` maximizers across `n=43,44,45,46,47`,
all at `h=14`. The near-ground layers are already capped by `n=45`, so the
finite object we can currently track exactly is the ground-face expansion, not
full D2 ratios.

## Claim Boundary

Safe internal claim:

> The PMF transfer-state API ports from Sidon to the `B_2[2]` bounded-sum
> deformation, and exact finite scouts show persistent field-sensitive
> ground-face selection across `12 <= n <= 40`, with one non-split row.

Safe extension:

> The same qualitative field-sensitive selection persists in the `B_2[3]`
> capacity lane across `12 <= n <= 30`, with one non-split row.

Still unsafe:

> PMF proves Sidon.

Still unsafe:

> #755 is solved.

Still unsafe:

> The `B_2[2]` finite split pattern proves an asymptotic theorem.

## Next Gate

The next honest move is scale and near-ground structure, not theorem language.

Run one of:

```text
B_2[3] ground-only at n = 47
B_2[2] ground-only at n = 48
```

The theorem-language candidate remains narrow:

> finite exact maximizer faces exhibit lattice-gas ground-state degeneracy with
> field-sensitive handoffs under Sidon-adjacent bounded-sum laws.
