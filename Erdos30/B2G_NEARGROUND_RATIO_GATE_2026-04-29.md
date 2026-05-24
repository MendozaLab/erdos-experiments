# B_2[g] Near-Ground Ratio Gate

Date: 2026-04-29
Problem: Erdős #755 candidate lane
Status: DERIVED_ANALYSIS_FROM_VERIFIED_PACKETS

## Question

Do normalized near-ground ratios separate the `B_2[2]` and `B_2[3]` capacity
lanes?

## Source Packets

- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-30-32-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-27-30-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-33-34-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-31-32-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-35-36-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-33-34-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-37-38-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-35-36-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-39-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-37-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-NEARGROUND-D2-40-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-NEARGROUND-D2-38-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-GROUNDONLY-40-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-IMPORTED-GROUND-LAYER-CAPPED-D2-40-FIXED-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-43-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-44-2026-04-29`
- `EXP-MM-755-PMF-B2G3-TRANSFER-IMPORTED-GROUND-LAYER-CAPPED-D2-41-2026-04-29`
- `EXP-MM-755-PMF-B2G2-TRANSFER-LAYER-CAPPED-D2-45-2026-04-29`

All listed packets have verified SHA256 sidecars.

## Result

Yes, but not as a monotone law.

The absolute near-ground layers grow in both capacity lanes. The normalized
ratios, however, fall sharply as the ground face becomes less sparse. At the
shared point `n = 30`, `B_2[2]` has a tiny ground face (`6` states), so its
near-ground bath is enormous relative to ground. `B_2[3]` has a much larger
ground face (`52,958` states), so its normalized bath ratio is far smaller even
though the absolute bath is larger.

The direct-overlap extensions add the important twist: normalized bath
compression resets when `h(n)` jumps. In `B_2[3]`, `h(n)` jumps from `14` to
`15` at `n = 31`, and the ground face collapses from `52,958` states at
`n = 30` to `36` states at `n = 31`. In `B_2[2]`, the same reset appears when
`h(n)` jumps from `12` to `13` at `n = 36`, where the ground face collapses
from `24,728` states at `n = 35` to `2` states at `n = 36`. The ratios spike
again. The next extension repeats the pattern again: `B_2[3]` jumps from
`h=15` to `h=16` at `n = 35`, where the ground face drops from `49,534` states
at `n = 34` to `14` states at `n = 35`, and the ratio spikes again. So the
finite rule is not "capacity monotonically compresses the bath"; it is:

```text
inside a fixed-h plateau -> ground face expands -> relative bath compresses
at an h-jump -> new sparse ground face -> relative bath resets high
```

| g | n | h | ground | h-1 | h-2 | (h-1)/ground | ln ratio h-1 | (h-2)/ground | ln ratio h-2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 30 | 12 | 6 | 27,856 | 1,195,604 | 4,642.667 | 8.443 | 199,267.333 | 12.202 |
| 2 | 31 | 12 | 40 | 89,612 | 2,311,480 | 2,240.300 | 7.714 | 57,787.000 | 10.965 |
| 2 | 32 | 12 | 218 | 240,288 | 4,195,380 | 1,102.239 | 7.005 | 19,244.862 | 9.865 |
| 2 | 33 | 12 | 1,372 | 614,630 | 7,427,762 | 447.981 | 6.105 | 5,413.821 | 8.597 |
| 2 | 34 | 12 | 6,014 | 1,390,680 | 12,532,552 | 231.240 | 5.443 | 2,083.896 | 7.642 |
| 2 | 35 | 12 | 24,728 | 3,027,438 | 20,703,726 | 122.430 | 4.808 | 837.258 | 6.730 |
| 2 | 36 | 13 | 2 | 81,762 | 6,091,794 | 40,881.000 | 10.618 | 3,045,897.000 | 14.929 |
| 2 | 37 | 13 | 30 | 255,392 | 11,888,022 | 8,513.067 | 9.049 | 396,267.400 | 12.890 |
| 2 | 38 | 13 | 238 | 684,478 | 21,837,994 | 2,875.958 | 7.964 | 91,756.277 | 11.427 |
| 2 | 39 | 13 | 1,570 | 1,768,784 | 39,247,586 | 1,126.614 | 7.027 | 24,998.462 | 10.127 |
| 2 | 40 | 13 | 7,400 | 4,096,168 | 67,361,252 | 553.536 | 6.316 | 9,102.872 | 9.116 |
| 3 | 27 | 14 | 48 | 139,090 | 2,846,606 | 2,897.708 | 7.972 | 59,304.292 | 10.990 |
| 3 | 28 | 14 | 752 | 495,394 | 6,316,830 | 658.769 | 6.490 | 8,400.040 | 9.036 |
| 3 | 29 | 14 | 8,652 | 1,614,908 | 13,413,866 | 186.651 | 5.229 | 1,550.377 | 7.346 |
| 3 | 30 | 14 | 52,958 | 4,417,694 | 26,461,952 | 83.419 | 4.424 | 499.678 | 6.214 |
| 3 | 31 | 15 | 36 | 278,242 | 11,381,826 | 7,728.944 | 8.953 | 316,161.833 | 12.664 |
| 3 | 32 | 15 | 548 | 1,075,494 | 26,212,384 | 1,962.580 | 7.582 | 47,832.818 | 10.775 |
| 3 | 33 | 15 | 7,164 | 3,807,894 | 57,852,068 | 531.532 | 6.276 | 8,075.386 | 8.997 |
| 3 | 34 | 15 | 49,534 | 11,161,238 | 118,073,652 | 225.325 | 5.418 | 2,383.689 | 7.776 |
| 3 | 35 | 16 | 14 | 300,424 | 30,868,234 | 21,458.857 | 9.974 | 2,204,873.857 | 14.606 |
| 3 | 36 | 16 | 154 | 1,317,756 | 75,689,678 | 8,556.857 | 9.054 | 491,491.416 | 13.105 |
| 3 | 37 | 16 | 2,520 | 5,263,480 | 177,483,736 | 2,088.683 | 7.644 | 70,430.054 | 11.162 |
| 3 | 38 | 16 | 20,824 | 17,159,770 | 382,339,504 | 824.038 | 6.714 | 18,360.522 | 9.818 |

The `B_2[3] n=40` imported-ground follow-up is not part of the exact D2 ratio
table because its `h-2` layer is capped. It does, however, confirm the next
sparse reset:

| g | n | h | ground | exact h-1 | h-2 lower | (h-1)/ground | ln ratio h-1 | h-2 exact? |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 3 | 40 | 17 | 8 | 810,894 | >= 1,000,000 | 101,361.750 | 11.526 | no |
| 2 | 43 | 14 | 18 | 352,990 | >= 1,000,000 | 19,610.556 | 9.884 | no |
| 2 | 44 | 14 | 122 | 997,280 | >= 1,000,000 | 8,174.426 | 9.009 | no |
| 3 | 41 | 17 | 246 | >= 1,000,000 | >= 1,000,000 | >= 4,065.041 | >= 8.310 | no |
| 2 | 45 | 14 | 724 | >= 1,000,000 | >= 1,000,000 | >= 1,381.215 | >= 7.231 | no |

## Interpretation

The ratio gate does not prove an asymptotic law. It does show a useful finite
diagnostic:

```text
fixed-h plateau -> ground face expands -> relative near-ground bath compresses
h-jump -> sparse new ground face -> relative near-ground bath resets
```

This is the mature Atlas framing in miniature. The physics analogy proposes a
low-description-length pattern, reset followed by plateau compression. The
mathematical object is the exact finite transfer state space, and the report
only promotes what the packets actually show. No physics term is allowed to
outrank the counts.

That is better than the first reading. It is more lattice-gas-like: the system
has plateaus and handoffs, not a single smooth curve. Raising capacity does not
erase field-sensitive exact-face selection. But the relative thermal bath is
controlled by where the row sits in the current `h(n)` plateau. The repeated
reset now appears three times: `B_2[3]` at `n=31`, `B_2[2]` at `n=36`, and
`B_2[3]` again at `n=35`. That is the strongest finite evidence so far for a
handoff mechanism in the #755 lane.

The safer language is entropy-gap language, not raw ratio language. The log
ratios are still positive and large, so the near-ground bath remains broad; it
is just less dominant relative to ground in the higher-capacity lane.

Engineering note: the `B_2[3]`, `n=38`, D2 packet required `8,217,427,559`
second-pass nodes and retained `399,520,098` terminal near-ground states. The
scanner now supports `--progress-every N`, which emits ground-pass and
near-ground-pass counters to stderr. Any larger D2 row should be run as a
single-row packet with this option enabled, and the next engineering target is
a cheaper diagnostic that avoids retaining the full D2 bath.

The first capped scout is recorded in
`B2G_CAPPED_NEARGROUND_SCOUT_2026-04-29.md`. It confirms that the next rows hit
a 10,000,000 retained-state D2 cap, but global caps are order-biased and do not
preserve exact compression ratios. The layer-capped follow-up fixes that
engineering problem enough to provide balanced lower bounds: both `h-1` and
`h-2` hit a `1,000,000` per-layer cap for `B_2[2] n=41`, `B_2[2] n=42`, and
`B_2[3] n=39`. The `B_2[3] n=40` ground-only row found the next sparse reset:
`h=17` with only `8` exact maximizers. The corrected imported-ground D2 packet
then showed `810,894` exact first-excited states and at least `1,000,000`
second-excited states. That row is a lower-bound scout, not an exact D2 ratio
row, but it strengthens the handoff story: the ground face resets sparse while
the near-ground bath remains enormous. The next `B_2[2]` row, `n=43`, repeats
that same reset shape: `h=14`, only `18` exact maximizers, `352,990` exact
first-excited states, and a capped second-excited shell. At `B_2[2] n=44`,
the ground face expands to `122` and the normalized first-excited ratio falls
from `19,610.556` to `8,174.426`, matching the plateau-compression reading.
At `B_2[3] n=41`, both retained near-ground layers hit the per-layer cap, so
the ratio row is censored: the lower bound is already `>= 4,065.041`, but the
exact compression rate is not known from this packet.
At `B_2[2] n=45`, the same thing happens in the lower-capacity lane. Exact D2
ratio tracking is therefore clean through `n=44`; beyond that, the current
million-per-layer scout gives lower bounds unless the cap is raised.

The plateau-slope extraction is recorded in
`B2G_PLATEAU_COMPRESSION_SLOPES_2026-04-29.md`. The cleanest finite signal is
`B_2[3]`: across the observed `h=14,15,16` plateaus, the slope of
`ln((h-1)/ground)` stays near `-1.1` to `-1.2` per added site.

## Claim Boundary

Safe internal claim:

> In exact finite `B_2[g]` transfer-state packets, increasing capacity from
> `g=2` to `g=3` preserves exact-face field sensitivity. Within fixed-`h`
> plateaus, the normalized near-ground bath ratios shrink as the ground face
> expands; at `h(n)` jumps, those ratios reset upward.

Unsafe:

> The ratio trend proves a #755 theorem.

Unsafe:

> The observed finite entropy gaps establish asymptotic thermodynamics.

## Next Gate

Run a direct overlap window so the comparison is less confounded by different
`n` ranges:

```text
B_2[2] layered-capped D2 at n = 43
B_2[3] ground-only at n = 41, or B_2[2] layered-capped D2 at n = 44
```

The question is no longer whether the reset/plateau pattern exists; it repeats.
The next question is whether its compression rate has a stable finite scaling
law, or whether the observed slopes are only local to each plateau.
