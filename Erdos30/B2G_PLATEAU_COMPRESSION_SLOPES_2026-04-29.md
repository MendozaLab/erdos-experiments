# B_2[g] Plateau Compression Slopes

Date: 2026-04-29
Problem: Erdős #755 candidate lane
Status: DERIVED_ANALYSIS_FROM_VERIFIED_PACKETS

## Question

Once the plateau/reset mechanism is visible, do the compression slopes look
stable inside fixed-`h` plateaus?

## Source

This is derived from the exact D2 near-ground packets summarized in
`B2G_NEARGROUND_RATIO_GATE_2026-04-29.md`.

## Result

Within each fixed-`h` plateau, the log near-ground ratios decrease roughly
linearly in `n`.

The strongest finite regularity is in `B_2[3]`: the slope of
`ln((h-1)/ground)` stays near `-1.1` to `-1.2` across the observed `h=14`,
`h=15`, and `h=16` plateaus.

`B_2[2]` also compresses inside plateaus, but its observed slope changes more
between the `h=12` and `h=13` plateaus.

| g | h plateau | n range | rows | slope ln((h-1)/ground) / n | slope ln((h-2)/ground) / n | slope ln(ground) / n |
|---:|---:|---|---:|---:|---:|---:|
| 2 | 12 | 30..35 | 6 | -0.739732 | -1.102777 | 1.671374 |
| 2 | 13 | 36..40 | 5 | -1.062657 | -1.438920 | 2.038981 |
| 3 | 14 | 27..30 | 4 | -1.190453 | -1.601916 | 2.346097 |
| 3 | 15 | 31..34 | 4 | -1.191181 | -1.644171 | 2.425124 |
| 3 | 16 | 35..38 | 4 | -1.118923 | -1.630749 | 2.470947 |

## Interpretation

The plateau/reset picture now has two levels:

```text
h-jump -> sparse new ground face -> bath ratio resets high
inside plateau -> ground face expands exponentially -> bath ratio compresses
```

The `g=3` slopes are the cleaner signal. Across three consecutive plateaus,
`ln((h-1)/ground)` compresses at about one natural-log unit per added site.
That is finite evidence for a stable local transfer-operator mechanism, not an
asymptotic theorem.

The `g=2` slopes are still consistent with compression, but two plateaus are
not enough to infer stability.

## Claim Boundary

Safe internal claim:

> Exact finite `B_2[g]` transfer packets show repeated plateau/reset behavior.
> In the observed `B_2[3]` plateaus, normalized near-ground bath ratios compress
> at a roughly stable log-linear rate.

Unsafe:

> The slope is the asymptotic exponent for #755.

Unsafe:

> The slope proves a phase transition.

## Next Gate

Do not push full exact D2 enumeration blindly. The `B_2[3]`, `n=38`, D2 packet
already required `8,217,427,559` second-pass nodes and retained `399,520,098`
near-ground terminal states.

Next engineering move:

```text
add a per-layer capped near-ground diagnostic that records lower bounds and cap-hit flags
```

Then run capped single-row scouts:

```text
B_2[2] layered-capped D2 at n = 42
B_2[3] ground-only or layered-capped D2 at n = 40
```

The capped scout cannot replace exact counts. It only answers whether the next
row is already large enough to preserve the qualitative plateau-compression
direction before paying for full enumeration.
