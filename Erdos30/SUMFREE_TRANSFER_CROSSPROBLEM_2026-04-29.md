# Sum-Free Transfer Cross-Problem Gate

Date: 2026-04-29
Problem: Erdős #166 candidate lane
Reference problem: Erdős #30 PMF transfer operator
Status: API_PORT_PASS, HANDOFF_SIGNATURE_NOT_SEEN

## Question

Is the field-sensitive handoff a Sidon-only artifact, or does the same PMF
state-machine signature appear under a changed local exclusion rule?

## Answer

The transfer-state API ports cleanly. The field-sensitive handoff signature does
not appear in this first sum-free ground-face test.

That is a useful negative result. It means the PMF lane is not just producing
the same story everywhere. Under the sum-free exclusion rule, the exact
maximizer face is tiny and rigid, and the prefix, mass, and joint fields all
select the same winner throughout the tested window.

## Operator Port

The new binary lives in the same Rust crate:

```text
erdos-experiments/Erdos30/rust-transfer-operator/src/bin/sumfree_transfer.rs
```

State:

```text
State = {
  occupied_mask: u128,
  cardinality: u8
}
```

Changed local exclusion rule:

```text
occupy x iff no occupied a,b satisfy a + b = x
```

The scan uses the positive lattice `1..=n`. The parity reference is the
classical finite maximum size:

```text
h(n) = ceil(n/2)
```

## Packets

- `EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-20-40-2026-04-29`
- `EXP-MM-166-PMF-SUMFREE-TRANSFER-D0-41-60-2026-04-29`

Both SHA256 sidecars verified `OK`.

## Result Summary

Combined window:

```text
20 <= n <= 60
checked rows: 41
h(n)=ceil(n/2) formula matches: 41 / 41
frontier split count: 0
combined runtime: 3.566481459 sec
```

Selected rows:

| n | h(n) | ground degeneracy | entropy ln | terminal retained | pruned states | joint winner |
|---|---:|---:|---:|---:|---:|---|
| 20 | 10 | 3 | 1.098612 | 3 | 1,241 | `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19]` |
| 40 | 20 | 3 | 1.098612 | 3 | 186,075 | `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39]` |
| 60 | 30 | 3 | 1.098612 | 3 | 19,472,383 | `[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59]` |

## Interpretation

The API generalizes; the Sidon handoff does not automatically generalize.

For Sidon #30, the exact maximizer face is large enough for different fields to
select different witnesses. For sum-free #166 in this finite model, the ground
face is essentially rigid: the same odd-set witness wins prefix, mass, and joint
fields across `20 <= n <= 60`.

That makes the PMF story more credible, not less. A real collider should produce
different phase behavior under different local exclusion laws.

## Claim Boundary

Safe internal claim:

> The PMF transfer-state API ports from Sidon to a changed additive exclusion
> rule, and the first #166 ground-face scan cleanly distinguishes rigid
> sum-free faces from Sidon field-sensitive handoffs.

Unsafe:

> The #166 scan proves the Sidon PMF theorem.

Unsafe:

> Field-sensitive handoffs are universal across additive combinatorics.

## Next Gate

The better next cross-problem target is #755 `B_h[g]`, not deeper #166. It is
closer to Sidon because it changes "no repeated sums/differences" into "bounded
repeated sums" rather than collapsing to the rigid odd/upper-half extremal
families.

