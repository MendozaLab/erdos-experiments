# EXP-MM-030-PMF-DIFFSET-SKELETON-WINDOW-56-58-2026-04-30

Date: 2026-04-30
Problem: Erdos #30
Status: DERIVED_SKELETON_DIAGNOSTIC / EXACT_PACKET_BACKED / INTERPRETIVE

## Source

Derived from `EXP-MM-030-PMF-MAXWELL-FACE-HANDOFF-56-58-V2-2026-04-30`.

## Question

Does the `n = 57` face handoff correspond to a split in the positive-difference skeleton, or only to a field/embedding switch inside the same skeleton?

## Result

The `n = 57` handoff is not a difference-set skeleton split. All six exported exact maximizers at `n = 57` share one positive-difference skeleton. Their site sets can be far apart, but their difference-memory content is identical.

| n | exact maximizers | unique difference skeletons | field split | max site distance | max diffset distance |
|---:|---:|---:|---|---:|---:|
| 56 | 4 | 1 | false | 18 | 0 |
| 57 | 6 | 1 | true | 18 | 0 |
| 58 | 10 | 2 | true | 20 | 16 |

## Interpretation

The Maxwell Face-Handoff is narrower than the previous structural guess. At `n = 57`, the missing information is not a new difference-memory skeleton. It is the full embedding/field-response geometry of a single skeleton. In this window, a second difference skeleton first appears at `n = 58`, after the handoff has already occurred.

## Claim Boundary

This is a finite derived diagnostic. It strengthens the mechanism by falsifying an over-broad reading: the `57` handoff is field/embedding selection inside one skeleton, not a proved spectral phase transition or theorem-level Sidon result.
