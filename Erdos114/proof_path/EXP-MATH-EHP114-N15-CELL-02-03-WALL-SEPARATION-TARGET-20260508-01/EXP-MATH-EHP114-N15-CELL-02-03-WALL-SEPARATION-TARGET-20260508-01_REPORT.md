# EXP-MATH-EHP114-N15-CELL-02-03-WALL-SEPARATION-TARGET-20260508-01

## Meaning

This packet turns the latest `CELL-02-03` adaptive slice into a narrow target:
the blocker is now wall separation, and the wall inequality is dominated by the
center-strip absolute bound rather than by the third-order collar remainder.

This is review-only. It is not a degree-15 certification and it does not change
any existing Erdős #114 packet.

## Source

- Source experiment: `EXP-MATH-EHP114-N15-BOUNDARY-SLICE-CELL-02-03-ADAPTIVE-REFINE-20260508-03`
- Source SHA-256: `5f4e6cadcee4133347f27755c026bc1bbf43411600b3e52b4d9ed52e59ae93bf`
- Source verdict: `BOUNDARY_SLICE_PARTIAL_NEEDS_REDUCTION`
- Remaining unresolved reasons: `{'third_order_critical_exclusion_failed': 69, 'third_order_wall_separation_failed': 108}`

## Wall-Separation Diagnosis

- Wall failures: `108`
- Critical-exclusion failures still present: `69`
- Dominant mode: `CENTER_STRIP_DOMINATED_WALL_SEPARATION`
- Required lhs factor, median: `0.008720274680208527`
- Center-strip share of lhs, median: `0.9988204069392397`
- Normal-remainder share of lhs, median: `0.0011795930607602328`

The practical read is simple: deeper interval subdivision helped move the first
blocker from critical exclusion to wall separation, but the present wall test is
too conservative because it pays for `sup_r |F(0,r)|` on the center strip. The
normal remainder is small. The next useful experiment is therefore not broad
compute; it is a signed/cancellation-aware wall test.

## Top Owner Groups

- `4488:3`: 16 failures, max gap 3.2771, median required factor 0.0129025
- `4571:2`: 11 failures, max gap 1.6483, median required factor 0.00390399
- `2484:4`: 16 failures, max gap 1.04721, median required factor 0.0287227

## Target Inequalities

1. `WS-01-CENTER-STRIP-CANCELLATION`: replace the raw absolute center-strip
   bound with a certified signed or cancellation-aware bound.
2. `WS-02-SIGNED-WALL-ENDPOINTS`: record signed normal-wall endpoint intervals
   for the current wall-failure boxes.
3. `WS-03-DOMINANT-OWNER-LOCAL-MODEL`: explain the largest-gap owner group
   before broadening this method to other cells.

## Recommended Next Experiments

- `EXP-MATH-EHP114-N15-CELL-02-03-CENTERLINE-SIGN-MODEL-20260508-01`: Replace raw center-strip absolute bound with a signed or cancellation-aware centerline model.
- `EXP-MATH-EHP114-N15-CELL-02-03-WALL-ENDPOINT-SPLIT-20260508-01`: Record direct signed endpoint separation on normal walls for the 108 current wall failures.
- `EXP-MATH-EHP114-N15-CELL-02-03-OWNER-LOCAL-MODEL-20260508-01`: Analyze the worst source_ownership_key group before broadening to CELL-07-03 or CELL-06-03.

## Safety

Forbidden writes remain D1, curated morphisms, public pages, proof registries,
and existing #114 finite packets. This artifact is an auto-research target
packet only.
