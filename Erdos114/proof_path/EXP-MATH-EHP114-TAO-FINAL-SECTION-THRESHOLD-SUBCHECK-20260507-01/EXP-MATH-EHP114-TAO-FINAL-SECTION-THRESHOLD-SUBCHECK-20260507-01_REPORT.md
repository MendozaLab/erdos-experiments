# EHP114 L40B Final-Section Threshold Subcheck

Experiment: `EXP-MATH-EHP114-TAO-FINAL-SECTION-THRESHOLD-SUBCHECK-20260507-01`

## Verdict

Status: `FINAL_SECTION_REMAINS_OPAQUE`.

The subcheck only computes row-level thresholds for dependencies whose constants are already explicit. Since the selected final-section rows still have named blockers, no partial maximum threshold is emitted.

## Claim Ceiling

Final-section threshold subcheck only. It emits no global candidate_N0 and authorizes no higher-degree computation.
