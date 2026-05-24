# EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-W3N8-20260505-02 Report

## Status

Rust per-core closure runner for Erdos #20. This is an internal diagnostic
artifact. It is not theorem progress, not lower-bound progress, and not a
Leg-4 pass. The claim ceiling remains: this is a shadow signature, not
universal law.

## Verdict

- Runner status: `PASS`
- Classification: `RUST_PER_CORE_SIGNAL_PRESENT`
- Rust engine: `erdos-experiments/Erdos20/rust_core_closure`

The runner confirms that Rust is now useful for this lane. Python was enough to
define the observable; Rust is the right layer for larger per-core sweeps.

## Executed Targets

| w | n | families/status | target m | aggregate I | strongest s | strongest local I | seconds |
|---|---|---:|---:|---:|---:|---:|---:|
| 3 | 8 | 148790380 | 8 | 5.87333 | 1 | 0.393312 | 29.2811 |

## Interpretation

The important comparison is not just the aggregate closure pressure. The Rust
runner records fixed-core channels at the selected target size, so the question
becomes whether particular core sizes carry a repeatable local closure cost.

This still does not execute the full Leg-4 test. The missing pieces are a
predeclared floor-normalized numerator, Abbott-Hansen-Sauer baseline controls,
and a symmetry-reduced transfer formulation that can push beyond exact
enumeration.

## Next Rust Target

The next queued run is `w=3,n=8`; it has about `148790380` sunflower-free
families in the saved aggregate artifact. The runner is ready, but that run
should be treated as a heavier execution packet rather than mixed into this
calibration artifact.

## Source Boundary

No scorecard, D1, public document, git, email, CLAUDE.md, or AGENTS.md was
changed.
