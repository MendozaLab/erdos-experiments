# Interval-arithmetic certificates for the three test-measure inequalities in Tao's Erdős #1038 supremum proof

## What this is

Terence Tao's December 2025 note *Sublevel Sets of Logarithmic Potentials* proves
that the supremum in Erdős problem #1038 equals 2√2 (Theorem 2.1). After two
potential-theory lemmas (2.1 duality, 2.2 rearrangement), the proof reduces to
three explicit one-variable inequalities — (2.4), (2.5), (2.6) — validating the
test measures in each of three cases; the note verifies these numerically
("see Figures 2/3/4"). This bundle replaces those numerical verifications with
machine-checked certificates computed in interval arithmetic with directed
(outward) rounding, so every certified bound is a rigorous enclosure rather than
a floating-point evaluation.

The test-measure constants are due to AlphaEvolve, as reported in Tao's note.
What this bundle adds is the rigor: a certificate that each inequality genuinely
holds on its whole domain.

## Certified results (directed-rounding run — canonical)

| Inequality | Domain | Certified bound |
|---|---|---|
| (2.6) Case 3 | [0.7624, 2.7987] (full domain) | max LHS ≤ −1.37×10⁻⁴ |
| (2.5) Case 2 | [0.7987, 4], tail t ≥ 4 closed by sign + monotonicity (full domain) | max LHS ≤ −1.38×10⁻⁸ |
| (2.4) Case 1, interior | [√2+10⁻³, 1.7624] | min D ≥ +5.57×10⁻¹¹ |
| (2.4) Case 1, corner | [√2, √2+10⁻³] | D″ ∈ [0.6698, 0.8045] ⇒ D″ > 0; with the exact identities D(√2) = D′(√2) = 0 (from (√2+1)(√2−1) = 1), Taylor's theorem closes (√2, √2+10⁻³] |

The margins are tight — a quadratic touch at √2, and order 10⁻⁸ near t ≈ 2.47 —
which is exactly why plots are not proofs at this step and certified computation
earns its keep.

Machine-readable record: `results/supremum_testmeasure_inequalities_DIRECTED_ROUNDING_RESULTS.public.json`.

## Honest scope (read this before citing)

This certifies the **computational inequalities only** — the step Tao's note
verifies by plot. It does **not** provide Tao's Lemmas 2.1 / 2.2 (the analysis
core, and the bulk of any Lean formalization of
`FormalConjectures/ErdosProblems/1038.lean parts.ii`). The supremum value 2√2 is
**Tao's** result, already proven; this bundle adds machine-checked rigor to one
step of it. No claim is made about the open infimum, no SOTA advance, no #1038
resolution, and this is not a Lean proof.

## Redactions

Selected implementation-detail fields (certification-method internals,
subdivision granularity, code references) are redacted or removed from the
public results file pending IP review; every certified bound, domain, and the
corner mathematics are unmodified. SHA-256 seals of the unredacted internal
originals are restated in `SEALS.md`. Verification code is planned for a
subsequent version pending review.
