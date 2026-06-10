# Erdős #1038 — interval certificates for the supremum test-measure inequalities (Tao, 2√2)

Machine-checked interval-arithmetic certificates (directed/outward rounding) for
the three explicit one-variable inequalities — (2.4), (2.5), (2.6) — that
validate the test measures in Terence Tao's December 2025 proof that the Erdős
#1038 supremum equals 2√2. Tao's note verifies these inequalities numerically
("see Figures 2/3/4"); this bundle replaces the plots with rigorous certified
enclosures.

## Contents

| File | What it is |
|---|---|
| `CERTIFIED_RESULTS_REPORT.md` | Certified bounds, domains, the corner argument, honest scope, redaction note |
| `results/supremum_testmeasure_inequalities_DIRECTED_ROUNDING_RESULTS.public.json` | Machine-readable certified results (canonical directed-rounding run; selected implementation-detail fields redacted pending review) |
| `SEALS.md` | SHA-256 seals of the unredacted internal originals |
| `SHA256SUMS` | Hashes of the files in this bundle |

## Scope, in one paragraph

This certifies only the computational step of Tao's already-proven supremum
result. It does not provide Tao's Lemmas 2.1/2.2 (the analysis core), says
nothing about the open infimum in Erdős #1038, and is not a formal (Lean) proof.
The test-measure constants are AlphaEvolve's, per Tao's note. Verification code
is planned for a subsequent version pending review.

## AI acknowledgment

The research direction, problem identification, alternative framings, and
error-catching are the author's work. AI tools were used inside an integrity
environment designed and maintained by the author (independent-checker gates,
fail-closed interval verification, no-overwrite artifact provenance). Within
that environment: Claude (Anthropic, Fable 5 model family, 2026-05 – 2026-06)
for code co-development and certification-pipeline orchestration, alongside
Python toolchains. All results were independently verified by the author. No AI
system is an author.

— Kenneth A. Mendoza
