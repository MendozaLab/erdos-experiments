# Zenodo strategy for the v13 submission packet

## Decision: v3.1.0 version published manually

A certificate row changed (n=13 was re-run on 2026-05-07 with the verdict-bug-patched binary and the new artifact landed in `results/erdos-114/`). Per the v12 strategy memo's own conditions, this is exactly the trigger for a Zenodo v6 mint:

> Mint a new Zenodo version (v6) only if any of the following hold:
> - The `n = 13` row is re-run with a populated branch-and-bound trace and the new JSON+SHA replaces the existing v5 row. ✓ (this happened 2026-05-07)
> - A new degree (`n = 15` or higher) is added.
> - The certificate format changes (e.g., new `inari` version, new schema fields).
> - A row's `l_star_lower` or `l_star_upper` changes.

The trigger condition is met. A new version is appropriate.

## What happened

Per `CLAUDE.md` § "GitHub → Zenodo Auto-Archive (erdos-experiments repo)":

> **MendozaLab/erdos-experiments is connected and green on Zenodo. Every GitHub Release automatically triggers a Zenodo snapshot and mints a versioned DOI.**

The 2026-05-07 GitHub release [`v3.1.0`](https://github.com/MendozaLab/erdos-experiments/releases/tag/v3.1.0) ("EHP n=13 — corrected B&B certificate") did not appear publicly on Zenodo after re-check. On 2026-05-08, a new version was manually created from latest record `19480329`, preserving the existing concept DOI family and publishing version DOI `10.5281/zenodo.20087919`. Downstream artifacts that need an immutable certificate-set pin should cite the version DOI; broader references can cite concept DOI `10.5281/zenodo.19184467`, which resolves to latest.

## Description used for v3.1.0

The published version description records:

```
Reproducible IEEE-1788 interval certificates for the lemniscate length
of monic complex polynomials of degree 3 ≤ n ≤ 14, computed in Rust
using the `inari` IEEE-1788 interval-arithmetic crate, in the context of
the Erdős–Herzog–Piranian conjecture (Erdős Problem 114).

Each n is certified by a separate row in
results/erdos-114/EXP-MM-EHP-007-n{N}-inari_RESULTS.json with a
SHA-256 sidecar at the same path with the .sha256 extension.

Companion small-degree variant (1 ≤ n ≤ 14) records the dependency split:

  - n = 1: direct analytic, translated unit circle.
  - n = 2: literature row, citing G. R. MacLane, "On a conjecture of
    Erdős, Herzog, and Piranian," Michigan Math. J. 2 (1953/54),
    147-148, doi:10.1307/mmj/1028989918, and Alexandre Eremenko and
    Walter K. Hayman, "On the length of lemniscates," Michigan Math.
    J. 46 (1999), 409-415; arXiv:0805.2295.
  - 3 ≤ n ≤ 14: the certificate JSONs in this record.

Note on n = 13: the v5 release (2026-04-09) carried an n=13 row with
bb_total_evals = 0 due to a verdict-logic bug at high reduced dimension
(the engine could not distinguish "BB exhaustively eliminated" from
"BB had nothing to evaluate"). The bug was patched in commit dae62b8
of MendozaLab/erdos-experiments and the n=13 row was re-run on
2026-05-07 (commit f89597e), producing 197,132,288 box evaluations and
closing the proof at level 0 in 52.96 minutes wall-clock on a 32-CPU
Modal worker. The new n=13 artifact replaces the prior version in this
record. Released as v3.1.0.

The all-degree Erdős #114 conjecture is not asserted by this record,
and no threshold is extracted from Tao's sufficiently-large-n theorem
(arXiv:2512.12455).

The research direction, problem identification, alternative framings,
and error-catching are the author's work. AI tools were used inside an
author-architected integrity environment (formal-verification gates,
Mathlib API correctness checks, verdict-logic guards, sorry-detection,
axiom audit, prepub-redteam audit folder). Within that environment:
Claude (Anthropic, claude-opus-4-7) for code co-development and
exploratory analysis; Perplexity (web app, quorum mode: claude-opus-4-7
+ gemini-3-pro-deep-think + gpt-5-pro) for adversarial pre-publication
review. All mathematical results were independently verified by the
author against the cited literature and SHA-checked certificate
artifacts.

Source repository: https://github.com/MendozaLab/erdos-experiments
GitHub release v3.1.0: https://github.com/MendozaLab/erdos-experiments/releases/tag/v3.1.0
Erdős Problem #114 reference: https://www.erdosproblems.com/114
```

## Concept DOI vs versioned DOI

- Version DOI for corrected v3.1.0 packet: [10.5281/zenodo.20087919](https://doi.org/10.5281/zenodo.20087919)
- Concept DOI (always latest): [10.5281/zenodo.19184467](https://doi.org/10.5281/zenodo.19184467)

External citations should target the concept DOI when long-term traceability is wanted and the versioned DOI when the audit must pin to a specific certificate set. The DeepMind formal-conjectures PR #3712's existing citation of v5 DOI 10.5281/zenodo.19480329 is **not invalidated** by the v3.1.0 mint — versioned DOIs are stable. The new finite-variant PR (this packet's `_FORMAL_CONJECTURES_PACKET.md`) should cite the latest versioned DOI 10.5281/zenodo.20087919 when it opens or is updated.

## Process gate

Before any new Zenodo description edit:

1. **prepub-redteam round on the v13 frozen post + Zenodo description text** — the description text counts as a public-facing artifact. New round (`init … --new-round`), import the Perplexity quorum critique, freeze, verify.
2. **Publisher Gate** over the description text.
3. **Crackpot-Scrub** over the description text.

The auto-minted v6 description (if it was minted from the GitHub release notes) does not need editing if the release notes already carry the Block-G AI-acknowledgment per the other agent's session report.

## What does not change

- DOIs already published (v5 versioned, v4 versioned, etc.) remain stable.
- The certificate JSONs at `results/erdos-114/EXP-MM-EHP-007-n{3..12,14}-inari_RESULTS.json` are unchanged from v5.
- The DeepMind formal-conjectures PR #3712 citation of v5 DOI is intact.
- The version DOI 10.5281/zenodo.20087919 pins the corrected v3.1.0 packet.
- The concept DOI 10.5281/zenodo.19184467 always points to latest.
