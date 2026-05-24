# SOTA Crosswalk - 2026-05-02

## Sources Checked

- Tao/erdosproblems AI contributions tracker: https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erd%C5%91s-problems
- Aristotle #728 writeup: https://arxiv.org/abs/2601.07421 and https://arxiv.org/html/2601.07421
- Alexeev/Putterman/Sawhney/Sellke/Valiant I: https://arxiv.org/abs/2603.29961
- Alexeev/Putterman/Sawhney/Sellke/Valiant II: https://arxiv.org/abs/2604.06609
- Harmonic Aristotle report: https://harmonic.fun/pdf/Aristotle_IMO_Level_Automated_Theorem_Proving.pdf

## Current Public Facts

The current external field is stronger than the Atlas on theorem resolution. The Tao tracker lists many AI contributions to Erdos problems, including full Lean solutions, partial results, incorrect proofs, and literature-dependent cases. That matters because the comparison set is not hypothetical anymore.

The #728 writeup states that the result was regarded as the first Erdos problem fully resolved autonomously by an AI system, using GPT-5.2 Pro and Harmonic's Aristotle, with a final Lean proof translated into informal mathematics. It also records Boris Alexeev helping simplify the proof. This is direct solver/formalizer performance, not just infrastructure.

The Alexeev/Putterman/Sawhney/Sellke/Valiant March 31, 2026 paper gives three short proofs answering Erdos questions and says each proof was due entirely to an internal OpenAI model. The April 8, 2026 sequel gives five more proofs from questions posed by Erdos, again attributed to an internal OpenAI model.

The Harmonic Aristotle report states a high bar for solved status: a complete Lean 4 proof against Mathlib, without gaps or unsound axioms such as `sorryAx`.

## Dimension Comparison

| Dimension | Alexeev / Aristotle / AI-Erdos ecosystem | ErdosAtlas current audit |
|---|---|---|
| Theorem resolution | Stronger. Public record includes full solutions and Lean formalizations. | Weaker. No claim that #30 is solved/proved/closed. |
| Formal verification | Stronger on autonomous Lean resolution in selected cases. | Real local Lean portfolio, but many artifacts are formalization/certificate/infrastructure and some use axioms. |
| Systematic deployment | Strong: repeated attempts across Erdos problems with public tracker and papers. | Partial: local map-check-formalize pipeline exists, but public-safe framework and final public packet are not yet shipped. |
| Reproducibility | Mixed but public: arXiv, Lean files, tracker entries. | Strong locally where SHA sidecars exist; not yet public enough. |
| Corpus coverage | External ecosystem is broad in attempts/results. | Atlas has broad indexing and 14,095 edges via live API, but method internals are not public. |
| Public artifacts | Strong: arXiv papers, tracker, Lean artifacts. | Partial: live API/workbench exists, but method-safe public package is not yet filed/published. |
| Platform differentiator | Solver/formalizer pipeline. | Map-check-formalize instrument with explicit falsification/demotion and evidence manifests. |

## Verdict On The Phrase

`Boris Alexeev-style systematic deployment`: PARTIAL.

Supported only if the phrase means systematic, repeated deployment across an Erdos corpus surface. It is an overclaim if it implies comparable autonomous theorem-resolution capability. The safe comparison is dimensional:

> ErdosAtlas is not comparable to Aristotle/Alexeev-style work on autonomous theorem resolution. It is comparable only as a systematic deployment instrument, and its differentiator is the map-check-formalize layer with explicit demotions and evidence manifests.

## Strongest Safe Rewrite

> ErdosAtlas is not yet a theorem-resolution platform at the level of autonomous solvers. Its defensible lane is a systematic deployment instrument: it maps putative structural correspondences, tests them with reproducible finite packets, demotes failures, and routes survivors toward formal proof targets.
