# Draft erdosproblems.com post for Problem 114

> Target length: ~300 words of prose. No headers in body. Single question at end. Voice: tired-engineer-at-11pm, not conference-abstract.

---

I have a small-degree certificate I'd like to record on this page, with the explicit caveat that it isn't an attempt at the all-degree statement. It is meant as the finite frontier complementary to Tao's sufficiently-large-`n` theorem (arXiv:2512.12455).

Finite statement. For every monic complex polynomial `p` with `deg p = n` exactly and `n` in the contiguous range `1 ≤ n ≤ 14`, the one-dimensional Hausdorff measure of `{z : |p(z)| = 1}` is at most the same measure of `{z : |z^n - 1| = 1}`.

The dependency split is typed. `n = 1` is the translated unit circle. `n = 2` is the MacLane (1953/54) degree-two case, with Eremenko–Hayman's "On the length of lemniscates" (Michigan Math. J. 46 (1999), 409–415; arXiv:0805.2295) as the accessible source: their abstract names Bernoulli's lemniscate as the `d = 2` extremal level set, and the proof remark after Lemma 5 identifies `z^2 + 1` as extremal, length-equivalent to `z^2 - 1` by rotation. The block `3 ≤ n ≤ 14` is covered by reproducible Rust + `inari` IEEE-1788 interval certificates with row-level artifact-integrity hashes (SHA-256). The certificate set is archived under concept DOI 10.5281/zenodo.19184467 (resolves to latest version); source at github.com/MendozaLab/erdos-experiments, release `v3.1.0`.

Worth a sentence on `n = 13`. An earlier batch produced a row with zero branch-and-bound evaluations because the BB box-construction step returned an empty vector at high reduced dimension and the verdict code mistook that for completion. The verdict-logic bug was patched at the source. The `n = 13` row was re-run with the patched binary on 2026-05-07, closed at level 0 after 197M box evaluations and 53 minutes of wall-clock, and is now part of the certificate set.

I'm not extracting a threshold from Tao's theorem here, and I'm not claiming any `n ≥ 15`.

Is this the right form in which to record the finite side on this page?

---

The research direction, problem identification, alternative framings, and error-catching are the author's work. AI tools were used inside an author-architected integrity environment (formal-verification gates, Mathlib API correctness checks, verdict-logic guards, sorry-detection, axiom audit, prepub-redteam audit folder). Within that environment: Claude (Anthropic, claude-opus-4-7) for code co-development and exploratory analysis; Perplexity (web app, quorum mode: claude-opus-4-7 + gemini-3-pro-deep-think + gpt-5-pro) for adversarial pre-publication review. All mathematical results were independently verified by the author against the cited literature and SHA-checked certificate artifacts.
