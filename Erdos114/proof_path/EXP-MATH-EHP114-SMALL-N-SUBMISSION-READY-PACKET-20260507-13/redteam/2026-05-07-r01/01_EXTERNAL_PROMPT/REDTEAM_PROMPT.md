# Pre-Publication Red-Team Request — Draft erdosproblems.com post for Problem 114

**Round:** 2026-05-07-r01
**Generated:** 2026-05-07T21:20:47
**Target audience:** all

I am preparing a preprint for external deposit (Zenodo + arXiv) and want an adversarial red-team critique BEFORE I freeze the artifacts. Please act as an independent reviewer with no prior context on this work. The artifact is reproduced below in full so you can critique without external lookup.

---

## Title

Draft erdosproblems.com post for Problem 114

## Abstract

> Target length: ~300 words of prose. No headers in body. Single question at end. Voice: tired-engineer-at-11pm, not conference-abstract.

## Body excerpt

```
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

```

## Companion Lean 4 source (excerpt)

```lean

```

---

## Red-team requests — answer ALL eight points

Provide a structured critique. Be specific, cite exact lines, propose verbatim rewrites where applicable.

### 1. Counterexamples or corner cases
Does the main result fail for any concrete instance the proof misses? Be specific. Consider degenerate cases (empty inputs, rank-deficient inputs, numerical-tolerance edge cases). If you find a counterexample, state it as `n = ?, C = ?, expected = ?, observed = ?`.

### 2. Implicit assumptions
Does the proof depend on hypotheses not stated explicitly? List each one with the exact line of the proof where it's implicitly used. Suggest the literal hypothesis text to add.

### 3. Crackpot-Index inflators (Baez)
Score the writing against Baez's published Crackpot Index (https://math.ucr.edu/home/baez/crackpot.html). Flag any phrasing that scores >0 with the exact line and a verbatim suggested rewrite. Common failure modes: grandiose claims, unprovable assertions about novelty, NotebookLM-flavored title bombast, "I have proved" without "kernel-checked".

### 4. Missing prior art
Is there a published paper, theorem, or preprint that anticipates or strictly subsumes the result? Cite by DOI, arXiv ID, or full bibliographic reference. Distinguish anticipating prior art (problematic) from supporting context (worth citing).

### 5. Adversarial-counsel critique
If you were a senior patent counsel trying to invalidate any patent claim that cites this preprint as a §101 anchor, what specific weakness would you attack? What narrowing would you demand? Quote the exact preprint sentence(s) you'd target.

### 6. Suggested narrowings or strengthenings
Quote the exact sentence(s) that should be edited. Provide verbatim replacement text. Distinguish narrowings (defensive, reduce overclaim risk) from strengthenings (offensive, sharpen the result).

### 7. Formal-verification audit
If a Lean source is included above:
- Identify any `axiom` declaration whose citation does not match the prose claim
- Identify any theorem signature that does not match the prose claim
- Flag any `sorry` or unsupported tactic that hides a substantive gap
- Note any Mathlib API name that may have been renamed in current `master`

### 8. Plain-English bottom line
One paragraph (~150 words) the inventor can act on at 11pm. Lead with the verdict in three categories:
- **publish-as-is** (only if zero meaningful findings in 1-7)
- **narrow-and-publish** (specific edits required, but the result holds)
- **hold-pending-revision** (substantive issue requires re-derivation or re-proof)

State the verdict, then the top 1-3 specific actions.

---

Be rigorous and adversarial. This critique becomes part of the permanent audit trail of the publication and will be reviewed by downstream institutional readers (Cooley patent counsel, NIH program officers, FDA reviewers, Mathlib reviewers, future-self diligence).

Format your response as a structured Markdown document with one section per numbered request above. The inventor will paste your response verbatim into the audit folder; preserve your formatting.
