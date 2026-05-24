# Pre-Publication Red-Team Request — Draft erdosproblems.com post for Problem 114

**Round:** 2026-05-07-r01
**Generated:** 2026-05-07T12:28:05
**Target audience:** all

I am preparing a preprint for external deposit (Zenodo + arXiv) and want an adversarial red-team critique BEFORE I freeze the artifacts. Please act as an independent reviewer with no prior context on this work. The artifact is reproduced below in full so you can critique without external lookup.

---

## Title

Draft erdosproblems.com post for Problem 114

## Abstract

> Recommended length: under 300 words of prose. No bullet stacks, no inline tables, no bold headers in the body. One question at the end.

## Body excerpt

```
# Draft erdosproblems.com post for Problem 114

> Recommended length: under 300 words of prose. No bullet stacks, no inline tables, no bold headers in the body. One question at the end.

---

I have a small-degree certificate I'd like to record on this page, with the explicit caveat that it is not an attempt at the all-degree statement and is meant as the finite frontier complementary to Tao's sufficiently-large-`n` theorem.

The finite statement: for every monic complex polynomial `p` of degree `n` with `1 <= n <= 12` and `n = 14`, the lemniscate length `L({|p|=1})` is at most the lemniscate length of `z^n - 1`. The dependency split is typed. `n = 1` is the translated unit circle. `n = 2` is the MacLane (1953/54) degree-two case, with Eremenko–Hayman's "On the length of lemniscates" (Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295) as the accessible source: their abstract names the Bernoulli lemniscate as the `d = 2` extremal level set, and the proof remark after Lemma 5 identifies `z^2 + 1` as extremal, length-equivalent to `z^2 - 1` by rotation. The remaining degrees `3..12` and `14` are covered by reproducible Rust + `inari` IEEE-1788 interval certificates with row-level SHA-256 sidecars; the certificate set is archived under DOI 10.5281/zenodo.19480329 and the source repository is github.com/MendozaLab/erdos-experiments.

The `n = 13` row exists in the certificate archive but is quarantined from the public table for this submission because the underlying run reports zero branch-and-bound evaluations while neighboring rows report tens to hundreds of millions. I'd rather flag that and re-run the row than ship a non-uniform table; treating `n = 13` as a TODO seemed more honest than burying it.

I am not extracting a threshold from Tao's theorem here, and I am not claiming any degree `n >= 15`.

I'd appreciate guidance on whether this is the right form in which to record the finite side on this page.

---

Tooling disclosure: this draft was prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature and SHA-checked certificate artifacts.

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
