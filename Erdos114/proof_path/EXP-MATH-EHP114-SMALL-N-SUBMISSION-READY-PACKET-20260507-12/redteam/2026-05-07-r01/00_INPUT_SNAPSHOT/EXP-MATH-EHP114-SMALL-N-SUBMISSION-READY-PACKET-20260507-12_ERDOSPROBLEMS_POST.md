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
