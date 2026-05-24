<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

### 1. Where Models Agree

| Finding | GPT-5.5 Thinking | Claude Opus 4.7 Thinking | Gemini 3.1 Pro Thinking | Evidence |
| :-- | :-- | :-- | :-- | :-- |
| The draft needs tighter *self-contained* definitions (esp. what “length” means) |  | ✓ | ✓ | Both flag missing/implicit definition of lemniscate length `L(·)` and edge behavior at singular level sets.[^1] |
| The `n=13` “zero evaluations” anomaly is a serious credibility/reproducibility risk that must be explained |  | ✓ | ✓ | Both emphasize that “0 branch-and-bound evaluations” could indicate a harness/logic failure, not merely a missing row.[^1] |
| The “Companion Lean 4 source” section being empty should be removed or populated |  | ✓ | ✓ | Both identify the empty Lean code fence as misleading and likely to be flagged.[^1] |
| The Tao result must be cited with a concrete bibliographic pointer (arXiv ID), not just “Tao’s theorem” |  | ✓ | ✓ | Both note the prose alludes to Tao but lacks an explicit citation; Tao has a relevant high-degree preprint.[^2][^3] |
| The post should add missing core prior-art pointers around Problem 114 and lemniscate-length bounds |  | ✓ | ✓ | Both recommend adding key contextual references beyond the one cited paper; Tao’s paper itself surveys prior bounds.[^4][^5] |


***

### 2. Where Models Disagree

| Topic | GPT-5.5 Thinking | Claude Opus 4.7 Thinking | Gemini 3.1 Pro Thinking | Why They Differ |
| :-- | :-- | :-- | :-- | :-- |
| Publish verdict | — | **narrow-and-publish** | **hold-pending-revision** | Claude Opus treats fixes as packaging/clarity issues; Gemini treats the `n=13` anomaly as potentially invalidating nearby rows (esp. `n=14`) until debugged.[^1] |
| How to interpret `n=13` “0 evaluations” | — | Could be bug; ask to rerun and explain; doesn’t necessarily taint other degrees | Likely singularity/NaN/early-exit; could indicate systemic flaw affecting other degrees | Different assumptions about failure modes of interval branch-and-bound and whether they can silently “pass” other degrees.[^1] |
| “Crackpot index” risk level | — | Low (≈ +1), mostly harmless phrasing | Higher (≈ +10), criticizes “SHA-256 sidecars” as techno-bombast | Claude Opus interprets hashes as standard reproducibility; Gemini weights “crypto talk” as signaling-risk unless framed carefully.[^1][^6] |
| Patent-counsel (§101) angle | — | Manageable if you identify a concrete technical contribution; otherwise vulnerable | Broadly vulnerable as “generic math on a computer”; would demand tying to a specific tech application | Claude Opus focuses on narrowing claims / wording; Gemini focuses on §101 “abstract idea” framing and demands application linkage.[^1] |


***

### 3. Unique Discoveries

| Model | Unique Finding | Why It Matters |
| :-- | :-- | :-- |
| Claude Opus 4.7 Thinking | Suggests a precise rewrite: define `L` as 1D Hausdorff measure and explicitly say `deg p = n` (exact) | Prevents reviewer attacks on ambiguity (“degree ≤ n?” “what notion of length?”) and makes the claim formally legible.[^1] |
| Gemini 3.1 Pro Thinking | Hypothesis: `n=13` zero-eval could come from singularities causing NaN/exception propagation in bounding logic | Points to a concrete class of bugs to test (critical points on `|p|=1`, interval extension blowups), not just “rerun it.”[^1] |


***

### 4. Comprehensive Analysis

**High-Confidence Findings.**
Claude Opus 4.7 Thinking and Gemini 3.1 Pro Thinking converge on a core point: the current draft is *not yet self-contained* in the way institutional readers (arXiv moderators, journal referees, or downstream auditors) will expect.[^1] The most important missing piece is a one-line definition of what you mean by lemniscate length `L({|p|=1})`—in practice, specifying a standard notion like the 1-dimensional Hausdorff measure (or arclength of a rectifiable curve, with a brief note about singular level sets) removes ambiguity and preempts corner-case attacks involving cusps/self-intersections when `p'(z)=0` on `|p(z)|=1`.[^1] Closely related, both models want you to clarify that “degree `n`” means `deg p = n` exactly, not `≤ n`, because optimization statements often hide this landmine.[^1]

Both models also agree that the empty “Companion Lean 4 source (excerpt)” block is actively harmful. Even if you are not claiming formal verification, an empty Lean fence reads like “there exists Lean code but it’s omitted,” which is exactly the kind of thing a skeptical reader will pounce on. Removing it entirely is better than leaving it blank; alternatively, include *at least* a theorem statement plus an honest “no formal proof included” note.[^1]

Finally, both models flag citation hygiene: you reference “Tao’s sufficiently-large-`n` theorem” but do not give an arXiv identifier or formal citation. Since Tao’s December 2025 preprint explicitly establishes the conjecture for all sufficiently large `n`, citing it precisely is mandatory. (And it also provides a natural place to cite other historical bounds, because it surveys them in the introduction/table.)[^4][^2][^1]

**Areas of Divergence.**
The real split is the publication verdict. Claude Opus 4.7 Thinking says “narrow-and-publish,” viewing the remaining work as wording, scoping, and presentation fixes that do not undermine the stated finite-`n` claim. Gemini 3.1 Pro Thinking says “hold-pending-revision,” because the `n=13` anomaly (“zero branch-and-bound evaluations”) is not just an omitted row; it is evidence that the computational pipeline can fail in a way you don’t yet understand, which *might* cast doubt on neighboring degrees, including `n=14`.[^1]

This disagreement is about risk tolerance and about how branch-and-bound + interval arithmetic typically fails. If your harness can ever “exit early” due to NaNs, exceptions, or logic short-circuits, then a missing `n=13` explanation can be interpreted as a possible systemic flaw. If, on the other hand, you can show (in logs + a pinned commit + a rerun) that `n=13` is an isolated bug (e.g., an overflow, a wrong parameter-space bound, or a degree-specific special case) and that `n=14`’s run is fully traced and independently reproducible, then Claude Opus’s “publish after tightening” stance becomes much safer.[^1]

There’s also a softer disagreement about “SHA-256 sidecars” and signaling. Claude Opus treats the hashes as standard reproducibility infrastructure and mostly fine. Gemini worries that leading with cryptographic jargon looks like substituting “integrity” for “mathematical correctness,” and recommends reframing hashes as merely pinning artifacts, not as proof. Practically, this is easy to resolve: keep the hashes, but phrase them as “artifact integrity / exact-bitstring pinning,” and separately explain what the certificates *mathematically* certify.[^1]

**Unique Insights Worth Noting.**
Claude Opus’s concrete rewrite recommendation is worth adopting because it resolves multiple reviewer questions at once: define `L` precisely, state `deg p = n`, list the certified degrees as a set `{1,…,12,14}`, and separate the classical (`n=1,2`) cases from the computationally certified (`3..12,14`) ones. Gemini’s speculation that the `n=13` anomaly might be triggered by singularities (critical points on the level set) is also actionable: it suggests exactly what to test and what invariants to log (NaN checks, exception traps, interval bounds at root node, etc.).[^1]

**Recommendations.**
Treat `n=13` debugging as a release blocker *unless* you can provide a short, audit-friendly explanation (root cause + fix or rerun evidence) showing it cannot affect the `n=14` certificate. Then: (i) add a one-clause definition of `L(·)` and specify `deg p = n` exactly, (ii) replace “Tao’s theorem” with an explicit citation to Tao’s arXiv preprint, and (iii) delete the empty Lean block (or populate it with a minimal honest stub).[^2][^1]
<span style="display:none">[^10][^11][^12][^13][^14][^15][^16][^17][^18][^19][^20][^21][^22][^23][^24][^25][^26][^27][^28][^29][^30][^31][^32][^33][^34][^35][^36][^37][^38][^39][^40][^41][^42][^43][^44][^45][^46][^47][^48][^49][^50][^51][^52][^53][^54][^55][^56][^57][^58][^59][^7][^8][^9]</span>

<div align="center">⁂</div>

[^1]: EHP114_v12_REDTEAM_PROMPT.md

[^2]: https://arxiv.org/abs/2512.12455

[^3]: https://terrytao.wordpress.com/2025/12/15/the-maximal-length-of-the-erdos-herzog-piranian-lemniscate-in-high-degree/

[^4]: https://arxiv.org/pdf/2512.12455.pdf

[^5]: https://arxiv.org/abs/0805.2295

[^6]: https://www.sciencedirect.com/science/article/abs/pii/S0022247X24004931

[^7]: https://www.erdosproblems.com/114

[^8]: https://www.erdosproblems.com/forum/thread/114

[^9]: https://arxiv.org/html/2512.12455v1

[^10]: https://www.reddit.com/r/askmath/comments/1qn9akb/can_anyone_explain_why_the_problem_of_the_maximum/

[^11]: https://arxiv.org/html/0805.2295v2

[^12]: https://www.math.purdue.edu/~eremenko/dvi/erdos23.pdf

[^13]: https://www.erdosproblems.com/latex/114

[^14]: https://arxiv.org/pdf/2407.14610.pdf

[^15]: https://terrytao.wordpress.com/tag/erdos/

[^16]: https://arxiv-math.livejournal.com/25709605.html

[^17]: https://arxiv.org/pdf/0805.2295.pdf

[^18]: https://arxiv.org/abs/0808.0717

[^19]: https://ui.adsabs.harvard.edu/abs/2008arXiv0808.0717F/abstract

[^20]: https://ui.adsabs.harvard.edu/abs/arXiv:2512.12455

[^21]: https://www.jstor.org/stable/2160803

[^22]: https://projecteuclid.org/download/pdf_1/euclid.mmj/1028998227

[^23]: https://app.icerm.brown.edu/materials/Slides/sp-f16-w2/Random_lemniscates_%5D_Erik_Lundberg,_Florida_Atlantic_University.pdf

[^24]: https://www.ams.org/proc/1995-123-03/S0002-9939-1995-1223265-3/S0002-9939-1995-1223265-3.pdf

[^25]: https://projecteuclid.org/journals/michigan-mathematical-journal/volume-46/issue-2/On-the-length-of-lemniscates/10.1307/mmj/1030132418.pdf

[^26]: https://terrytao.wordpress.com/2025/12/

[^27]: https://archive.org/download/commencement19531953univ/commencement19531953univ.pdf

[^28]: https://digitalcommons.conncoll.edu/context/alumnews/article/1108/viewcontent/CCAlumnaeNews_December_1953.pdf

[^29]: https://scholarship.law.cornell.edu/cgi/viewcontent.cgi?article=3849\&context=clr

[^30]: https://arxiv.org/pdf/2503.18270.pdf

[^31]: https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Novikov-Cobordism.pdf

[^32]: https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli

[^33]: https://www.ams.org/books/surv/003/surv003-endmatter.pdf

[^34]: https://www.math.purdue.edu/~eremenko/dvi/lempert.pdf

[^35]: https://mathcurve.com/courbes2d.gb/lemniscate/lemniscate.shtml

[^36]: https://www.nationalacademies.org/read/11540/chapter/23

[^37]: https://www.sciencedirect.com/science/article/pii/S2405844024101260

[^38]: https://csclub.uwaterloo.ca/~pbarfuss/mckean-moll2.pdf

[^39]: https://web.stanford.edu/~lindrew/18.218-2.pdf

[^40]: https://github.com/przchojecki/agentic-erdos

[^41]: https://en.wikipedia.org/wiki/Polynomial_lemniscate

[^42]: http://e.math.hr/aggregator/categories/1

[^43]: https://www.semanticscholar.org/paper/THE-ARC-LENGTH-OF-THE-LEMNISCATE-|-w-n-+-c-|-=-1-Wang-Peng/70131e5fdf9ab4de4fc2f5f518df1b516fce82fe

[^44]: https://www.erdosproblems.com

[^45]: https://b.hatena.ne.jp/entry/s/terrytao.wordpress.com/2009/01/09/245b-notes-3-lp-spaces/

[^46]: https://mathstodon.xyz/@kenmendoza

[^47]: https://zenodo.org/records/10034352

[^48]: https://zenodo.org/records/16946279

[^49]: https://zenodo.org/records/7065654

[^50]: https://zenodo.org/records/15084347

[^51]: https://github.com/topics/combinatorics?o=asc\&s=forks

[^52]: https://github.com/erdos-project/erdos-experiments

[^53]: https://zenodo.org/records/8358829

[^54]: https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erdős-problems

[^55]: https://zenodo.org/records/15080343

[^56]: https://mendozalab.io

[^57]: https://github.com/neelsomani/gpt-erdos

[^58]: https://zenodo.org/records/13329929

[^59]: https://news.ycombinator.com/item?id=46560445

