# Round 22 — Methodology Notes

Mode 1 dispatches don't produce computational outputs to audit — there is no sweep code, no numerical results, no reproducibility check. The methodology questions here are about the quality of PC's literature selection, morphism design, and schema engineering rather than numerical correctness. On those axes, Round 22 is mostly sound, with two targeted follow-up notes.

## §1 — Literature scout quality

The 24-citation bibliography is well-calibrated for the endpoint-limit gate. The coverage is appropriate across four distinct zones:

**The hyperelliptic period-theory core (Part A)** is correctly anchored by Deift et al. 2001, Kuijlaars-Mo 2009, Frauendiener-Klein 2015. These are the right foundational papers; PC correctly notes that Frauendiener-Klein was already in the assimilation record (arXiv:1408.2201, cited in Round 19) and treats the 2015 journal version as the canonical reference. The addition of Wolper 2021 (period constraints from branch point geometry) and Shinomiya 2021 (explicit degeneration formulas in the boundary limit) fills genuine gaps in the prior citation set.

**The endpoint universality family (Part B)** is exactly the right territory for the endpoint-limit gate. Lubinsky 2020 (B1) is the most directly relevant theorem: it provides the Bessel-kernel universality limit at endpoints of the support for regular measures, which is the mathematical model for the boundary-near atom. The distinction between Airy (soft edge, density vanishing like square root at a free endpoint) and Bessel (hard edge, endpoint of the physical domain) is crucial and PC gets it right — including citing both the 2020 hard-edge result (B1) and the 2010 soft-edge result (B2) to bracket the classification question.

**The Eynard-Orantin section (Part E)** is valuable primarily for the exclusion results (E4) rather than for positive morphisms. PC correctly identifies that the smooth-ramification hypothesis makes topological recursion inapplicable at the hard-wall boundary of [-1,1]. This is a non-trivial exclusion that saves local effort — the local agent should not invest in adapting the recursion kernel for the endpoint case without a modification to handle the hard constraint.

**The Stahl-Totik section (Part C)** gives the right potential-theoretic framing. The Eichinger-Lukić 2025 entry (C4) is the most novel citation here — a 2025 paper extending Stahl-Totik regularity to continuum Schrödinger operators via Martin compactification, directly connecting the Robin constant to FIXED_CLOUD_BOUND. This is a useful and non-obvious citation.

**One citation to verify before citing downstream:** F2 (Kandaian 2026) is cited from Semantic Scholar with a URL-level identifier rather than a DOI or arXiv number. The field is actively developing and a 2026 paper may not yet have a stable DOI. Before citing this in any public artifact, the local agent should confirm the paper exists and the identifier resolves correctly.

## §2 — Morphism T4: bare a-period matrix vs. symplectic period matrix

The test for M4 (Deift et al. theta-function parametrix → normal theta-row) builds the period matrix by integrating holomorphic differentials `x^(j-1) dx / y` over support components. This constructs the a-period matrix M (real a-cycle integrals), which is the same object Round 20's sweep computed. It is *not* the same as the full symplectic period matrix Ω of the hyperelliptic surface, which includes both a-period (A-matrix) and b-period (B-matrix) entries, and satisfies Ω = A^{-1} B (in one convention).

The Deift et al. theta function `θ(n A(z) + d)` uses the Abel map A(z), which depends on a normalized basis of holomorphic differentials (normalized so that a-period integrals form the identity matrix after a change of basis). The "period row" in the theta-function context is a row of the normalized period matrix Ω = A^{-1}B, not a row of M.

This matters for T4 because:
1. The condition number of M at the boundary configuration is the Round 20 quantity (cond ~ 2.1e11 for monomial, hoping for better with Chebyshev rescaling per G2.5).
2. The condition number of the normalized Ω is a separate object that depends on the b-period integrals as well. These are harder to compute but are the mathematically correct object for testing the Deift et al. morphism.

For the toy genus-2 test that T4 specifies, the a-period matrix is 2×2 and the b-period integrals can be computed numerically. The local agent should run T4 against the full symplectic Ω rather than just M to get an answer that speaks to the Deift et al. theta-function structure. Running T4 against M alone answers a weaker question (linear independence of the a-period row) — still useful, but should be labeled correctly.

This is the same conceptual layer issue that appears in the Round 23 FT-02A diagnostic (referenced in the dispatch context). The pattern: bare a-period matrix M is what local computation produces easily; the theta-function structure depends on the normalized Ω. These come apart exactly in the high-genus, near-degenerate regime that #1038 lives in.

## §3 — Exclusion documentation: notably strong

PC's treatment of X1 and X2 is the best exclusion documentation in the #1038 assimilation record so far. Previous rounds have either omitted exclusion results or buried them in passing remarks. Here, PC gives each exclusion its own titled section, states the consequence for #1038 explicitly, and identifies what pre-check is required before the positive morphisms can be applied.

X1 (CD formula fails for biorthogonal structure) is particularly well-handled. PC correctly identifies it as a *prerequisite check* rather than a simple negative result: if the boundary component has critical density vanishing (Kuijlaars-McLaughlin higher order than square root), the entire family M1–M4 must be revised to use the Claeys-Wang framework instead. This is the right framing — it doesn't block the positive morphisms, it conditions them on a testable property of the specific cloud. T1's density-vanishing-exponent fit is designed exactly to run this prerequisite check.

X2 (Eynard-Orantin smooth-ramification hypothesis) closes off a potential time sink. Without this exclusion, a local agent might invest effort in adapting the topological recursion for the endpoint-limit case only to find the theorem doesn't apply. PC saves that effort.

## §4 — HTTP 404 reading NEXT_LOCAL_GATE.md: worth investigating

PC's MANIFEST.json notes: "HTTP 404 at EXTERNAL-REVIEW-ASSIMILATION-ROUND_19/NEXT_LOCAL_GATE.md — path may have moved."

The file exists locally at `Erdos1038/agent-work/problem-at-hand/EXTERNAL-REVIEW-ASSIMILATION-ROUND_19/NEXT_LOCAL_GATE.md` (verified by `ls`). The Round 20 MANIFEST.json confirms PC successfully read the Round 19 branch at `94606179a6deb9c60fbdd6d44498a7470d67da3f` during that round, so PC's sandbox does have file-access capability on the branch.

Three plausible explanations:

1. **Branch head difference**: Round 20 PC read branch head `94606179`; Round 22 PC read head `33109f0`. The Round 19 assimilation directory was committed to the branch at some point; if the commit history between these two heads moved the directory or changed the path casing, the 404 would occur.

2. **Path case-folding**: macOS is case-insensitive by default; PC's sandbox may run on a case-sensitive Linux filesystem. If any path component in `EXTERNAL-REVIEW-ASSIMILATION-ROUND_19/` differs by case between the filesystem and the reference, the 404 would occur on Linux but not locally.

3. **PC accessing origin/main instead of the feature branch**: PC confirmed branch checkout, but if there was a fallback in the file-read mechanism (not the git checkout mechanism), the read could have been against main, which may not have the assimilation directory. PC correctly identified this risk in its manifest note.

**Assessment**: this is not a blocking issue. PC correctly pivoted to the playback-row tail for the information it needed (Task C deferred status). The investigation is worth doing before Round 23 dispatch if any file-reads are expected from older assimilation directories — route the specific file path to a fresh `git show` or `git ls-files` check against the branch head to confirm visibility from the remote.

## §5 — Score-card

| confound | Status entering Round 22 | Round 22 outcome |
|---|---|---|
| C1 (tautological QR diagnostic) | RESOLVED (Round 20) | Carry-forward: unchanged |
| C2 (unprobed g=24) | RESOLVED (Round 20) | Carry-forward: unchanged |
| C3 (f64-only sampling) | STILL_OPEN | Carry-forward: unchanged |
| C4 (silent main-fallback) | RESOLVED_BY_TEMPLATE (Round 20) | Carry-forward: no fallback in R22 |
| C5 (boundary component rank-drop) | Not yet identified | **NEW — OPEN**: surfaced by Bogatyrev M3 |

Round 22 adds one new confound (C5) and provides literature tools to address it (T3, and the Round 21 G2.5 cond(M_T) sweep that is already in flight). The confound is real but falsifiable locally. If G2.5 shows cond(M_T) well below 1e10 for the specific #1038 cloud's boundary component width, C5 resolves. If cond(M_T) is still high, C5 becomes the next methodology gap to close.

No altitude movement. No claim promotion. Six receipts still absent.
