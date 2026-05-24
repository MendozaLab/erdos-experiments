# Edit Log — 2026-05-07-r01

Append-only log of edits applied in response to external critique.
Each row: timestamp, description, optional diff path.

---


## 2026-05-07T21:19:08 — v12 -> v13 evolution: incorporated the n=13 re-run that landed in commit f89597e and was released as v3.1.0. Specific edits in the post draft: (1) theorem statement updates from '1 <= n <= 12, n = 14' to contiguous '1 <= n <= 14'; (2) n=13 paragraph rewrites from quarantine narrative to verdict-bug-fix-and-rerun narrative (cause traced to ehp_general_ieee1788.rs:987, patched in dae62b8, re-run 2026-05-07 with 197M B&B evals over 53min on Modal 32-CPU); (3) certificate set citation updates to concept DOI 10.5281/zenodo.19184467 (resolves to latest, expected v6 from auto-archive of v3.1.0 release); (4) Block B AI-acknowledgment carries through unchanged — first sentence still canonical, still names claude-opus-4-7 and the Perplexity quorum models. The Perplexity quorum critique (claude-opus-4-7 + gemini-3-pro-deep-think + gpt-5-pro) imported from v12's redteam round; its findings remain applicable since they were content-level (L definition, deg p exact, Tao arXiv ID, Lean fence) and were already addressed in v12 with carry-through to v13.


## 2026-05-07T21:20:28 — Quorum-critique carry-forward verification for v13. Confirmed each Perplexity quorum (Claude Opus 4.7 + Gemini 3.1 Pro + GPT-5.5) high-confidence finding from v12 is preserved in the v13 post text: (a) L is defined as 1D Hausdorff measure inline ('the one-dimensional Hausdorff measure of {z : |p(z)| = 1}'); (b) deg p = n exactly is stated explicitly; (c) certified degrees listed as the contiguous range '1 <= n <= 14' (no exclusion in v13); (d) Tao arXiv:2512.12455 cited inline; (e) 'artifact-integrity hashes (SHA-256)' phrasing replaces the prior crypto-bombast wording. Block B AI-acknowledgment first sentence is canonical-exact ('research direction, problem identification, alternative framings, and error-catching'). v13 carries no fresh review findings beyond the v12 quorum (no new external review run for v13 since the content evolution is purely incorporating the n=13 re-run that landed post-quorum).

