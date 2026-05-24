# Edit Log — 2026-05-07-r01

Append-only log of edits applied in response to external critique.
Each row: timestamp, description, optional diff path.

---


## 2026-05-07T13:07:26 — Revise post in response to Perplexity quorum critique (Claude Opus 4.7 + Gemini 3.1 Pro + GPT-5.5) and the verdict-bug remediation completed today. Specific changes: (1) define lemniscate length explicitly as 1D Hausdorff measure (Claude Opus 4.7 + Gemini high-confidence); (2) state degree as 'deg p = n exactly' rather than the implicit 'monic of degree n' (Claude Opus 4.7 unique); (3) list certified degrees as the explicit set {1,2,3,...,12,14} (Claude Opus 4.7 unique); (4) add explicit Tao arXiv ID 2512.12455 inline (Claude Opus 4.7 + Gemini high-confidence); (5) replace 'SHA-256 sidecars' with 'artifact-integrity hashes (SHA-256)' (Gemini crypto-bombast finding); (6) rewrite n=13 paragraph: now traced to a verdict-logic bug at high reduced dimension, patched at source today (commit dae62b8 in MendozaLab/erdos-experiments); (7) replace ad-hoc tooling sentence with Block B SHORT FORM AI acknowledgment per AI_ACKNOWLEDGMENT_TEMPLATES.md, naming claude-opus-4-7 and the Perplexity quorum models explicitly. Empty Lean code fence in the redteam prompt template is a separate item (post itself never had Lean source).


## 2026-05-07T13:12:07 — Canonical-phrasing patch on AI-disclosure first sentence per AI_ACKNOWLEDGMENT_TEMPLATES.md Block B binding template. Replaced 'research direction, problem framing, and error-catching' with the canonical exact-match form 'research direction, problem identification, alternative framings, and error-catching'. No semantic change to the disclosure; this aligns the post with the binding template's first-sentence requirement (which is the only mandatory non-droppable element of Block B). Triggered by Publisher gate finding (4) earlier this turn.

