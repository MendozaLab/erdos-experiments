# EHP114 v13 submission packet — narrative report

## Why this version exists

The v12 packet was built around a quarantine of the n=13 certificate row. Between v12 freezing on 2026-05-07 and v13 being assembled later the same day, three things happened: (1) the verdict-logic bug that produced the v12 n=13 anomaly was traced to `ehp_general_ieee1788.rs:987` and patched in commit `dae62b8` of MendozaLab/erdos-experiments; (2) the patched binary was used to re-run n=13 cleanly (commit `f89597e`, 197M B&B evaluations, 53 minutes wall-clock on a 32-CPU Modal worker); (3) a `v3.1.0` GitHub release shipped the corrected certificate, with a Block-G AI-acknowledgment in the release notes. v13 incorporates that re-run across all packet artifacts and updates the finite theorem statement from `1 ≤ n ≤ 14, n ≠ 13` to the contiguous `1 ≤ n ≤ 14`.

## What the finite theorem now says

For every monic complex polynomial `p` of degree `n` with `1 ≤ n ≤ 14`, the lemniscate length `L({z : |p(z)| = 1})` is at most the lemniscate length of `z^n - 1`. There is no degree exclusion.

Three kinds of evidence stand behind that statement, exactly as in v12 with one column refilled:
- `n = 1` is a direct analytic calculation
- `n = 2` is a literature citation, MacLane (1953/54) historical, Eremenko–Hayman (1999) accessible source pin
- `3 ≤ n ≤ 14` are reproducible IEEE-1788 interval certificates produced by a Rust + `inari` pipeline, archived under concept DOI 10.5281/zenodo.19184467 with row-level SHA-256 sidecars. The n=13 row in this set is the 2026-05-07 re-run.

## What the finite theorem still does not say

This statement is not a partial proof of the all-degree Erdős #114 conjecture. It is not in scope of Tao's sufficiently-large-`n` theorem (arXiv:2512.12455). It does not extract a threshold from Tao. It does not authorize any computation for `n ≥ 15`. It is a finite frontier — the small side of a small-and-large pincer that has not been closed.

## Why v12 stays on disk

v12 was FROZEN-locked under its 2026-05-07-r01 redteam round before the n=13 re-run completed. Per the H² portfolio's never-delete rule, v12 is preserved as the historical record of "what we would have shipped if n=13 hadn't been re-run." The FROZEN.lock is intact; the audit folder is unchanged. v13 is not an edit of v12; it is a sibling packet that supersedes v12 forward.

## Why each per-degree certificate stays a Lean axiom

The formal-conjectures style (per H² Formalization Integrity Protocol and DeepMind upstream) is to make external dependencies explicit. A tactic-closed lemma whose statement requires an external interval-arithmetic certificate is a hidden sorry. Naming each row as an axiom — one per degree, with the SHA-256 sidecar in the docstring — keeps the audit surface explicit. v13's Lean stub now contains 12 axioms (n=3 through n=14) plus the two literature/analytic axioms for n=1 and n=2, with the n=13 axiom carrying a provenance note about the verdict-bug-fix-and-re-run history.

## Zenodo record update

A certificate row changed (n=13). v12's Zenodo strategy memo explicitly listed this exact case as the trigger for a new version. The GitHub→Zenodo auto-archive did not appear publicly after re-check, so on 2026-05-08 the corrected v3.1.0 packet was manually published as version DOI `10.5281/zenodo.20087919` under the existing concept DOI family `10.5281/zenodo.19184467`.

## Remaining process gates before forum/social posting

In order:

1. **prepub-redteam** on the v13 frozen post — `init … --new-round` (since v12 round froze a stale post), re-import the Perplexity quorum critique (still applicable; the L-definition / deg-p-exact / Tao-citation findings are content-level and not n=13-specific), log edits, freeze, verify.
2. `anthropic-skills:publisher` over the v13 frozen post and the assembled paper at `Math/preprints/ehp114-finite/`.
3. `crackpot-scrub` over the same artifacts.
4. Post to erdosproblems.com using the version DOI `10.5281/zenodo.20087919` and the open-conjecture caveat.
5. Use any social pointer posts only as short links to the Zenodo/GitHub/forum record.

The certificate packet is now deposited on Zenodo. The local v13 editorial packet remains the source for the forum wording and downstream reviewer-facing text.

## Tooling disclosure

This packet, including the post draft, the formal-conjectures packet, the Lean 4 stub, the Zenodo strategy memo, and this report, was prepared with AI-assisted tooling. The mathematical claims rest on cited literature and SHA-checked certificate artifacts only. No AI output is used as proof evidence. The author takes responsibility for the mathematical claims, axiomatization choices, and submission wording.
