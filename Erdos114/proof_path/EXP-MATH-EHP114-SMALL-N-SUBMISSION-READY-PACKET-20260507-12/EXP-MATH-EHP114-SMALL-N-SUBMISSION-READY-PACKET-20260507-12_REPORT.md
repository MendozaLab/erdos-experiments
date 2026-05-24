# EHP114 v12 submission packet — narrative report

## Why this version exists

The v11 packet had three reviewer-facing problems flagged in pre-submission review on 2026-05-07: an unannotated `route exception flagged` annotation on the `n = 13` row, a `LOCAL_INTERVAL_CERTIFIED` status string that reads as a hedge to outside readers, and an erdosproblems.com post draft formatted as a filing packet rather than a problem-page note. v12 fixes each of those and adds a Lean 4 stub plus a Zenodo strategy memo.

## What the finite theorem actually says

For every monic complex polynomial `p` of degree `n` with `1 <= n <= 14` and `n != 13`, the lemniscate length `L({z : |p(z)| = 1})` is at most the lemniscate length of `z^n - 1`. The exclusion of `n = 13` is documented and quarantined; see "n=13 route note" in the README.

Three kinds of evidence stand behind that statement. `n = 1` is direct: a monic linear polynomial gives a translated unit circle of length `2π`. `n = 2` is a literature citation, with MacLane (1953/54) as the historical row and Eremenko–Hayman (Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295) as the accessible source pin. `n = 3..12` and `n = 14` are reproducible IEEE-1788 interval certificates produced by a Rust + `inari` pipeline, archived under DOI 10.5281/zenodo.19480329 with row-level SHA-256 sidecars.

## What the finite theorem does not say

This statement is not a partial proof of the all-degree Erdős #114 conjecture. It is not in scope of Tao's sufficiently-large-`n` theorem (arXiv:2512.12455). It does not extract a threshold from Tao. It does not authorize any computation for `n >= 15`. It is a finite frontier — the small side of a small-and-large pincer that has not been closed.

## Why `n = 13` is excluded

The certificate JSON `EXP-MM-EHP-007-n13-inari_RESULTS.json` reports `bb_total_evals = 0` and `bb_level_count = 0`, while neighboring rows (`n = 12`: 4.5×10⁷, `n = 14`: 8.5×10⁸) are populated. The interval bounds `L*_lower = 28.85923995588822`, `L*_upper = 28.859239955888235` are present in the row, but the branch-and-bound search structure is not. Two interpretations are consistent with the artifact: a different evaluator pass was used, or the row is an artifact-of-tooling. Either way, shipping a non-uniform table to a public table risks "what is this row?" as the first reviewer comment. Quarantining `n = 13` is the cheaper move.

The honest options are (a) re-run `EXP-MM-EHP-007-n13-inari` with the same evaluator profile as `n != 13`, then unblock; or (b) leave `n = 13` quarantined and ship `n != 13`. v12 is built around option (b).

## Why the Zenodo record is being edited rather than versioned

The certificate JSONs that v12 cites are exactly the ones on v5 of the Zenodo record. No artifact is changing. What is changing is the description text and the dependency disclosure. A description edit fits that; a v6 mint would imply new evidence and break the existing PR-#3712 citation pattern.

`ZENODO_STRATEGY.md` records the proposed v5 description text and the conditions under which v6 would be appropriate (n=13 re-run, new degrees, schema change, or any row's `l_star_*` changing).

## Why each per-degree certificate is being recorded as a Lean axiom

The formal-conjectures style (per H² Formalization Integrity Protocol and DeepMind upstream) is to make external dependencies explicit. A tactic-closed lemma whose statement requires an external interval-arithmetic certificate is a hidden sorry. Naming each row as an axiom — one per degree, with the SHA-256 sidecar in the docstring — keeps the audit surface explicit. Reviewers can see in one pass exactly what the finite theorem is assuming.

The `n = 13` axiom is intentionally absent from the Lean stub. The finite theorem statement excludes `n = 13` rather than asserting an axiom and never using it.

## Process gates that must run before any public deposit

In order:

1. Re-run `EXP-MM-EHP-007-n13-inari` (only if shipping the n=13 row). v12 is built to ship without this; if the row is re-run, update v12 → v13 with `n = 13` re-included.
2. `prepub-redteam init → prompt → import-critique → log-edit → freeze → verify` over the post draft, the formal-conjectures PR body, and the Zenodo description edit. Produces `FROZEN.lock`; the deposit step must refuse without it.
3. `anthropic-skills:publisher` over the same three artifacts.
4. `crackpot-scrub` over the same three artifacts.
5. Edit Zenodo v5 description.
6. Post to erdosproblems.com.
7. Open the formal-conjectures PR (or issue, per upstream preference) with the Lean 4 stub.

The packet itself is local and editorial; nothing in v12 has been deposited.

## Tooling disclosure

This packet, including the post draft, the formal-conjectures packet, the Lean 4 stub, the Zenodo strategy memo, and this report, was prepared with AI-assisted tooling. The mathematical claims rest on cited literature and SHA-checked certificate artifacts only. No AI output is used as proof evidence. The author takes responsibility for the mathematical claims, axiomatization choices, and submission wording.
