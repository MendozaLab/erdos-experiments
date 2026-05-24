# Formal-conjectures PR packet for Erdős Problem 114 (finite variant)

## Intended PR title

`Erdos114: add finite n <= 14 certificate variant; do not change main status`

## Intended status change

**Do not** change the status of `Erdos114.erdos_114`. It stays open. The PR adds a sibling theorem `erdos_114_finite_le_14` to the same file or namespace, plus its 14 per-degree dependencies, with the explicit understanding that the finite variant is not a partial proof of the open conjecture but a separately scoped statement.

## What we have today

The existing `FormalConjectures/ErdosProblems/114.lean` already cites our v5 Zenodo DOI 10.5281/zenodo.19480329 (per DeepMind PR #3712). The main `erdos_114` theorem in that file remains an open `sorry`. This PR does not modify or amend that theorem or that PR — it adds a sibling theorem that covers the explicit finite block.

## Proposed file: `FormalConjectures/ErdosProblems/114.lean` (additions to existing file) or `FormalConjectures/ErdosProblems/Finite/114.lean` (new file in same namespace)

A complete Lean 4 stub is included alongside this packet at `EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-13_LEAN_STUB.lean`. The stub:

- Keeps `Erdos114.erdos_114` exactly as upstream defines it (untouched).
- Adds `Erdos114.erdos_114_finite_le_14` as a separate theorem statement.
- Records the per-degree dependencies as explicit `axiom` declarations rather than hiding them in a tactic block.
- Tags the literature row (`n = 2`) as an axiom that names MacLane and Eremenko–Hayman in the docstring.
- Tags every certified row (`n = 3..14`) as an axiom that points to the SHA-256 sidecar of the corresponding result JSON.

The intent is to make every external dependency a named axiom rather than a tactic-closed lemma, so reviewers can audit what is being assumed without reading the body of any proof term.

## PR body (proposed)

This PR adds a finite variant of Erdős–Herzog–Piranian (Problem 114) covering every monic complex polynomial of degree `n` with `1 ≤ n ≤ 14`. It does not change the status of the main conjecture `erdos_114`, which remains open.

The finite variant is dependency-typed:

- `n = 1`: direct analytic — a monic linear polynomial has a translated unit circle as its unit lemniscate.
- `n = 2`: literature row, citing G. R. MacLane, "On a conjecture of Erdős, Herzog, and Piranian," Michigan Math. J. 2 (1953/54), 147-148, doi:10.1307/mmj/1028989918, with Alexandre Eremenko and Walter K. Hayman, "On the length of lemniscates," Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295 as the accessible source pin. The Eremenko–Hayman abstract names Bernoulli's lemniscate as the `d = 2` extremal level set; the proof remark after Lemma 5 identifies `z^2 + 1` as extremal, length-equivalent to `z^2 - 1` by rotation.
- `3 ≤ n ≤ 14`: each degree is recorded as an axiom referring to a DOI-archived Rust + `inari` IEEE-1788 interval certificate, with the result JSON SHA-256 sidecar named in the axiom docstring.

The submission is intentionally not a partial step toward the all-degree conjecture; it is a finite frontier whose boundary is set by the current certified-degree limit. Tao's "Lemniscate length and the Erdős–Herzog–Piranian conjecture" (arXiv:2512.12455) covers all sufficiently large `n`. The remaining bridge question — whether Tao's threshold can be made effective at or below the present finite ceiling — is not addressed in this PR.

The certificate set is archived at concept DOI [10.5281/zenodo.19184467](https://doi.org/10.5281/zenodo.19184467) (resolves to latest); source repository at <https://github.com/MendozaLab/erdos-experiments>, release [`v3.1.0`](https://github.com/MendozaLab/erdos-experiments/releases/tag/v3.1.0).

## Suggested docstring for `erdos_114_finite_le_14`

```
/--
For every monic complex polynomial `p` of degree `n` with `1 ≤ n ≤ 14`,
the lemniscate length of `{|p| = 1}` is at most the lemniscate length
of `z^n - 1`.

Dependencies:
- `n = 1` is the translated unit circle.
- `n = 2` is MacLane (Michigan Math. J. 2 (1953/54), 147-148,
  doi:10.1307/mmj/1028989918), with Eremenko-Hayman (Michigan Math. J. 46
  (1999), 409-415; arXiv:0805.2295) as accessible pin.
- `3 ≤ n ≤ 14` are recorded as axioms backed by Rust + inari
  IEEE-1788 interval certificates; SHA-256 sidecars named per axiom.
- `n = 13` is the row that was re-run on 2026-05-07 (197M B&B
  evaluations, 53min wall-clock on Modal 32-CPU) after a verdict-logic
  bug was patched in commit dae62b8.

This statement is separate from the open all-degree conjecture
`erdos_114` and from Tao's sufficiently-large-`n` theorem
(arXiv:2512.12455). No threshold is extracted from Tao's theorem here,
and no claim is made for any `n ≥ 15`.

Certificate set archived at concept DOI 10.5281/zenodo.19184467.
-/
```

## Lean 4 stub file

See `EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-13_LEAN_STUB.lean` in this packet.

## Review questions (one)

Is this the right shape — keep `erdos_114` open, add a separately named finite theorem `erdos_114_finite_le_14` with each per-degree dependency as an explicit axiom, contiguous coverage `1 ≤ n ≤ 14`?

## Claim ceiling

This PR does not assert the all-degree Erdős #114 conjecture, does not authorize any `n ≥ 15` computation, and does not extract a threshold from Tao (arXiv:2512.12455).

## AI assistance disclosure (Block H form)

Lean 4 code in this PR was drafted with assistance from Claude (Anthropic, claude-opus-4-7, via the Claude Code CLI). Citation verification used Perplexity (web app, quorum mode including claude-opus-4-7, gemini-3-pro-deep-think, gpt-5-pro). I selected the theorem statements, reviewed every proof body, and verified the axiom citations against published sources (MacLane 1953/54, Eremenko–Hayman 1999) and the DOI-archived interval certificates. All non-axiom proofs are machine-checked by Lean 4 against the toolchain pinned in `lean-toolchain` and `lake-manifest.json`.

## Tooling disclosure

This packet and the proposed PR body were prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature and SHA-checked certificate artifacts; no AI output is used as proof evidence. The author takes responsibility for the statement, axiomatization choices, and PR wording.
