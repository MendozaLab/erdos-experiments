# Formal-conjectures packet for Erdos Problem 114

## Intended PR title

`Erdos114: add finite n < 15 certificate theorem`

## Intended status change

Do **not** mark the main theorem `erdos_114` as closed. Keep the main statement open unless and until the Tao threshold bridge is explicit.

## Proposed formal-conjectures shape

- Keep namespace: `Erdos114`.
- Keep the main theorem as `@[category research open, AMS 30] theorem erdos_114 ... := by sorry`.
- Add or retain the finite variant `erdos_114_finite_lt_15`.
- Mark the finite variant with the repository's resolved-research category only if the certificate axioms are accepted as explicit external dependencies.
- Keep every machine-certified row as an explicit axiom or external-certificate lemma; do not hide computational certification inside a proof term.

## Suggested docstring wording

"Finite certified range for Erdos Problem 114. For `1 <= n < 15`, the Erdos-Herzog-Piranian lemniscate inequality holds, using the direct `n=1` case, the Eremenko-Hayman `n=2` case, and DOI-backed interval certificates for `3 <= n <= 14`. This finite statement is separate from the open all-degree conjecture and from Tao's sufficiently-large-`n` theorem."

## PR body

This PR records a finite variant of Erdos Problem 114 rather than changing the status of the main conjecture. The finite theorem covers exactly `1 <= n < 15`.
The proof dependencies are typed explicitly:

- direct analytic lemma for `n=1`;
- literature-only Eremenko-Hayman input for `n=2`;
- DOI-backed interval certificates for `3 <= n <= 14`;
- Tao's sufficiently-large-`n` theorem is mentioned only as context for the remaining bridge question.

The public certificate record is [10.5281/zenodo.19480329](https://zenodo.org/records/19480329); the code/results repository is [https://github.com/MendozaLab/erdos-experiments](https://github.com/MendozaLab/erdos-experiments).  Local preflight checks report missing degree rows `0`, source SHA failures `0`, and `full_ehp114_claim = false`.

## Review question

Is this the right formal-conjectures shape: keep `erdos_114` open, add the resolved finite variant, and leave the Tao-threshold bridge outside this PR until an explicit threshold is available?

## Claim ceiling

Submission preflight packet only. It prepares erdosproblems.com and Google DeepMind formal-conjectures wording for the finite n<15 result and the conditional Tao-threshold bridge; it does not assert the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.
