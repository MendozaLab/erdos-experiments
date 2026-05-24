# Zenodo strategy for the v12 submission packet

## Decision: edit the v5 description, do not mint v6

The v5 record (DOI 10.5281/zenodo.19480329) holds the n=3..14 certificate JSONs and SHA-256 sidecars. The v12 packet does not change those artifacts; it changes presentation, attribution wording, and the `n = 13` route disclosure.

A new Zenodo version is appropriate when the underlying artifacts change. A description edit is appropriate when only narrative changes. Mixing those signals confuses what the DOI is certifying.

The v5-pinned citation in DeepMind formal-conjectures `FormalConjectures/ErdosProblems/114.lean` continues to resolve to the same artifact set after a description edit; minting v6 would not break that citation but would orphan it slightly behind the latest concept-DOI tip.

## Description edit (proposed text for v5)

Replace the current description with:

```
Reproducible IEEE-1788 interval certificates for the lemniscate length
of monic complex polynomials of degree 3 <= n <= 14, computed in Rust
using the `inari` IEEE-1788 interval-arithmetic crate, in the context of
the Erdős–Herzog–Piranian conjecture (Erdős Problem 114).

Each n is certified by a separate row in
results/erdos-114/EXP-MM-EHP-007-n{N}-inari_RESULTS.json with a
SHA-256 sidecar at the same path with the .sha256 extension.

Companion small-degree variant (1 <= n <= 14, excluding n = 13) records
the dependency split:

  - n = 1: direct analytic, translated unit circle.
  - n = 2: literature row, citing G. R. MacLane, "On a conjecture of
    Erdős, Herzog, and Piranian," Michigan Math. J. 2 (1953/54),
    147-148, doi:10.1307/mmj/1028989918, and Alexandre Eremenko and
    Walter K. Hayman, "On the length of lemniscates," Michigan Math.
    J. 46 (1999), 409-415; arXiv:0805.2295.
  - 3 <= n <= 12, n = 14: the certificate JSONs in this record.

Note on n = 13: the row EXP-MM-EHP-007-n13-inari_RESULTS.json reports
bb_total_evals = 0 and bb_level_count = 0, while neighboring rows
report 4.5×10^7 (n=12) and 8.5×10^8 (n=14). The row is preserved in
this record for transparency and is excluded from the small-degree
variant theorem until the row is re-run with the same evaluator
profile as n != 13.

The all-degree Erdős #114 conjecture is not asserted by this record,
and no threshold is extracted from Tao's sufficiently-large-n theorem
(arXiv:2512.12455).

Source repository: https://github.com/MendozaLab/erdos-experiments
Erdős Problem #114 reference: https://www.erdosproblems.com/114
```

## When v6 would be appropriate

Mint a new Zenodo version (v6) only if any of the following hold:

- The `n = 13` row is re-run with a populated branch-and-bound trace and the new JSON+SHA replaces the existing v5 row.
- A new degree (`n = 15` or higher) is added.
- The certificate format changes (e.g., new `inari` version, new schema fields).
- A row's `l_star_lower` or `l_star_upper` changes.

Editorial-only changes (axiomatization choices, formal-conjectures wording, presentation packets) do not justify a new version.

## Concept DOI vs versioned DOI

- Concept DOI (always latest): [10.5281/zenodo.19184467](https://doi.org/10.5281/zenodo.19184467)
- v5 DOI (current): [10.5281/zenodo.19480329](https://doi.org/10.5281/zenodo.19480329)

External citations should target the concept DOI when long-term traceability is wanted and the versioned DOI when the audit must pin to a specific certificate set. Past practice in the formal-conjectures PR has used the versioned DOI; that's appropriate when the citing artifact is itself version-pinned.

## Process gate

The Zenodo description edit must wait until prepub-redteam, Publisher Gate, and Crackpot-Scrub have all run on the edit text. The Zenodo edit endpoint is treated as a public publication surface for purposes of the Publisher Gate.
