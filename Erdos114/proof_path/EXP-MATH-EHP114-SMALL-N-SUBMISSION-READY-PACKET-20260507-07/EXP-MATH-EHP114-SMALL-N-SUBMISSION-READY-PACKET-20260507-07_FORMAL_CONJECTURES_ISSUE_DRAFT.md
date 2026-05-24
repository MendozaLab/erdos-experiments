# Draft GitHub Issue: Formal-conjectures shape for Erdos 114 finite variant

## Title

Question: preferred shape for an Erdos 114 finite `n < 15` variant?

## Body

I am preparing a small, finite-range entry for Erdős Problem 114 and would like maintainer guidance before opening a PR.

The intended theorem is not the all-degree conjecture. It is the finite variant:

`For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the length for p(z)=z^n-1.`

Proposed shape:

- keep the main `erdos_114` theorem open;
- add a named finite theorem, e.g. `erdos_114_finite_lt_15`;
- keep the dependencies explicit as external axioms/certificate lemmas rather than hiding the computational part in proof terms;
- use the direct `n=1` row, the MacLane / Eremenko-Hayman `n=2` literature row, and DOI-backed Rust/inari IEEE-1788 certificate rows for `3 <= n <= 14`.

The public certificate record for `3 <= n <= 14` is [10.5281/zenodo.19480329](https://zenodo.org/records/19480329). Row-level result JSON files and SHA-256 sidecars are available in [https://github.com/MendozaLab/erdos-experiments](https://github.com/MendozaLab/erdos-experiments). The `n=2` literature row is historically MacLane, with accessible source pin Eremenko-Hayman, *On the length of lemniscates*, Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295.

Important non-claim: No explicit threshold is extracted from Tao's sufficiently-large-n theorem in this finite packet, and no degree n >= 15 is claimed here.

Question for maintainers: would this be welcome as a finite variant alongside the open all-degree statement, or would you prefer the computational certificate dependencies to be represented differently?

Local preflight status: missing degree rows `0`, dependency SHA failures `0`, full all-degree claim `false`.
