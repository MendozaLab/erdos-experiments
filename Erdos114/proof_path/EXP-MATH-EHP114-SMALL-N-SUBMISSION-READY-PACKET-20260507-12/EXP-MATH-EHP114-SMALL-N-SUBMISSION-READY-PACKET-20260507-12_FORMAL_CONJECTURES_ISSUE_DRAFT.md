# Formal-conjectures issue draft (alternative to PR)

> Use this only if upstream maintainers prefer a discussion issue before a PR.

## Title

`Discussion: finite variant of Erdos #114 (n <= 14 except n=13)?`

## Body

I have a finite-degree variant of Erdős–Herzog–Piranian (Problem 114) that I'd like to contribute as a separate theorem statement, without changing the status of the open `erdos_114`.

The finite statement would cover monic complex polynomials of degree `n` with `1 <= n <= 14` and `n != 13`. The dependencies are typed: `n = 1` is the translated unit circle, `n = 2` is MacLane (1953/54) with Eremenko–Hayman (1999) as the accessible source, and `3 <= n <= 12, n = 14` are reproducible Rust + `inari` IEEE-1788 interval certificates with row-level SHA-256 sidecars archived under DOI [10.5281/zenodo.19480329](https://zenodo.org/records/19480329).

I'm proposing to register the per-degree dependencies as explicit axioms — one per certified degree — rather than tactic-closing them, so reviewers can audit which artifacts each axiom is relying on. A draft Lean 4 file is in the accompanying packet.

The `n = 13` row has a route exception (the underlying experiment artifact reports zero branch-and-bound evaluations against tens-to-hundreds of millions in neighboring rows) and is intentionally excluded from the finite statement until the row is re-run.

The all-degree conjecture and Tao's sufficiently-large-`n` theorem (arXiv:2512.12455) are deliberately out of scope.

Question: would you prefer this as a PR adding `Erdos114.erdos_114_finite_le_14_except_13` alongside the open `erdos_114`, or as a separate file under a `Finite/` subdirectory? And do you want each per-degree certificate as a named axiom (current draft) or as a single hypothesis with row-level data?

## Tooling disclosure

Draft prepared with AI-assisted tooling. The mathematical content rests on the cited literature and SHA-checked certificate artifacts only; no AI output is used as proof evidence.
