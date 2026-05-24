# Formal-conjectures issue draft (alternative to PR)

> Use this only if upstream maintainers prefer a discussion issue before a PR.

## Title

`Discussion: finite variant of Erdos #114 (1 ≤ n ≤ 14)?`

## Body

I have a finite-degree variant of Erdős–Herzog–Piranian (Problem 114) that I'd like to contribute as a separate theorem statement, without changing the status of the open `erdos_114`.

The finite statement covers monic complex polynomials of degree `n` with `1 ≤ n ≤ 14` contiguously. The dependencies are typed: `n = 1` is the translated unit circle, `n = 2` is MacLane (1953/54) with Eremenko–Hayman (1999) as the accessible source, and `3 ≤ n ≤ 14` are reproducible Rust + `inari` IEEE-1788 interval certificates with row-level SHA-256 sidecars archived under concept DOI [10.5281/zenodo.19184467](https://doi.org/10.5281/zenodo.19184467).

I'm proposing to register the per-degree dependencies as explicit axioms — one per certified degree — rather than tactic-closing them, so reviewers can audit which artifacts each axiom is relying on. A draft Lean 4 file is in the accompanying packet.

The `n = 13` row was re-run on 2026-05-07 after a verdict-logic bug at high reduced dimension was patched in commit `dae62b8` of MendozaLab/erdos-experiments. The clean re-run completed with `bb_total_evals = 197,132,288` in 52.96 minutes on a 32-CPU Modal worker; released as `v3.1.0`.

The all-degree conjecture and Tao's sufficiently-large-`n` theorem (arXiv:2512.12455) are deliberately out of scope.

Question: would you prefer this as a PR adding `Erdos114.erdos_114_finite_le_14` alongside the open `erdos_114` in the existing `FormalConjectures/ErdosProblems/114.lean`, or as a separate file under `FormalConjectures/ErdosProblems/Finite/114.lean` in the same namespace? And do you want each per-degree certificate as a named axiom (current draft) or as a single hypothesis with row-level data?

## Tooling disclosure

Draft prepared with AI-assisted tooling. The mathematical content rests on the cited literature and SHA-checked certificate artifacts only; no AI output is used as proof evidence.
