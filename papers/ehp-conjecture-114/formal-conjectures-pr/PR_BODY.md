## Summary

Adds `FormalConjectures/ErdosProblems/114.lean` — the Erdős–Herzog–Piranian (EHP) conjecture (#114) with a computationally certified small-*n* result.

- **`erdos_114`**: Full open conjecture (all *n*)
- **`erdos_114_small_n`**: Solved for 3 ≤ *n* ≤ 14 via IEEE 1788-rigorous interval arithmetic

## Context

Moritz Firsching confirmed a standalone PR is welcome:
[Zulip message](https://leanprover.zulipchat.com/#narrow/channel/524981-Formal-conjectures/topic/How.20to.20proceed.20on.20.22All.201179.20Erd.C5.91s.20conjectures.20formalized.22/near/584182899)

This PR is intentionally separate from #3422.

## Certificates

Each *n* case is verified independently by branch-and-bound over the compact
parameter space of monic degree-*n* polynomials, with IEEE 1788-2015 certified
interval arithmetic bounds:

- *n* = 3–12: Python/mpmath
- *n* = 13–14: Rust/inari (MPFR-backed)

All results deposited with SHA-256 checksums:
**[doi:10.5281/zenodo.19322367](https://doi.org/10.5281/zenodo.19322367)**

## Definitions

| Name | Description |
|------|-------------|
| `levelCurveUnit p` | Lemniscate {z ∈ ℂ : ‖p(z)‖ = 1} |
| `arcLength p` | 1-dimensional Hausdorff measure of the level curve |

## Future work

- Uniqueness theorem (`erdos_114_small_n_unique`) as a follow-up
- Extension to *n* = 15+ as compute budget allows
