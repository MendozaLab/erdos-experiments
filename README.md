# Erdős Problem Computational Experiments

Computational solvers and experimental results for open Erdős problems.

**Author**: Kenneth A. Mendoza ([ORCID 0009-0000-9475-5938](https://orcid.org/0009-0000-9475-5938))

---

## Problem #114 — Erdős–Herzog–Piranian Lemniscate Conjecture

**Statement**: Among all monic polynomials of degree n, the lemniscate {z ∈ ℂ : |p(z)| = 1} has maximal arc length when p(z) = zⁿ − 1.

**Prize**: $250 · **Status**: Open for small n (solved for sufficiently large n by Tao, Dec 2025)

### Background

Tao (arXiv:2512.12455) proved EHP for all sufficiently large n. The small-n regime — the original computational challenge — remains open. Eremenko–Hayman settled n = 2 analytically in 1999. This repository addresses n = 3 through n = 8.

### Results

| n | Verdict | Method | Margin over nearest competitor | Samples |
|---|---------|--------|-------------------------------|---------|
| 3 | **VERIFIED** | Branch-and-bound proof, symmetry-reduced | 15.2% | 200,000 |
| 4 | Supported | B&B eliminated 99.6% of search space | 20.6% | 100,000 |
| 5 | Supported | Feasibility mode | 39.3% | 50,000 |
| 6 | Supported | Feasibility mode | 52.6% | 20,000 |
| 7 | Supported | Feasibility mode | 61.7% | 10,000 |
| 8 | Supported | Feasibility mode | 61.7% | 5,000 |

**Zero counterexamples across 435,000 total samples for n = 3–8.**

The n = 3 result is the first computationally verified case for any specific n > 2.

**Additional finding**: The circular symmetrization hypothesis (P114-3) was computationally falsified — symmetrization does not preserve lemniscate length. Tao's proof does not use symmetrization, independently confirming this.

### Exact Formula (verified)

For the conjectured maximizer:

    L(zⁿ − 1) = 2^{1/n} · √π · Γ(1/(2n)) / Γ(1/(2n) + 1/2)

Verified against Richardson extrapolation (5 resolutions: 200–3200) for n = 3–8.

### Methodology

- **Symmetry reduction**: Fix constant term a₀ real ≥ 0 via rotation z → e^{it/n}z, reducing search from 2(n−1) to 2n−3 real parameters
- **Branch-and-bound (n=3)**: Two-pass — coarse threshold at 50% of L*, then per-box Lipschitz upper bounds. 17,632 B&B evaluations; only the extremizer box survives
- **Lemniscate computation**: Marching squares on 2D grid, level-set arc length
- **Hessian check**: Central finite differences at extremizer (resolution 1600); all eigenvalues negative
- **Rigor note**: Conservative floating-point, not IEEE 1788 interval arithmetic. Margins are 100×–1000× larger than discretization errors, making conclusions robust

### Preprint

[Computational Verification of the EHP Conjecture for 3 ≤ n ≤ 8 (PDF)](papers/EHP_Bounds_Computational_Draft.pdf)

### Code and Results

```
scripts/erdos-114/
├── Cargo.toml
├── src/
│   ├── main.rs                  # n=3 branch-and-bound
│   └── bin/
│       ├── ehp_general.rs       # General n=3..8 verifier
│       ├── level2.rs            # Level 2 refinement
│       ├── level3_ia.rs         # Level 3 interval arithmetic
│       └── hessian_check.rs     # Hessian negativity at extremizer

results/erdos-114/
├── EHP_GENERAL_SUMMARY.json     # Full n=3..8 results with margins and verdicts
├── EHP_N3_LEVEL3_RESULTS.json   # n=3 branch-and-bound detail
└── EXP-MM-EHP-007-n{3..8}-inari_RESULTS.json  # Per-degree results with checksums
```

### Lean 4

A Lean 4 formalization sketch is in `scripts/erdos-114/EHP_N3.lean4`. Full 0-sorry formalizations for related bounds are available on request.

### References

- Tao, T. (2025). "The maximal length of the EHP lemniscate in high degree." arXiv:2512.12455
- Eremenko, A. & Hayman, W. (1999). On the length of lemniscates. *Michigan Math. J.* 46.
- Fryntov, A. & Nazarov, F. (2008). Local maximality of zⁿ − 1.
- Krishnapur, M., Lundberg, E. & Ramachandran, K. (2025). arXiv:2503.18270 (minimal area, dual problem)

---

## Other Experiments

| Problem | Script | Prize | Key Finding |
|---------|--------|-------|-------------|
| #86 C₄-free hypercube | `scripts/erdos-86/` | $100 | Exact ILP for n=2–5; ratio ex(Q_n,C₄)/e(Q_n) declines 0.75→0.583 for n=2–10 |
| #123 d-completeness | `scripts/erdos-123/` | $250 | Coprime triples containing 2 are d-complete to N=10,000; triples without 2 fail |
| #161 Discrepancy jumps | `scripts/erdos-161/` | $500 | F^(2)(n,α) step function; universal jump at α ≈ 1/3 across n=4–8 |
| #30 Sidon sets | `scripts/erdos-30/` | $1,000 | Erdős–Turán 1941 upper bound formalized in Lean 4, 0 sorry |

All results are in `results/` as JSON with SHA-256 checksums.

---

## Citation

```bibtex
@misc{mendoza2026erdos-experiments,
    author    = {Mendoza, Kenneth A.},
    title     = {Computational Experiments for Erd\H{o}s Problems},
    year      = {2026},
    url       = {https://github.com/MendozaLab/erdos-experiments},
    note      = {ORCID: 0009-0000-9475-5938}
}
```

## License

Apache 2.0 — see [LICENSE](LICENSE). Contributions welcome under the same license.
