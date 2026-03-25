# Erdős Problem Computational Experiments

Computational solvers and experimental results for open Erdős problems.

**Author**: Kenneth A. Mendoza ([ORCID 0009-0000-9475-5938](https://orcid.org/0009-0000-9475-5938))

---

## Problem #114 — Erdős–Herzog–Piranian Lemniscate Conjecture

**Statement**: Among all monic polynomials of degree n, the lemniscate {z ∈ ℂ : |p(z)| = 1} has maximal arc length when p(z) = zⁿ − 1.

**Prize**: $250 · **Status**: Open for small n (solved for sufficiently large n by Tao, Dec 2025)

### Background

Tao (arXiv:2512.12455) proved EHP for all sufficiently large n with an effective but tower-exponential threshold N₀ far beyond current computational reach. Eremenko–Hayman settled n = 2 analytically in 1999. This repository addresses n = 3 through n = 10.

### Results

All cases use IEEE 1788 certified interval arithmetic throughout via the `inari` Rust crate. Certified enclosures for L(zⁿ − 1) have widths ≲ 2.5 × 10⁻¹⁴.

| n | Dim. | L(zⁿ − 1) certified interval | Margin over 2nd best | B&B evals | Time |
|---|------|------------------------------|----------------------|-----------|------|
| 3 | 1 | [9.17972422234315, 9.17972422234317] | 17.1% | 8 | < 1 s |
| 4 | 3 | [11.0700205172566, 11.0700205172566] | 28.9% | 13,504 | 254 s |
| 5 | 5 | [13.0068113819187, 13.0068113819187] | 45.4% | 32 | < 1 s |
| 6 | 7 | [14.9657321896586, 14.9657321896586] | 54.3% | 128 | < 1 s |
| 7 | 9 | [16.9369006482509, 16.9369006482509] | 62.3% | 512 | 8 s |
| 8 | 11 | [18.9155531362862, 18.9155531362863] | 67.3% | 2,048 | 1 s |
| 9 | 13 | [20.8991118016671, 20.8991118016671] | 69.6% | 8,192 | 4 s |
| 10 | 15 | [22.8860603281654, 22.8860603281654] | 71.4% | 32,768 | 12 s |

**Dim.** is the reduced parameter space dimension 2n − 5 (symmetry-reduced from 2(n−1)).

**Margins increase monotonically** from 17.1% (n=3) to 71.4% (n=10), consistent with Tao's asymptotic result. Zero counterexamples found across the full search.

### Exact Formula (proven)

For the conjectured maximizer:

    L(zⁿ − 1) = 2^{1/n} · √π · Γ(1/(2n)) / Γ(1/(2n) + 1/2)

Verified against Python/mpmath at 50+ decimal digits for n = 3–10. Asymptotic: L(zⁿ − 1) = 2πn + 4 log 2 + O(1/n), consistent with Fryntov–Nazarov.

### Methodology

- **Symmetry reduction**: Translation (fix aₙ₋₁ = 0) + rotation (fix a₀ ∈ ℝ≥0) reduces search from 2(n−1) to 2n−5 real dimensions — a 1D search at n=3, 15D at n=10
- **Branch-and-bound**: Certified upper bound per box via Lipschitz constant + diameter; comparison u_i.sup() < L*.inf() uses directed rounding throughout. All non-extremizer boxes eliminated at level 0 for n ≥ 5
- **IEEE 1788 arithmetic**: Full certificate chain — Lipschitz bounds, fill distances, grid error envelopes, elimination comparisons — via `inari::Interval` with directed rounding
- **Hessian check**: All eigenvalues strictly negative at zⁿ − 1 (verified at four step sizes), confirming strict local optimality
- **Boundary verification**: L(p) < L(zⁿ − 1) confirmed on the frontier of the feasible region

### Computational Gap

Our results cover n ∈ {3, …, 10}. Tao's proof covers n ≥ N₀ (tower-exponential, N₀ ≫ 10¹⁰⁰). The gap n ∈ {11, …, N₀ − 1} remains open. Practical limit on current hardware: n ≲ 15–18.

### Preprint

[Computational Verification of the EHP Conjecture for 3 ≤ n ≤ 10 (PDF)](papers/EHP_Bounds_Computational_Draft.pdf)

### Code and Results

```
scripts/erdos-114/
├── Cargo.toml
├── src/
│   ├── main.rs                      # n=3 branch-and-bound
│   └── bin/
│       ├── ehp_general.rs           # General n=3..8 verifier
│       ├── ehp_general_ieee1788.rs  # IEEE 1788 certified verifier (n=3..10)
│       ├── ehp_n3_ieee1788.rs       # n=3 certified B&B
│       ├── level2.rs / level3_ia.rs # Refinement stages
│       └── hessian_check.rs         # Hessian negativity at extremizer

results/erdos-114/
├── EHP_GENERAL_SUMMARY.json         # n=3..8 results
├── EXP-MM-EHP-007-n{3..8}-inari_RESULTS.json  # Per-degree certified results
```

Reproducible via: `cargo run --release --bin ehp_general_ieee1788`

### References

- Tao, T. (2025). "The Erdős–Herzog–Piranian conjecture for large polynomial degrees." arXiv:2512.12455
- Eremenko, A. & Hayman, W. (1999). On the length of lemniscates. *Michigan Math. J.* 46.
- Fryntov, A. & Nazarov, F. (2025). On the local extremality of lemniscates of zⁿ − 1.
- Krishnapur, M., Lundberg, E. & Ramachandran, K. (2025). arXiv:2503.18270 (minimal area, dual problem)
- Tanaka, M. inari: IEEE 1788 interval arithmetic for Rust. https://crates.io/crates/inari

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
    title     = {Computational Verification of the {Erd\H{o}s--Herzog--Piranian} Conjecture
                 for Degrees $3 \leq n \leq 10$},
    year      = {2026},
    url       = {https://github.com/MendozaLab/erdos-experiments},
    note      = {ORCID: 0009-0000-9475-5938}
}
```

## License

Apache 2.0 — see [LICENSE](LICENSE). Contributions welcome under the same license.
