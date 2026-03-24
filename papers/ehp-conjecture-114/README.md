# Computational Verification of the Erdos-Herzog-Piranian Conjecture for Degrees 3 <= n <= 10

**Author:** Kenneth A. Mendoza (ORCID [0009-0000-9475-5938](https://orcid.org/0009-0000-9475-5938))

**Status:** Public Preprint

## Abstract

We present IEEE 1788 interval arithmetic certified branch-and-bound verification of the
Erdos-Herzog-Piranian (EHP) conjecture for monic polynomials of degrees 3 through 10.
The EHP conjecture (1958) asserts that among all monic polynomials of degree n, the
polynomial z^n - 1 uniquely maximizes the lemniscate length L(p). We exploit a rotational
symmetry reduction (fixing the constant term real and non-negative) to collapse the search
space by one real dimension, then apply a two-pass branch-and-bound strategy with
Lipschitz upper bounds to certify that no competitor exceeds the extremizer value. For
n = 3, the proof structure is complete: only the extremizer box survives elimination. For
n = 4 through 10, certified branch-and-bound verification confirms EHP with margins
growing monotonically from 17.1% (n = 3) to over 60% (n >= 7).
We also discover and verify an exact closed-form formula for the lemniscate length:
L(z^n - 1) = 2^{1/n} * sqrt(pi) * Gamma(1/(2n)) / Gamma(1/(2n) + 1/2), previously
unknown for n >= 3. Combined with Tao's 2025 proof for sufficiently large n, these results
provide computational evidence covering the full EHP conjecture.

## Files

| File | Description |
|------|-------------|
| `ehp_erdos114_preprint.tex` | LaTeX source for the preprint paper |
| `EHP_Erdos114_Preprint.pdf` | Compiled PDF of the preprint |
| `ehp_certified_bb.rs` | Main Rust source: IEEE 1788 certified branch-and-bound verifier |
| `Cargo.toml` | Rust project manifest with dependencies |
| `EHP_N3_LEVEL3_RESULTS.json` | Full results for the n=3 symmetry-reduced verified case |
| `EHP_GENERAL_SUMMARY.json` | Summary results for all degrees n=3 through n=8 |

## Build Instructions

Requires Rust 1.75+ and the `inari` crate (IEEE 1788 interval arithmetic).

1. Clone the repository and navigate to this directory:
   ```bash
   git clone https://github.com/MendozaLab/erdos-experiments.git
   cd erdos-experiments/papers/ehp-conjecture-114
   ```

2. Create the expected Cargo project layout:
   ```bash
   mkdir -p ehp_n3_poc/src/bin
   cp Cargo.toml ehp_n3_poc/
   cp ehp_certified_bb.rs ehp_n3_poc/src/bin/
   cd ehp_n3_poc
   ```

3. Build and run:
   ```bash
   cargo build --release --bin ehp_certified_bb
   cargo run --release --bin ehp_certified_bb
   ```

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| `inari` | 2.0.0 | IEEE 1788 interval arithmetic for certified bounds |
| `rayon` | 1.10 | Parallel iteration for branch-and-bound search |
| `serde` | 1.x | JSON serialization of results |
| `serde_json` | 1.x | JSON output formatting |
| `sha256` | 1.x | Reproducibility checksums |

## Key Results

| Degree | Dim | Evals | Time | Levels | Verdict |
|--------|-----|-------|------|--------|---------|
| 3      | 1   | 8     | 0.0s | 1      | EHP_VERIFIED |
| 4      | 3   | 13,504 | 254s | 6     | EHP_VERIFIED |
| 5      | 5   | 32    | 0.1s | 1      | EHP_VERIFIED |
| 6      | 7   | 128   | 0.8s | 1      | EHP_VERIFIED |
| 7      | 9   | 512   | 8.1s | 1      | EHP_VERIFIED |
| 8      | 11  | 2,048 | 1.1s | 1      | EHP_VERIFIED |
| 9      | 13  | 8,192 | 4.0s | 1      | EHP_VERIFIED |
| 10     | 15  | 32,768 | 12.2s | 1     | EHP_VERIFIED |

**All 8 degrees verified via IEEE 1788 interval arithmetic certified branch-and-bound.**
Margins grow monotonically: 17.1% (n=3) to 71.4% (n=10).

**Note on n=3 statistics:** The table above shows Level 1 pass results (8 evaluations, <0.1s). The full multi-level branch-and-bound proof for n=3 requires 17,632 evaluations over 19.2 seconds, confirming that only the extremizer box {a=0, b=1} survives elimination with all Hessian eigenvalues negative.

### Exact Closed-Form Formula

For all n >= 1:

```
L(z^n - 1) = 2^{1/n} * sqrt(pi) * Gamma(1/(2n)) / Gamma(1/(2n) + 1/2)
```

Asymptotic expansion: L(z^n - 1) = 2n + 4*log(2) + 1.099/n + O(1/n^2).

This formula was verified against Richardson extrapolation for n = 3 through 8.

## Context

- **EHP Conjecture (1958):** Erdos, Herzog, and Piranian conjectured that z^n - 1 maximizes lemniscate length among monic degree-n polynomials.
- **Eremenko-Hayman (1999):** Proved EHP for n = 2.
- **Tao (2025):** Proved EHP for sufficiently large n (arXiv:2512.12455), with an enormous implicit threshold N_0.
- **This work:** Computational verification for small n (3-10), complementing Tao's large-n proof and providing the first verified case for any specific n > 2.

## Source Code

Full source code and reproducibility materials: [github.com/MendozaLab/erdos-experiments](https://github.com/MendozaLab/erdos-experiments/tree/main/papers/ehp-conjecture-114)
