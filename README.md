# Erdős Problem Computational Experiments

Computational solvers and experimental results for open Erdős problems.

**Author**: Kenneth A. Mendoza ([ORCID 0009-0000-9475-5938](https://orcid.org/0009-0000-9475-5938))

## Preprints / Manuscripts

- **[Draft Paper] Computational Verification of the EHP Conjecture for 3 ≤ n ≤ 8** 
  [Read PDF Draft Here](papers/EHP_Bounds_Computational_Draft.pdf)

## Experiments

| Problem | Script | Prize | Key Finding |
|---------|--------|-------|-------------|
| #114 EHP Conj Bounds | `scripts/erdos-114/` | $1,000 | Computationally verified n=3 through n=8 bounds |
| #86 C4-free hypercube | `scripts/erdos-86/` | $100 | Ratio ex(Q_n,C4)/e(Q_n) trends toward 1/2 |
| #123 d-completeness | `scripts/erdos-123/` | $250 | All coprime triples with 2 are d-complete to N=10000 |
| #161 Discrepancy jumps | `scripts/erdos-161/` | $500 | F^(2)(n,α) is a step function with universal jump at α=1/3 |
| #30 Sidon sets | `scripts/erdos-30/` | $1,000 | Greedy Sidon, Singer construction, WalkSAT for W(k) in Rust |

## Results

All results in `results/` as JSON with SHA-256 checksums.

## Lean 4 Proofs

Formal proofs are in the companion repo: [MendozaLab/erdos-lean4](https://github.com/MendozaLab/erdos-lean4)

## Citation

```bibtex
@misc{mendoza2026erdos-experiments,
  author = {Mendoza, Kenneth A.},
  title = {Computational Experiments for Erdős Problems},
  year = {2026},
  url = {https://github.com/MendozaLab/erdos-experiments},
  note = {ORCID: 0009-0000-9475-5938}
}
```

## License

Apache 2.0
