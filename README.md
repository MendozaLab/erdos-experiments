# Erdős Problem Computational Experiments

Computational solvers and experimental results for open Erdős problems.

**Author**: Kenneth A. Mendoza ([ORCID 0009-0000-9475-5938](https://orcid.org/0009-0000-9475-5938))

## Motivation

Many classical Erdős problems have resisted proof for decades, but their small cases are computationally tractable. This repository collects experimental evidence — verified bounds, ratio estimates, and structural data — for five open problems whose prize values range from $100 to $1,000. The goal is to expose patterns that guide future proofs and to provide reproducible baselines for the community.

## Preprints / Manuscripts

- **[Draft Paper] Computational Verification of the EHP Conjecture for 3 ≤ n ≤ 8**
  [Read PDF Draft Here](papers/EHP_Bounds_Computational_Draft.pdf)

## Experiments

| Problem | Link | Script | Prize | Key Finding |
|---------|------|--------|-------|-------------|
| #114 EHP Conj Bounds | [erdosproblems.com/114](https://www.erdosproblems.com/114) | `scripts/erdos-114/` | $1,000 | Computationally verified n=3 through n=8 bounds |
| #86 C4-free hypercube | — | `scripts/erdos-86/` | $100 | Ratio ex(Q_n,C4)/e(Q_n) trends toward 1/2 |
| #123 d-completeness | — | `scripts/erdos-123/` | $250 | All coprime triples with 2 are d-complete to N=10000 |
| #161 Discrepancy jumps | — | `scripts/erdos-161/` | $500 | F^(2)(n,α) is a step function with universal jump at α=1/3 |
| #30 Sidon sets | — | `scripts/erdos-30/` | $1,000 | Greedy Sidon, Singer construction, WalkSAT for W(k) in Rust |

## Getting Started

### Requirements

- **Rust** ≥ 1.75 (for Sidon set solvers and performance-critical code)
- **Python** ≥ 3.10 (for analysis scripts and plotting)
- Standard crates and pip packages are listed in each script directory

### Running an experiment

```bash
# Example: run the EHP 114 solver
cd scripts/erdos-114/
cargo run --release   # Rust solver
python analyze.py     # Python post-processing
```

## Results

All results are in `results/` as JSON with SHA-256 checksums for reproducibility. Each JSON file contains the problem number, parameters, computed values, and a timestamp.

## Lean 4 Proofs

Formal machine-verified proofs for selected bounds are being developed in a companion repository (currently private — available on request). These use Lean 4 with Mathlib and target zero `sorry` stubs.

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

Apache 2.0 — see [LICENSE](LICENSE) for details. Contributions are welcome under the same license.
