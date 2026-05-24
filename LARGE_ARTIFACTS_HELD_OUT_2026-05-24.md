# Large Artifacts Held Out Of Normal Git — 2026-05-24

This preservation branch uses a 95 MB normal-Git cutoff. The artifacts below are
kept on local disk but intentionally not committed as ordinary Git blobs. Route
them through Git LFS, R2, Zenodo, or another approved artifact store before any
cleanup that would remove local copies.

| SHA-256 | Bytes | Path | Reason |
|---|---:|---|---|
| `05b69c24376abbe7c9026b38fc6e4891bc880823dbf08a8012037623363e897f` | 172,748,872 | `Erdos114/EXP-MATH-RH-BN-HETEROGENEOUS-COMPUTE-MDL-HARDENING-20260507-01_RESULTS.json` | Over 95 MB |
| `07a0e77c2f564b0cc9bc4bf195a3eb90ae91a6415b425d6d77fe6e4670a6473e` | 116,219,324 | `Erdos114/validated_length/EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.json` | Over 95 MB |
| `604bb09c51df0c6bc5bd44afc7a43574a585f00694668f955a97e88072519f59` | 190,080,142 | `Erdos114/validated_length/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.json` | Over 95 MB |
| `6e087a3804de56db0e926556337c71d9b09c66584c5105b63bcbaa6528e54c6d` | 467,854,755 | `results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-BRANCH-59-64-2026-05-01_RESULTS.json` | Over 95 MB |
| `0c8505515f3b783b14e6642226ce9b9d459d2113eb75f3c59029e116d680036b` | 150,550,855 | `results/erdos-30/EXP-MM-030-PMF-GROUNDFACE-BRANCH-71-2026-05-01_RESULTS.json` | Over 95 MB |
| `fce4298447ef3f13ea0009fda6930d8ca94818454ce6b2e7f1a61c808f83d6bc` | 122,069,987 | `scripts/erdos-114/bridge-diagnostic-taylor-lipschitz-output-z16/EXP-MATH-EHP114-N14-EPS01-WORST-SUBCELL-BRIDGE-DIAGNOSTIC-20260505-01_RESULTS.json` | Over 95 MB |

Generated Lean/Lake cache packs under `.lake/` were also found above 95 MB.
Those are build caches, already ignored, and are not preservation artifacts.
