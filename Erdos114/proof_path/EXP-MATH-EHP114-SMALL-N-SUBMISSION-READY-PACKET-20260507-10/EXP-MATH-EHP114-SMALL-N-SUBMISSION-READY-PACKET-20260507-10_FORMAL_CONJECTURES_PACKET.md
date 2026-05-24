# Formal-conjectures packet for Erdos Problem 114

## Intended PR title

`Erdos114: add finite n < 15 certificate theorem`

## Intended status change

Do **not** mark the main theorem `erdos_114` as closed. Keep the main statement open unless and until the Tao threshold bridge is explicit.

## Proposed formal-conjectures shape

- Keep namespace: `Erdos114`.
- Keep the main theorem as `@[category research open, AMS 30] theorem erdos_114 ... := by sorry`.
- Add or retain the finite variant `erdos_114_finite_lt_15`.
- Mark the finite variant with the repository's resolved-research category only if the certificate axioms are accepted as explicit external dependencies.
- Keep every machine-certified row as an explicit axiom or external-certificate lemma; do not hide computational certification inside a proof term.

## Suggested docstring wording

"Finite certified range for Erdos Problem 114. For `1 <= n < 15`, the Erdos-Herzog-Piranian lemniscate inequality holds, using the direct `n=1` case, the MacLane / Eremenko-Hayman degree-two case, and DOI-backed Rust/inari IEEE-1788 interval certificates for `3 <= n <= 14`. The historical `n=2` source is MacLane, and the accessible source pin is Alexandre Eremenko and Walter K. Hayman, "On the length of lemniscates", Michigan Math. J. 46 (1999), 409-415; arXiv:0805.2295. This finite statement is separate from the open all-degree conjecture and from Tao's sufficiently-large-`n` theorem. No explicit threshold is extracted from Tao's sufficiently-large-n theorem in this finite packet, and no degree n >= 15 is claimed here."

## PR body

This PR records a finite variant of Erdos Problem 114 rather than changing the status of the main conjecture. The finite theorem covers exactly `1 <= n < 15`.
The proof dependencies are typed explicitly:

- direct analytic lemma for `n=1`;
- literature-only MacLane / Eremenko-Hayman input for `n=2`, with MacLane as historical citation [G. R. MacLane, On a conjecture of Erdos, Herzog, and Piranian, Michigan Math. J. 2 (1953/54), 147-148, doi:10.1307/mmj/1028989918.](https://doi.org/10.1307/mmj/1028989918) and Eremenko-Hayman pinned to Michigan Math. J. 46 (1999), 409-415 / arXiv:0805.2295;
- DOI-backed Rust/inari IEEE-1788 interval certificates for `3 <= n <= 14`, with row-level result links and SHA-256 sidecars;
- Tao's sufficiently-large-`n` theorem is mentioned only as context for the remaining bridge question. No explicit threshold is extracted from Tao's sufficiently-large-n theorem in this finite packet, and no degree n >= 15 is claimed here.

## Reviewer-grade certificate table

| n | dependency | public source | SHA / pin | status |
|---:|---|---|---|---|
| 1 | direct analytic | translated unit circle | not applicable | covered |
| 2 | MacLane / Eremenko-Hayman literature row | [MacLane DOI](https://doi.org/10.1307/mmj/1028989918) / [Eremenko-Hayman arXiv](https://arxiv.org/abs/0805.2295) / [PDF](https://www.math.purdue.edu/~eremenko/dvi/erdos23.pdf) | G. R. MacLane, On a conjecture of Erdos, Herzog, and Piranian, Michigan Math. J. 2 (1953/54), 147-148, doi:10.1307/mmj/1028989918.; Eremenko-Hayman accessible pin: Michigan Math. J. 46 (1999), 409-415; degree-two extremal case pinned to Bernoulli lemniscate | covered |
| 3 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n3-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n3-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n3-inari_RESULTS.json) | `a884d1bfec1563f6e6f7ae4cbb2ec607b43be033d06c41d14782459e67ec2b95`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n3-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 4 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n4-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n4-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n4-inari_RESULTS.json) | `0924dd7424d2615099ff95d47cb4c120ba22e907adaa9af881cded9678241209`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n4-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 5 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n5-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n5-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n5-inari_RESULTS.json) | `21ca3c7607dc1fbb7b08982666f4620dff808fdc581eecae9967c51fafb05447`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n5-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 6 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n6-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n6-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n6-inari_RESULTS.json) | `41ac3027e9ae5add9e1208c0faa5897c36762d69b5eeccef96068c96af567b3d`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n6-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 7 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n7-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n7-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n7-inari_RESULTS.json) | `832ddaf219d717e275ee95c01f271dd3120e255cf81f3aa72ba3c2f56ff84054`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n7-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 8 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n8-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n8-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n8-inari_RESULTS.json) | `c7a1fd80fbfed1efd53eaa35283e467994e3a4175541f3817485d9551d14dcdc`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n8-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 9 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n9-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n9-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n9-inari_RESULTS.json) | `5bc2887826c9ef21752c115c9a4a2ab983f94eea06f0111ea98db454fe1358f4`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n9-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 10 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n10-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n10-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n10-inari_RESULTS.json) | `2b72e052aa7200f7ac5d40992843601988de234093c44860fde99a7871e19581`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n10-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 11 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n11-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n11-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n11-inari_RESULTS.json) | `67f20cce1d3d54cad2d6bc708ab9ec796c17cb4d34a9728f2954b8a6cbf7c89c`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n11-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 12 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n12-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n12-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n12-inari_RESULTS.json) | `42a517997445d158649feefae2b7287bc9b548c6391796df0c1246489c6aa064`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n12-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |
| 13 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n13-inari_RESULTS.json) | `a4e72a9be2811e9d2290c6cdd0f6a9f1dd17fa85790e380bd5af8b09179e02ac`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED route exception flagged |
| 14 | DOI-backed Rust/inari IEEE-1788 interval certificate | [result JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json) / [raw](https://raw.githubusercontent.com/MendozaLab/erdos-experiments/main/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json) / [Zenodo file](https://zenodo.org/records/19480329/files/EXP-MM-EHP-007-n14-inari_RESULTS.json) | `50b1c965c842ced25b2930c2b71ffb6e2da693872aa464a19fbd9d5d9efa0ca7`; [sha sidecar](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.sha256) | LOCAL_INTERVAL_CERTIFIED |

The public certificate record is [10.5281/zenodo.19480329](https://zenodo.org/records/19480329); the code/results repository is [https://github.com/MendozaLab/erdos-experiments](https://github.com/MendozaLab/erdos-experiments).  Local preflight checks report missing degree rows `0`, source SHA failures `0`, and `full_ehp114_claim = false`.

## Review question

Is this the right formal-conjectures shape: keep `erdos_114` open, add the finite `n < 15` variant, and leave the Tao-threshold bridge outside this PR until an explicit threshold is available?

## Claim ceiling

Submission preflight packet only. It prepares erdosproblems.com and Google DeepMind formal-conjectures wording for the finite n<15 result and the conditional Tao-threshold bridge; it does not assert the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.

## Tooling disclosure

Tooling disclosure: this packet and draft wording were prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature rows and SHA-checked certificate artifacts; no AI output is used as proof evidence.
