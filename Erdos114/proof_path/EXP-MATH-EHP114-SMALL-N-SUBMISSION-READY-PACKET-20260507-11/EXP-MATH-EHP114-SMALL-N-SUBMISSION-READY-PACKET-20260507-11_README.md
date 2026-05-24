# EHP114 finite n < 15 submission preflight packet

## Finite theorem

For every monic complex polynomial p of degree n with 1 <= n < 15, the lemniscate length is at most the length for p(z)=z^n-1.

## Public references

- Zenodo certificate record: [10.5281/zenodo.19480329](https://zenodo.org/records/19480329)
- Code and results repository: [https://github.com/MendozaLab/erdos-experiments](https://github.com/MendozaLab/erdos-experiments)
- Erdős Problems #114 page: [https://www.erdosproblems.com/114](https://www.erdosproblems.com/114)
- Google DeepMind formal-conjectures: [https://github.com/google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures)
- Tao large-degree paper: [https://arxiv.org/abs/2512.12455](https://arxiv.org/abs/2512.12455)

## Dependency split

- `n=1`: direct analytic linear case.
- `n=2`: MacLane / Eremenko-Hayman degree-two literature row. MacLane is the historical citation ([doi:10.1307/mmj/1028989918](https://doi.org/10.1307/mmj/1028989918)); Eremenko-Hayman is the accessible source pin at Michigan Math. J. 46 (1999), 409-415 / arXiv:0805.2295.
- `3 <= n <= 14`: DOI-backed Rust/inari IEEE-1788 interval-certificate rows with row-level result links and SHA-256 hashes.
- Tao bridge status: No explicit threshold is extracted from Tao's sufficiently-large-n theorem in this finite packet, and no degree n >= 15 is claimed here.

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

## Local preflight status

- Finite packet status: `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`
- Covered degrees: `14`
- Missing degrees: `0`
- Source SHA failures: `0`
- Full EHP114 claim: `false`
- External submission authorized: `false`

## Claim ceiling

Submission preflight packet only. It prepares erdosproblems.com and Google DeepMind formal-conjectures wording for the finite n<15 result and the conditional Tao-threshold bridge; it does not assert the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.

## Tooling disclosure

Tooling disclosure: this packet and draft wording were prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature rows and SHA-checked certificate artifacts; no AI output is used as proof evidence. The author takes responsibility for the mathematical claims, certificate selection, and submission wording.
