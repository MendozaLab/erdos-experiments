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
- `n=2`: Eremenko-Hayman literature row.
- `3 <= n <= 14`: DOI-backed interval-certificate rows.

## Local preflight status

- Finite packet status: `FINITE_N_LESS_15_PROOF_PACKET_PASS_NOT_FULL_EHP_PROOF`
- Covered degrees: `14`
- Missing degrees: `0`
- Source SHA failures: `0`
- Full EHP114 claim: `false`
- External submission authorized: `false`

## Claim ceiling

Submission preflight packet only. It prepares erdosproblems.com and Google DeepMind formal-conjectures wording for the finite n<15 result and the conditional Tao-threshold bridge; it does not assert the all-degree Erdos #114 statement and does not authorize n=15 or higher computation.
