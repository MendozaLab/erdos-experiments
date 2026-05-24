# EHP114 finite-degree submission packet (v13)

## Scope

Submission preflight packet for the finite-degree variant of the Erdős–Herzog–Piranian lemniscate-length conjecture. Covers `1 ≤ n ≤ 14` contiguously. Does not assert the all-degree Erdős #114 conjecture.

## What changed from v12

| Change | Reason |
|---|---|
| Theorem statement is contiguous (`1 ≤ n ≤ 14`, no `n ≠ 13` exclusion) | n=13 was re-run cleanly on 2026-05-07 (197M B&B evals, 53min on Modal 32-CPU); it's no longer quarantined |
| Added n=13 axiom to the Lean stub | Now part of the certified set |
| `_ZENODO_STRATEGY.md` flips: new Zenodo version **is** appropriate | A certificate row changed; v3.1.0 was manually published to Zenodo on 2026-05-08 as version DOI `10.5281/zenodo.20087919` after the GitHub auto-archive did not appear |
| `_N13_RERUN_RECORD.md` (replaces v12's `_N13_RERUN_COMMAND.md`) | Re-run was performed; this records what happened, not what to do next |
| Verdict-bug fix in the engine | Commit `dae62b8` — `EHP_N{n}_PROVEN` now requires `bb_total_evals > 0 && !level_log.is_empty()` |
| Sidecar format fix | n=13 sidecar rewritten to two-column form (commit `98ea20c`) |

v12 stays as the historical record (FROZEN-locked under its 2026-05-07-r01 redteam round). It is not deleted, not edited. v13 is the now-current packet.

## Finite theorem

For every monic complex polynomial `p` with `deg p = n` exactly, where

```
n ∈ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
```

the lemniscate length `L({z ∈ ℂ : |p(z)| = 1})` is at most the lemniscate length of `z^n - 1`.

## Dependency split

| Range | Source | Note |
|---|---|---|
| `n = 1` | Direct analytic | Monic linear polynomial gives a translated unit circle of length `2π` |
| `n = 2` | MacLane (1953/54), accessible pin Eremenko–Hayman (1999) | Literature row, no local computation |
| `3 ≤ n ≤ 14` | Rust + `inari` IEEE-1788 interval branch-and-bound | DOI-backed result JSONs with SHA-256 sidecars; n=13 re-run on 2026-05-07 |

## n=13 narrative

The n=13 row in v12 reported `bb_total_evals = 0` and an empty `bb_levels` array — the symptom of a verdict-logic bug at high reduced dimension. The bug was traced to `ehp_general_ieee1788.rs:987` (the verdict could not distinguish "BB exhaustively eliminated everything" from "BB had no boxes to start") and patched in commit `dae62b8`. The patched binary was re-run against n=13 on 2026-05-07; it closed the proof at level 0 with 4,194,304 non-extremizer boxes eliminated and 0 surviving, after 197,132,288 box evaluations and 52.96 minutes of wall-clock on a 32-CPU Modal worker. Released as `v3.1.0`. Full audit trail in [`_N13_RERUN_RECORD.md`](EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-13_N13_RERUN_RECORD.md).

n=13 is part of the certified set in v13.

## Public references

- Zenodo certificate record: version DOI [10.5281/zenodo.20087919](https://doi.org/10.5281/zenodo.20087919); concept DOI [10.5281/zenodo.19184467](https://doi.org/10.5281/zenodo.19184467) resolves to the latest version. This v3.1.0 version was manually created from latest record `19480329` on 2026-05-08 after the GitHub auto-archive did not appear.
- Code and results repository: <https://github.com/MendozaLab/erdos-experiments>
- v3.1.0 release: <https://github.com/MendozaLab/erdos-experiments/releases/tag/v3.1.0>
- Erdős Problems #114 page: <https://www.erdosproblems.com/114>
- Google DeepMind formal-conjectures: <https://github.com/google-deepmind/formal-conjectures>
- Tao large-degree paper: [arXiv:2512.12455](https://arxiv.org/abs/2512.12455)

## Reviewer-grade certificate table

| n | dependency | public source | SHA-256 | status |
|---:|---|---|---|---|
| 1 | direct analytic | translated unit circle | n/a | covered |
| 2 | MacLane / Eremenko–Hayman literature row | [MacLane DOI](https://doi.org/10.1307/mmj/1028989918) / [Eremenko–Hayman arXiv](https://arxiv.org/abs/0805.2295) / [PDF](https://www.math.purdue.edu/~eremenko/dvi/erdos23.pdf) | n/a | covered |
| 3 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n3-inari_RESULTS.json) | `a884d1bfec1563f6e6f7ae4cbb2ec607b43be033d06c41d14782459e67ec2b95` | IEEE-1788 interval-certified |
| 4 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n4-inari_RESULTS.json) | `0924dd7424d2615099ff95d47cb4c120ba22e907adaa9af881cded9678241209` | IEEE-1788 interval-certified |
| 5 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n5-inari_RESULTS.json) | `21ca3c7607dc1fbb7b08982666f4620dff808fdc581eecae9967c51fafb05447` | IEEE-1788 interval-certified |
| 6 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n6-inari_RESULTS.json) | `41ac3027e9ae5add9e1208c0faa5897c36762d69b5eeccef96068c96af567b3d` | IEEE-1788 interval-certified |
| 7 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n7-inari_RESULTS.json) | `832ddaf219d717e275ee95c01f271dd3120e255cf81f3aa72ba3c2f56ff84054` | IEEE-1788 interval-certified |
| 8 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n8-inari_RESULTS.json) | `c7a1fd80fbfed1efd53eaa35283e467994e3a4175541f3817485d9551d14dcdc` | IEEE-1788 interval-certified |
| 9 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n9-inari_RESULTS.json) | `5bc2887826c9ef21752c115c9a4a2ab983f94eea06f0111ea98db454fe1358f4` | IEEE-1788 interval-certified |
| 10 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n10-inari_RESULTS.json) | `2b72e052aa7200f7ac5d40992843601988de234093c44860fde99a7871e19581` | IEEE-1788 interval-certified |
| 11 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n11-inari_RESULTS.json) | `67f20cce1d3d54cad2d6bc708ab9ec796c17cb4d34a9728f2954b8a6cbf7c89c` | IEEE-1788 interval-certified |
| 12 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n12-inari_RESULTS.json) | `42a517997445d158649feefae2b7287bc9b548c6391796df0c1246489c6aa064` | IEEE-1788 interval-certified |
| **13** | **IEEE-1788 interval certificate (re-run 2026-05-07)** | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.json) | `c06c633b4053cdf2c4c6003327f30ee2ce683e6a9da70e047d5c5d43e685fd17` | **IEEE-1788 interval-certified (197M evals, 53min)** |
| 14 | IEEE-1788 interval certificate | [JSON](https://github.com/MendozaLab/erdos-experiments/blob/main/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json) | `50b1c965c842ced25b2930c2b71ffb6e2da693872aa464a19fbd9d5d9efa0ca7` | IEEE-1788 interval-certified |

## Local preflight status

- Packet status: `SMALL_N_PRESENTATION_PACKET_READY_FULL_RANGE_1_TO_14`
- Covered degrees: `14` (all of n=1..14)
- Quarantined degrees: `0`
- Source SHA failures: `0`
- Full EHP114 claim: `false` (we still do not assert the all-degree conjecture)
- External submission authorized: `false` (manual user confirmation gates remain)

## Process gates required before public posting

In order:

1. **prepub-redteam** — `init … --new-round`, re-import the Perplexity quorum critique (still applicable), log edits, freeze, verify.
2. **Publisher Gate** — `anthropic-skills:publisher` over the v13 frozen post and the assembled paper at `Math/preprints/ehp114-finite/`.
3. **Crackpot-Scrub** — over the same artifacts.
4. **Zenodo state check** — complete. Version DOI `10.5281/zenodo.20087919` now pins the corrected v3.1.0 certificate packet; concept DOI `10.5281/zenodo.19184467` resolves to latest.

## Tooling disclosure

This packet and its draft wording were prepared with AI-assisted tooling. The mathematical claim rests only on the cited literature rows and SHA-checked certificate artifacts; no AI output is used as proof evidence. The author takes responsibility for the mathematical claims, certificate selection, and submission wording.
