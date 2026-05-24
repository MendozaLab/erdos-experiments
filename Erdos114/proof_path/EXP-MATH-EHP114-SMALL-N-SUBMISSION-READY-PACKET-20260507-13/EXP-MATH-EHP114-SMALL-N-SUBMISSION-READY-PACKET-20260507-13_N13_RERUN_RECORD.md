# n=13 re-run record

The re-run that v12's `_N13_RERUN_COMMAND.md` proposed was performed on 2026-05-07. This file records the outcome.

## What was found in v12

`EXP-MM-EHP-007-n13-inari_RESULTS.json` (2026-03-27 batch) reported `bb_total_evals = 0`, `bb_levels = []`, while passing the outer-domain and Hessian checks. Wall-clock was 18.16 s. The verdict was set to `EHP_N13_PROVEN` despite no branch-and-bound work having been done.

## Root cause

Logic bug in `scripts/erdos-114/src/bin/ehp_general_ieee1788.rs:987`. The verdict assembly required only `proof_complete && outer_safe && hess_all_neg`. `proof_complete` is set true in two places: when `boxes.is_empty()` at the START of an iteration (the buggy path — fires when `create_initial_boxes_recursive` returns empty), and when `ne_c == 0` after subdivision (the legitimate path). The original verdict logic could not distinguish the two.

## Bug fix

Commit [`dae62b8`](https://github.com/MendozaLab/erdos-experiments/commit/dae62b8) on `triage-recovered-2026-05-03`. The verdict guard now additionally requires `bb_total_evals > 0 && !level_log.is_empty()`. A distinct verdict `EHP_N{n}_INCOMPLETE_BB_NO_OP` is emitted when the BB phase short-circuits via empty initial boxes. Same commit moved the four 2026-03-27 zero-eval artifacts (`n=13, 14, 15, 16`) to `scripts/erdos-114/archive/exploratory_2026-03-27/` with a README explaining the scope.

## Re-run

Commit [`f89597e`](https://github.com/MendozaLab/erdos-experiments/commit/f89597e) on `triage-recovered-2026-05-03`. Released as [`v3.1.0`](https://github.com/MendozaLab/erdos-experiments/releases/tag/v3.1.0).

- Run host: Modal 32-CPU
- Wall-clock: **3,177.7 s** (52.96 minutes)
- Binary: `ehp_general_ieee1788` (post-`dae62b8` patched verdict logic)
- Config: degree 13, `coeff_bound = 3.0`, `grid_per_axis = 2`, `bb_res = 100`, `max_levels = 4`
- Output paths:
  - Canonical: `results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.json`
  - Versioned snapshot: `results/erdos-114/EXP-MM-EHP-007-n13-inari-RERUN-20260507-01/`

## Acceptance gate (from v12 `_N13_RERUN_COMMAND.md`)

| Field | Required | Actual | Pass |
|---|---|---|---|
| `verdict` | `EHP_N13_PROVEN` | `EHP_N13_PROVEN` | ✓ |
| `rigor` | `ieee_1788_interval_arithmetic_inari` | same | ✓ |
| `bb_proof_complete` | `true` | `true` | ✓ |
| `bb_total_evals` | > 0; in 10⁶–10⁹ window between n=12 (45M) and n=14 (855M) | **197,132,288** | ✓ |
| `bb_levels` length | > 0 | 1 (proof complete at level 0) | ✓ |
| `outer_domain_safe` | `true` | `true` | ✓ |
| `hessian_negative` | `true` | `true` | ✓ |
| `l_star_lower` | within `[28.85923995, 28.85924000]`, ≤ upper | `28.85923995588822` | ✓ |
| `l_star_upper` | within `[28.85923995, 28.85924000]`, ≥ lower | `28.859239955888235` | ✓ |

All nine fields pass. The level-0 trace shows 8,388,608 boxes, 4,194,304 eliminated, 4,194,304 extremizer-survivors, 0 non-extremizer survivors — the legitimate "all non-extremizer boxes eliminated at level 0" path that closes the proof.

## SHA-256

`c06c633b4053cdf2c4c6003327f30ee2ce683e6a9da70e047d5c5d43e685fd17`

Sidecar: `results/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.sha256` (two-column format, fixed in commit `98ea20c`).

## Audit trail

- Original anomalous artifact: `scripts/erdos-114/archive/exploratory_2026-03-27/EXP-MM-EHP-007-n13-inari_RESULTS.json` (preserved per Rule 3)
- Verdict-bug-fix: `dae62b8`
- Clean re-run promoted: `f89597e`
- GitHub release: `v3.1.0`
- Zenodo: corrected v3.1.0 packet published under version DOI 10.5281/zenodo.20087919 and existing concept DOI 10.5281/zenodo.19184467.

## Status

n=13 is fully certified and included in the v13 finite theorem statement (`1 ≤ n ≤ 14`, no exclusion). No further re-run is required.

## Tooling disclosure

This record was prepared with AI-assisted tooling. The mathematical claim rests only on the cited certificate artifact and SHA-256 sidecar. The author takes responsibility for the diagnosis, the verdict-bug-fix, and the re-run acceptance.
