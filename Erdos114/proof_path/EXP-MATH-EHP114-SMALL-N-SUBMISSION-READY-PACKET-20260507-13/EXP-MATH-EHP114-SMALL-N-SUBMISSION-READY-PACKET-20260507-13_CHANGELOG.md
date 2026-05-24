# v13 changelog (vs v12)

| File / area | v12 | v13 |
|---|---|---|
| Finite theorem hypothesis | `1 ≤ n ≤ 14, n ≠ 13` (n=13 quarantined) | `1 ≤ n ≤ 14` contiguous |
| Certified degrees | `{1, 2, 3, ..., 12, 14}` (13 degrees) | `{1, 2, 3, ..., 14}` (14 degrees, contiguous) |
| n=13 in public certificate table | quarantined; row preserved in Zenodo for transparency | included; SHA `c06c633b…b685fd17`, 197M B&B evals, 53min wall-clock |
| Lean stub | n=13 explicitly absent with TODO comment; `erdos_114_finite_le_14_except_13` | n=13 axiom present; theorem renamed `erdos_114_finite_le_14`; `(hp_deg_ne_13 : p.natDegree ≠ 13)` hypothesis dropped |
| `_N13_RERUN_COMMAND.md` | re-run protocol + acceptance gate (forward-looking) | replaced by `_N13_RERUN_RECORD.md` (records the actual re-run, audit trail, acceptance gate results) |
| Zenodo decision | edit v5 description; v6 only if certain conditions met | v6 mint **is** appropriate per the same memo's conditions (cert row changed); v3.1.0 GitHub release likely auto-triggered v6 mint already |
| Engine | bug at `ehp_general_ieee1788.rs:987` (verdict could not distinguish empty initial boxes from exhaustive elimination) | patched in commit `dae62b8`; verdict now requires `bb_total_evals > 0 && !level_log.is_empty()`; new `EHP_N{n}_INCOMPLETE_BB_NO_OP` verdict for the bug's failure mode |
| Anomalous artifact location | `scripts/erdos-114/EXP-MM-EHP-007-n{13,14,15,16}-inari_RESULTS.{json,sha256}` (8 files) | moved to `scripts/erdos-114/archive/exploratory_2026-03-27/` with README explaining scope (commit `dae62b8`) |
| Sidecar format | n=13 sidecar was bare hash (broke `shasum -c`) | rewritten to two-column format (commit `98ea20c`) |
| Public release | n/a | `v3.1.0` tag on `triage-recovered-2026-05-03`, GitHub release with Block-G AI-acknowledgment |
| File count in packet | 12 | 12 (one removed: `_N13_RERUN_COMMAND.md`; one added: `_N13_RERUN_RECORD.md`; net 0) |

## Issues v13 still does not resolve

- The `inari` IEEE-1788 interval-arithmetic implementation itself is not formally verified. A reviewer who does not trust the toolchain has no recourse from this packet alone.
- The Tao threshold extraction is not attempted. The bridge between the finite packet and the all-degree statement remains open.
- The deeper question of why `create_initial_boxes_recursive` returned empty at high reduced dimension (the upstream cause of the verdict bug's symptom) is unidentified. The verdict-bug fix prevents the spurious `_PROVEN` verdict but does not explain or fix the empty-boxes condition. Tracked in `scripts/erdos-114/archive/exploratory_2026-03-27/README.md`.
- No external endorser path is established for an arXiv submission; per memory, Moree is dead, Baez is gated.
- n=15 and n=16 remain out of production scope. Re-running them is not part of v13 because the box-construction issue at high d is unresolved.

## What v13 supersedes

v13 supersedes v12 forward — v12's narrative ("n=13 quarantine") is no longer the current claim. v12 stays on disk as the historical packet under its own FROZEN.lock; it is not edited or deleted.
