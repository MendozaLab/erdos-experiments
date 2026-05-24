# n=13 re-run command and acceptance gate

## Diagnosis

The current `EXP-MM-EHP-007-n13-inari_RESULTS.json` payload is:

```json
{
  "experiment": "EXP-MM-EHP-007-n13-inari",
  "degree": 13,
  "reduced_dim": 23,
  "verdict": "EHP_N13_PROVEN",
  "rigor": "ieee_1788_interval_arithmetic_inari",
  "l_star_lower": 28.85923995588822,
  "l_star_upper": 28.859239955888235,
  "bb_proof_complete": true,
  "bb_total_evals": 0,
  "bb_levels": [],
  "outer_domain_safe": true,
  "hessian_negative": true,
  "total_time_secs": 18.164062584
}
```

The row passes outer-domain and Hessian checks but reports zero branch-and-bound evaluations and empty `bb_levels`. Walking through `ehp_general_ieee1788.rs:710-825`, the only path that produces this combination is when `create_initial_boxes_recursive` returns an empty vector and the for-loop short-circuits at line 736-739 with `proof_complete = true`. A wall-clock of `18.16 s` is consistent with the BB phase being effectively skipped and the Hessian + outer-domain checks running standalone.

In other words: the n=13 row's BB step did not exercise the search domain. The verdict `EHP_N13_PROVEN` was set without the same global-elimination evidence that backs n=12 (4.5×10⁷ evaluations) and n=14 (8.5×10⁸ evaluations).

The validator at `ehp114_proof_path_packets.rs:454-458` and 476-484 hardcodes this as a route exception specifically for `degree == 13`, with the comment "flag for route reconciliation before any refreshed public theorem packet." The validator is telling us, in source, that this row needs to be re-run before public submission.

## Re-run command (in-place, write to scripts/erdos-114)

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114

# Back up the current n=13 artifact before overwriting.
mkdir -p archive
cp EXP-MM-EHP-007-n13-inari_RESULTS.json \
   archive/EXP-MM-EHP-007-n13-inari_RESULTS_PRE_RERUN_20260507.json
cp EXP-MM-EHP-007-n13-inari_RESULTS.sha256 \
   archive/EXP-MM-EHP-007-n13-inari_RESULTS_PRE_RERUN_20260507.sha256

# Release build, single-degree run, in-place output to `.`
cargo build --release --bin ehp_general_ieee1788
./target/release/ehp_general_ieee1788 --outdir . 13 13

# Verify the new SHA-256 sidecar matches the new JSON
shasum -a 256 -c EXP-MM-EHP-007-n13-inari_RESULTS.sha256
```

Notes:

- The two `13 13` arguments force the explicit-list code path at `ehp_general_ieee1788.rs:1192-1195` (a single arg is interpreted as "start at 13, run through `MAX_PRECOMPUTED = 16`"). If you want only `n = 13`, pass it twice.
- The Cargo manifest pins `inari = "2.0"` and the `release` profile uses `opt-level = 3, lto = true`. Match this on the re-run host.
- Expected wall-clock: a properly-exercised BB run for `n = 13` should land between the `n = 12` (45M evals, on the order of minutes) and `n = 14` (855M evals, on the order of hours) budgets. A re-run that completes in under ~30 seconds and reports zero evals again means the box-construction step is the actual issue, not the BB engine.

## Acceptance gate for the new n=13 artifact

The re-run output is acceptable for inclusion in the public table if **all** of:

| Field | Required value |
|---|---|
| `verdict` | `EHP_N13_PROVEN` |
| `rigor` | `ieee_1788_interval_arithmetic_inari` |
| `bb_proof_complete` | `true` |
| `bb_total_evals` | greater than `0` (target: in the 10⁶–10⁹ range, sandwiched by `n=12` and `n=14`) |
| `bb_levels` length | greater than `0` |
| `outer_domain_safe` | `true` |
| `hessian_negative` | `true` |
| `l_star_lower` | within `[28.85923995, 28.85924000]` and not greater than `l_star_upper` |
| `l_star_upper` | within `[28.85923995, 28.85924000]` and not less than `l_star_lower` |

If any field fails, **do not** ship n=13 in v13 of the public packet. Quarantine the row, escalate to the BB-engine maintainer, and either ship v13 with `n != 13` (current v12 shape) or block v13.

If all fields pass:

1. Update the SHA-256 sidecar (`shasum -a 256 EXP-MM-EHP-007-n13-inari_RESULTS.json > EXP-MM-EHP-007-n13-inari_RESULTS.sha256`).
2. Promote `scripts/erdos-114/EXP-MM-EHP-007-n13-inari_RESULTS.{json,sha256}` to `results/erdos-114/` (the public path referenced from Zenodo and GitHub).
3. Cut a new packet version `EXP-MATH-EHP114-SMALL-N-SUBMISSION-READY-PACKET-20260507-13` that:
   - Updates the n=13 row hash in `_README.md`, `_RESULTS.json`, `_SHA_MANIFEST.md`, `_RESULTS.sha256`.
   - Restores n=13 to the public table (drops the "quarantined for v12" annotation).
   - Adds a `lemniscateLength_le_referenceLength_n13` axiom to `_LEAN_STUB.lean` between `n12` and `n14`.
   - Adjusts the finite theorem hypothesis from `(hp_deg_ne_13 : p.natDegree ≠ 13)` to drop that hypothesis.
   - Updates the Zenodo strategy memo: a v6 mint **is** appropriate now because a certificate row changed.
4. Mint Zenodo v6 with the replaced n=13 artifact. The DOI for the new version becomes the citation for the `n != 13` axioms in the Lean stub once the v13 packet is updated to point at v6.
5. Re-run prepub-redteam, Publisher Gate, Crackpot-Scrub on the v13 deposit set.
6. Update the formal-conjectures PR draft (or open it from v13).

## What this does not address

- The reason `create_initial_boxes_recursive` was returning empty for `n = 13` in the first place. If the underlying issue is a degree-specific bug in the box-construction code (e.g., a `coeff_bound` or `grid_per_axis` choice that produces zero boxes for `d = 23`), the re-run will surface it. Fix the bug before promoting any new artifact.
- The n=15 and n=16 rows (`EXP-MM-EHP-007-n15-inari_RESULTS.json`, `EXP-MM-EHP-007-n16-inari_RESULTS.json`) are out of scope for this packet. v12 explicitly ceilings at `n <= 14`.

## Tooling disclosure

This re-run document was prepared with AI-assisted tooling. The Rust code paths cited above were read out of the working tree at `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/`. The author takes responsibility for the diagnosis and the acceptance criteria.
