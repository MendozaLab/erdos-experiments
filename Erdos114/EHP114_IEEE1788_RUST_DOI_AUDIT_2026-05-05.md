# EHP114 IEEE-1788 / Rust / Zenodo DOI Audit

Date: 2026-05-05

Scope: read-only audit of the existing EHP #114 interval-certificate surface:
Zenodo DOI, local Rust/inari proof engine, canonical result JSONs, and current
claim ceiling for the Hessian/Puiseux-completion program.

## Bottom Line

The strongest hard artifact is not the May 5 Hessian-fantasy diagnostic. It is
the already published Zenodo record:

- DOI: `10.5281/zenodo.19480329`
- Concept DOI: `10.5281/zenodo.19184467`
- Title: `Computational Verification of the Erdős-Herzog-Piranian Conjecture for Degrees 3 <= n <= 14`
- Publication date: `2026-04-09`
- Version: `n3-n14`
- Resource type: `Dataset`
- Access: `open`

The Zenodo API reports 41 files, including every
`EXP-MM-EHP-007-n{3..14}-inari_RESULTS.json` result file, every matching
SHA-256 sidecar, `EHP_GENERAL_SUMMARY.json`, and
`Mendoza_EHP_n3-14_April2026.pdf`.

Local canonical files under
`erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n{3..14}-inari_RESULTS.*`
match the Zenodo file MD5s byte-for-byte for all result JSONs and sidecars.
The local SHA-256 sidecars also match the local JSON payloads for every
degree n = 3 through n = 14.

## What This Means

For book/agent/science guidance, the Rust/inari/DOI surface is the rigorous
anchor. The Hessian/Puiseux-completion language should point back to it as the
certified finite-n base case, not replace it.

The phrase to preserve is:

> shadow signature, not universal law.

The defensible claim is:

> We have a DOI-backed, byte-reconciled IEEE-1788/Rust-inari computational
> certificate for EHP degrees 3 <= n <= 14. The May 5 Hessian/Puiseux work
> supplies a candidate analytic compression of that certificate and a path
> toward a smaller closure theorem. It does not yet replace the interval
> certificate and does not prove the all-n problem.

## Current Certified Frontier

The local canonical result files report:

| n | Verdict | Reduced dim | B&B complete | Evals | Outer safe | Hessian negative |
|---:|---|---:|---|---:|---|---|
| 3 | `EHP_N3_PROVEN` | 3 | true | 7,560 | true | true |
| 4 | `EHP_N4_PROVEN` | 5 | true | 42,656 | true | true |
| 5 | `EHP_N5_PROVEN` | 7 | true | 312,598 | true | true |
| 6 | `EHP_N6_PROVEN` | 9 | true | 135,936 | true | true |
| 7 | `EHP_N7_PROVEN` | 11 | true | 23,552 | true | true |
| 8 | `EHP_N8_PROVEN` | 13 | true | 110,592 | true | true |
| 9 | `EHP_N9_PROVEN` | 15 | true | 507,904 | true | true |
| 10 | `EHP_N10_PROVEN` | 17 | true | 2,293,760 | true | true |
| 11 | `EHP_N11_PROVEN` | 19 | true | 10,223,616 | true | true |
| 12 | `EHP_N12_PROVEN` | 21 | true | 45,088,768 | true | true |
| 13 | `EHP_N13_PROVEN` | 23 | true | 0 | true | true |
| 14 | `EHP_N14_PROVEN` | 25 | true | 855,638,016 | true | true |

The material stress case is n = 14, not n = 20. The canonical n = 14 record
has real branch-and-bound work:

- Reduced dimension: 25
- Initial boxes: 33,554,432
- Interval evaluations: 855,638,016
- `L*` interval: `[30.852910841548532, 30.852910841548546]`
- Max non-extremizer upper bound: `9.792200323344273`
- B&B complete: true
- Outer-domain safe: true
- Hessian negative: true
- Total reported time: `269202.107988042` seconds

This is the serious certificate. It is the one a hostile reader has to answer.

## Rust Engine

The proof engine is:

`erdos-experiments/scripts/erdos-114/src/bin/ehp_general_ieee1788.rs`

Cargo manifest:

`erdos-experiments/scripts/erdos-114/Cargo.toml`

The binary name is:

`ehp_general_ieee1788`

Build check on this machine:

```text
cargo build --release --bin ehp_general_ieee1788
Finished `release` profile [optimized] target(s) in 5m 36s
```

The build completed with three dead-code warnings for the unused interval
complex helper path (`IC`, helper methods, and `eval_poly_interval`). Those
warnings do not invalidate the existing result artifacts, but they should be
cleaned before the next DOI/versioned release if the source is presented as a
polished artifact.

The source describes a four-step computer-assisted proof:

1. Tight interval enclosure for `L(z^n - 1)`.
2. Branch-and-bound over the symmetry-reduced coefficient space.
3. Outer-domain exclusion.
4. Interval Hessian negativity at the extremizer.

The code routes interval operations through the Rust `inari` crate and labels
the rigor surface as:

`ieee_1788_interval_arithmetic_inari`

The source currently has `L*` intervals precomputed through n = 16 and
`MAX_PRECOMPUTED = 16`.

## Important Caveat: n = 15 and n = 16

There are script-local files:

- `erdos-experiments/scripts/erdos-114/EXP-MM-EHP-007-n15-inari_RESULTS.json`
- `erdos-experiments/scripts/erdos-114/EXP-MM-EHP-007-n16-inari_RESULTS.json`

They report `EHP_N15_PROVEN` / `EHP_N16_PROVEN`, but with
`bb_total_evals = 0` and empty `bb_levels`. Those are not currently part of
the Zenodo n3-n14 DOI record and should not be promoted as certified frontier
without reconciliation.

Interpretation: n = 15 and n = 16 are promising precomputed-script artifacts,
not publication-safe certificates yet.

## What To Do Next

Do not restart from Python. Use the Rust/inari engine as the certificate
backend and the Hessian/Puiseux diagnostics as the theorem-target finder.

Recommended closure sequence:

1. Reconcile the n = 13 zero-eval local canonical record and the script-local
   n = 15/n = 16 zero-eval records. Confirm whether zero evaluations are a
   legitimate shortcut path or a serialization bug.
2. Freeze n = 14 as the DOI-backed interval-hardening calibration anchor.
3. Make a new versioned experiment before extending beyond the DOI record.
4. Add `L*` intervals for n = 17 through n = 20 only after the reproducible
   precompute path is documented.
5. Run feasibility probes before any full n = 20 B&B attempt; the full n = 20
   coefficient-space B&B should be assumed expensive until proven otherwise.
6. In parallel, harden the analytic compression route:
   radial interval bound first, shape cone second, mixed remainder last.

## Claim Ceiling

Safe:

- DOI-backed certified computational verification for n = 3 through n = 14.
- Rust/inari IEEE-1788 interval arithmetic certificate surface.
- Hessian/Puiseux diagnostics suggest an analytic compression of the existing
  certificate.

Unsafe:

- Saying #114 is completely closed.
- Saying the all-degree theorem is established by this packet.
- `Lean formal proof of EHP`
- Saying the n = 20 endpoint has an interval certificate.
- Promoting n = 15 or n = 16 until the zero-evaluation records are reconciled
  and, if appropriate, versioned and published.

## Commands Run

```bash
python3 - <<'PY'
# Summarized local canonical and script-local EXP-MM-EHP-007 result JSONs.
PY

python3 - <<'PY'
# Queried https://zenodo.org/api/records/19480329 and listed record metadata
# plus file checksums.
PY

python3 - <<'PY'
# Verified local SHA-256 sidecars match local JSON files for n = 3..14.
PY

python3 - <<'PY'
# Verified local MD5s match Zenodo MD5s for every n = 3..14 result JSON and
# SHA sidecar.
PY

cargo build --release --bin ehp_general_ieee1788
```

## Files Read

- `ERDOS_MASTER_SCORECARD.md`
- `README.md`
- `erdos-experiments/Erdos114/STRATIFICATION_N3_N14_2026-05-02.md`
- `erdos-experiments/Erdos114/KOOPMAN_LIFT_DESIGN_2026-05-02.md`
- `erdos-experiments/Erdos114/TENSOR_CONE_DESIGN_2026-05-02.md`
- `erdos-experiments/Erdos114/TOEPLITZ_MOMENT_LIFT_DESIGN_2026-05-02.md`
- `erdos-experiments/scripts/erdos-114/src/bin/ehp_general_ieee1788.rs`
- `erdos-experiments/scripts/erdos-114/Cargo.toml`
- `erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n{3..14}-inari_RESULTS.json`
- `erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n{3..14}-inari_RESULTS.sha256`
