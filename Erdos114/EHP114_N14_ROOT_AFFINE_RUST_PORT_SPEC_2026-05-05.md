# EHP114 n=14 Root-Affine Rust Port Spec

Date: 2026-05-05

## Meaning

The Python low-dimensional certificate has found the right representation:
do interval arithmetic on affine root boxes, not on expanded coefficient boxes.
The next hardened engine should be Rust because the certificate layer needs
deterministic interval arithmetic, parallel subdivision, reproducible JSON
outputs, and reviewable performance.

This is still a shadow signature, not universal law. The Rust port is a
certificate-engine hardening step, not a proof of Erdős #114.

## Existing Rust Base

Use the existing crate:

```text
erdos-experiments/scripts/erdos-114/
```

It already has:

- `inari` interval arithmetic
- `rayon` parallelism
- `serde` / `serde_json` output
- `ehp114_batch_interval_lengths.rs`, the current pointwise interval
  lemniscate-length oracle
- `ehp114_n14_radial_compact_interval.rs`, the compact radial certificate
  style to imitate for claim ceilings and JSON artifacts

## New Binary Target

Add:

```text
src/bin/ehp114_n14_root_affine_cell.rs
```

and a `Cargo.toml` stanza:

```toml
[[bin]]
name = "ehp114_n14_root_affine_cell"
path = "src/bin/ehp114_n14_root_affine_cell.rs"
```

## Required Input

The binary should accept a JSON file with:

```json
{
  "experiment_id": "EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-REPRO-20260505-01",
  "degree": 14,
  "eps": 0.1,
  "extent": 3.0,
  "res": 220,
  "subdivision": 8,
  "lstar_lower": 30.852910841548532,
  "target": 10.180114778928864,
  "base_roots": [[1.0, 0.0]],
  "u0_direction": [[0.0, 0.0]],
  "u1_direction": [[0.0, 0.0]],
  "u0_interval": [-0.00175, 0.0],
  "u1_interval": [0.0, 0.00175]
}
```

The real input must contain all 14 complex roots/direction coordinates. The
stub above only shows shape.

## Required Computation

For each `8 x 8` subcell:

1. Build affine interval roots:

```text
r_i = base_i + eps^(1/28) * (a * u0_i + b * u1_i)
```

where `a` and `b` are subcell intervals.

2. Evaluate

```text
p(z) = prod_i (z - r_i)
```

with interval complex arithmetic.

3. Run the same interval marching-squares upper-bound functional as the Python
root-affine prototype.

4. Emit per-subcell rows:

```text
sub_i, sub_j, u0_interval, u1_interval,
length_upper, deficit_lower, margin_lower, pass,
active_cells, definite_case_cells, uncertain_corner_cells, ambiguous_avg_cells
```

## Reproduction Gate

Before generalizing, the Rust binary must reproduce the Python certificate:

```text
EXP-MATH-EHP114-N14-EPS01-ONE-CELL-ROOT-AFFINE-SUBDIV8-20260505-01
status = ONE_CELL_ROOT_AFFINE_SUBDIV8_PASS
failure_count = 0
max_length_upper <= 18.110795101362747 + tolerance
min_margin_lower >= 2.5620009612569206 - tolerance
```

Use a small tolerance only for harmless implementation-order differences in
floating-point-to-interval wrapping. The Rust output must remain conservative:
if it is tighter, explain why; if it is looser but still passes, record the new
margin.

## Output Contract

The Rust binary must write:

```text
*_RESULTS.json
*_REPORT.md
*_RESULTS.sha256
```

The report must state:

```text
This is a continuous grid-oracle cell certificate, not exact lemniscate
certification, not a Lean theorem, and not a proof of Erdős #114.
```

## Next Lift After Reproduction

Only after Rust reproduces the Python `SUBDIV8` PASS should it move to:

1. neighboring cells in the same `eps = 0.1` spectral-span disk;
2. all `164` candidate cells from the dense box-candidate run;
3. admissible-only epsilon slabs below `eps = 0.1`;
4. replacement of the marching-squares oracle functional by an exact or
   validated lemniscate-length enclosure.

The port is successful when it makes the certificate faster and more auditable
without raising the claim ceiling.

## Completion Note

The first Rust reproduction gate has been implemented and passed:

```text
binary = ehp114_n14_root_affine_cell
output = erdos-experiments/scripts/erdos-114/rust-root-affine-subdiv8-output-v2/
experiment_id = EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01
status = RUST_ROOT_AFFINE_SUBDIV8_REPRO_PASS
failure count = 0
rows = 64
max length upper = 18.110795101366655
min margin lower = 2.5620009612530126
```

Verified commands:

```bash
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114
cargo build --release --bin ehp114_n14_root_affine_cell
cd rust-root-affine-subdiv8-output-v2
shasum -a 256 -c EXP-MATH-EHP114-N14-EPS01-ROOT-AFFINE-RUST-SUBDIV8-REPRO-20260505-01_RESULTS.sha256
```

A clean `/tmp` rerun also reproduced the same PASS. This completes the
assigned Rust hardening gate for the selected cell. It does not complete exact
lemniscate certification.
