# EHP114 n=14 Global Atlas Certificate Gate

Experiment: `EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-02`

## Theorem-Facing Obligation

L33 asks whether the n=14 root-affine atlas has theorem-grade local certificates for every one of its `64` cells. Each accepted cell packet must pass checksum, coverage, ownership, source-filter, candidate-count, and length-budget checks. The global certificate can pass only if all cells are present and accepted.

## Verdict

- Status: `GLOBAL_ATLAS_PASS_NOT_FULL_EHP_PROOF`
- Global atlas certificate pass: `true`
- Root-affine skeleton cell count: `64`
- Certified theorem-grade cells: `64`
- Missing theorem-grade cells: `0`
- Discovered local packet files: `66`
- Selected current packet files: `64`
- Superseded packet files ignored: `2`
- Rejected discovered cell packets: `0`
- Duplicate local cell certificates: `0`
- Source SHA fail count: `0`
- Coverage failure count: `0`
- Max certified cell length upper: `20.66460519288245`
- Exact cap: `20.672796062619668`
- First failed condition: `none`

## Interpretation

This run is an atlas-level integration certificate. It does not recompute the local geometry; it checks that the root-affine skeleton enumerates exactly 64 cells and that the selected latest local packet for each cell passes checksum, coverage, ownership, source-filter, candidate-drift, and length-budget checks. Older packet versions are preserved as immutable evidence but ignored as superseded current coverage.

## Next Dependency

If this gate passes, the next dependency is an independent theorem-packet audit over the atlas statement and then finite-degree/high-degree bridge accounting. Do not promote a full EHP114 claim from this local n=14 atlas gate alone.

## Claim Ceiling

L33 global atlas integration gate only. This run is not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. It may pass only after all 64 local theorem-grade cell packets are present and accepted.
