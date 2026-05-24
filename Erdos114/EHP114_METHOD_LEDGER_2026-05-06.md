# EHP114 Method Ledger For Orthodox Review

Date: 2026-05-06

Scope: local Erdős #114 / EHP n=14 hard-cell proof-route ledger, centered on subcell `(6,4)` and the validated-length cap `20.672796062619668`.

Claim ceiling: this ledger is not a proof of Erdős #114 and not a global n=14 certificate. It is a rigor ledger: what was tried, why it was mathematically reasonable, what the artifacts showed, why a route was retired or retained, and what remains live.

## Verification Summary

| Check | Result |
|---|---|
| Ledger rows | 47 |
| Retired routes | 9 |
| Retained diagnostic routes | 35 |
| Superseded planning routes | 2 |
| Active next route | 1 |
| Checksum verification | PASS on 2026-05-06 for cited lanes EHP114-L01 through EHP114-L14 and EHP114-L16 through EHP114-L33L where SHA sources exist |
| Overclaim scan | PASS on 2026-05-06 for repo Markdown, Downloads Markdown, sidecar notes, L31/L32/L33A-L33L artifacts, patched upstream Rust binaries, and existing L29/RH-MDL artifacts |
| JSON parse | PASS on 2026-05-06 |

## Reader Contract

The useful standard here is not whether every route succeeded. The useful standard is whether every route left an auditable residue:

- artifact path;
- checksum status when available;
- exact failure or retention reason;
- what the result rules out;
- what the result does not rule out;
- next proof-facing dependency.

Negative rows are therefore evidence of method control, not embarrassment. A route marked `retired` is retired as a proof-facing route, not erased as historical evidence.

## At-A-Glance Route Map

| Lane | Route | Decision | Meaning |
|---|---|---|---|
| EHP114-L00 | Calibration ladder n=3 through n=14 | retained_as_diagnostic | Establishes computational context only; not a hard-cell certificate. |
| EHP114-L01 | Radial/Puiseux/hypergeometric scaffold | retained_as_diagnostic | Supplies analytic shape of the boundary term; closure still needs interval hardening. |
| EHP114-L02 | Hessian/shape cone/local Taylor stability | retained_as_diagnostic | Retires naive uniform shape positivity; keeps radial reserve vs shape softening as the right local story. |
| EHP114-L03 | Low-dimensional cone probes | retained_as_diagnostic | Supports search geometry; does not certify the continuous cone. |
| EHP114-L04 | Exact-length budget extraction | retained_as_diagnostic | Establishes the cap and budget bookkeeping; not an exact-length certificate. |
| EHP114-L05 | Direct validated-length patches | retired | Regularity and ownership worked, but the method counted a thick tube. |
| EHP114-L06 | Branch census and slab validator | retained_as_diagnostic | Single-branch-per-slab assumption retired; branch-atlas framing retained. |
| EHP114-L07 | Targeted tube resolver | retired | Recursive exclusion/chart bookkeeping did not close the unresolved tubes. |
| EHP114-L08 | Bernstein tube certificate | retired | Bernstein subdivision reached depth limits while leaving too many pieces unresolved. |
| EHP114-L09 | One-dimensional Krawczyk tube certificate | retired | Krawczyk could not certify tubes under the full root-affine uncertainty. |
| EHP114-L10 | Parameter-sliced Krawczyk | retired | Parameter slicing was the cheapest repair; pilot still certified zero pieces. |
| EHP114-L11 | Bivariate Bernstein/Krawczyk | retired | Zero exclusions and zero certifications; failure attributed to dependency/coefficient swelling, not curve geometry. |
| EHP114-L12 | Third-party review | retained_as_diagnostic | Independent review agrees with retiring bivariate Bernstein/Krawczyk and pivoting to affine/Taylor or root-collar theorem. |
| EHP114-L13 | Root-collar Taylor/affine hard-cell pilot | retained_as_diagnostic | Dependency recovered on the first 10 pieces, but no strict collar candidate yet. |
| EHP114-L14 | Validated Taylor root-collar pilot | retained_as_diagnostic | Root-affine Taylor intervals keep derivative signs on 10/10 pieces, but Krawczyk still is not a strict subset. |
| EHP114-L15 | Analytic collar or higher-order Taylor model | superseded | Planning row superseded by the executed third-order pilot and the narrower active analytic theorem lane. |
| EHP114-L16 | Branch isolation collar-atlas pilot | retained_as_diagnostic | Branch-atlas wrapper is in place, but depth-2 wall sign separation did not certify the pilot sample. |
| EHP114-L17 | Normal-collar critical-point exclusion pilot | retired | Rotating to gradient-normal coordinates still certified zero regions; first-order collar geometry is no longer the next route. |
| EHP114-L18 | Third-order root-collar remainder pilot | retired | Third-order remainder control improved diagnosis but certified zero regions; sampled local collars are no longer the proof-facing route. |
| EHP114-L19 | Analytic critical-point/root-collar theorem target | superseded | Planning row superseded by the executed branch-centered moving-frame pilot. |
| EHP114-L20 | Branch-centered moving-frame collar pilot | retired | Attacked the exact C0 weakness but failed the Class A pass gate; local collar ladder is retired. |
| EHP114-L21 | Global critical-point exclusion target diagnostic | retained_as_diagnostic | Splits the 64 L18 residuals into 8 regular x-chart regions and 56 critical candidates. |
| EHP114-L22 | Regular residual decomposition diagnostic | retained_as_diagnostic | Fx stayed stable on all 8 regular regions, but raw x-wall sign separation failed on all 8. |
| EHP114-L23 | Regular-slice Taylor-y wall separation diagnostic | retained_as_diagnostic | Taylor-y reduced wall width but did not certify branch collars. |
| EHP114-L24 | Regular-slice monotone Taylor exclusion diagnostic | retained_as_diagnostic | Closed all 8 regular regions as monotone exclusions; no length promoted. |
| EHP114-L25 | Critical-candidate affine gradient diagnostic | retained_as_diagnostic | First-order affine gradient product closed 0 of 56 candidates. |
| EHP114-L26 | Direct p-prime exclusion diagnostic, p=4 | retained_as_diagnostic | Direct p-prime invariant closed 8 of 56 candidates. |
| EHP114-L27 | Direct p-prime exclusion diagnostic, p=8 | retained_as_diagnostic | p=8 closed the same 8; remaining obstruction is 48 derivative-root candidates. |
| EHP114-L28 | p-prime root-location diagnostic | retained_as_diagnostic | Root-free disk test closes 31 of 48 remaining candidates; 17 still need sharper p-prime root-location. |
| EHP114-L29 | sharp p-prime root-location diagnostic | retained_as_diagnostic | Sharp Taylor/Rouche local splitting closes all 17 derivative-root-near candidates; integration remains. |
| EHP114-L30 | RH/MDL Beurling-Nyman quantization-stability finite audit | retained_as_diagnostic | Checks 98 finite rows and 312 quantized entries; method-shaping only, not RH evidence. |
| EHP114-L31 | residual-chain integration certificate | retained_as_diagnostic | L24/L27/L28/L29 compose: 64/64 residual buckets closed, with no ownership/source/count drift. |
| EHP114-L32 | local hard-cell certificate packet | retained_as_diagnostic | Local hard-cell theorem packet passed: L31 plus child hashes are bound to the accepted length and local cap. |
| EHP114-L33 | global n=14 atlas certificate | retained_as_diagnostic | Gate ran cleanly and failed for the precise reason: 63 theorem-grade local cell packets are missing. |
| EHP114-L33A | missing-cell local packet generation gate | retained_as_diagnostic | Worklist emitted, but batch generation is blocked because seven proof-facing binaries are hard-coded to `(6,4)`. |
| EHP114-L33B | parameterized local-pipeline smoke | retained_as_diagnostic | Removed the hard-coded local-pipeline blocker; smoke now blocks honestly on missing per-cell source artifacts. |
| EHP114-L33C | per-cell source-generation smoke | retained_as_diagnostic | Source-chain smoke now passes after `CELL-00-00` source generation; it is a software gate only. |
| EHP114-L33D | CELL-00-00 slab-source generation | retained_as_diagnostic | First real `CELL-00-00` slab source exists, but current bound is over cap and still unresolved. |
| EHP114-L33E | CELL-00-00 source sharpening | retained_as_diagnostic | Overage is concentrated in top high-slope slabs; targeted slab repair is now justified. |
| EHP114-L33F | CELL-00-00 high-slope slab repair | retained_as_diagnostic | Targeted repair certifies all three worst high-slope branches and drops adjusted diagnostic total below cap. |
| EHP114-L33G | CELL-00-00 repaired source overlay | retained_as_diagnostic | Replacement accounting is explicit, under cap, and source-chain smoke passes. |
| EHP114-L33H | CELL-00-00 local certificate packet | retained_as_diagnostic | First non-hard cell packet passes: residual chain closes 64/64 with source hashes passing. |
| EHP114-L33I | CELL-00-01 controlled next-cell generation | retained_as_diagnostic | Second non-hard cell packet passes: high-slope source repair transfers and residual chain closes 64/64. |
| EHP114-L33J | representative-cell harness | retained_as_diagnostic | `CELL-00-02` passes, but interior `CELL-03-03` exposes variable regular-region closure as the next blocker. |
| EHP114-L33K | variable regular-region closure | retained_as_diagnostic | `CELL-03-03` now passes: 12/12 regular regions closed, L31 closes 64/64, and local packet is under cap. |
| EHP114-L33L | representative-cell harness continuation | retained_as_diagnostic | `CELL-07-07` and `CELL-07-00` both pass after single high-slope source repairs; residual chains close 64/64. |
| EHP114-L33M | controlled remaining-cell batch generation | active_next | Start the remaining n=14 cell batch under source/subcell contracts, stopping on the first unrepairable blocker. |

## Detailed Register

### EHP114-L00 - Calibration Ladder n=3 Through n=14

- Date: 2026-05-02 and earlier local calibration runs.
- Method or route: EHP calibration ladder, including n=3 proof-of-concept runs and `EXP-MM-EHP-007-n3-inari` through `EXP-MM-EHP-007-n14-inari`.
- Artifact IDs: `EHP_N3_*`, `EXP-MM-EHP-007-n*-inari`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EHP_N3_LEVEL2_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EHP_N3_LEVEL3_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MM-EHP-007-n14-inari_RESULTS.json`
- Status: calibration artifacts present.
- What we tried: establish a computational ladder from toy and small-n EHP runs toward the n=14 hard-cell regime.
- Why it was reasonable: n=14 should not be attacked before the code path has lower-n sanity checks and reproducible artifact conventions.
- Observed result: the ladder exists and provides context, but it does not isolate or certify the hard n=14 branch pieces.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating the n=14 sprint as an ungrounded first run.
- What it does not rule out: it does not prove the n=14 hard-cell certificate or the global EHP statement.
- Next dependency: keep lower-n artifacts as calibration, but do proof-facing work on the n=14 hard-cell certificate.
- Claim ceiling: calibration only.
- Orthodox reader note: lower-n sanity checks are necessary hygiene, not theorem evidence for n=14.
- Verification commands: checksum files exist for the `EXP-MM-EHP-007-*` family; targeted checks should be run when citing any specific lower-n result.
- SHA status: not fully enumerated in this ledger version.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L01 - Radial/Puiseux/Hypergeometric Scaffold

- Date: 2026-05-02 to 2026-05-05.
- Method or route: radial hypergeometric calibration and Puiseux interval target.
- Artifact IDs:
  - `EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01`
  - `EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01`
  - `EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01`
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/erdos-114/EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01_RESULTS.json`
- Status: scaffold and target artifacts present.
- What we tried: isolate the radial boundary term and formulate the Puiseux/hypergeometric singularity model as the analytic reserve.
- Why it was reasonable: the hard-cell behavior appears boundary-dominated; a radial singular model is the natural analytic object before tangential shape modes are mixed in.
- Observed result: the scaffold produced a target and claim ceiling, but not a closed interval-hard theorem.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating the radial boundary term as merely metaphorical; there is an explicit target artifact.
- What it does not rule out: it does not close the radial theorem or certify the mixed local cone.
- Next dependency: interval-harden the radial/Puiseux singularity model and bind it to local shape remainders.
- Claim ceiling: shadow signature, not universal law; certificate target, not proof.
- Orthodox reader note: this is admissible as theorem planning only if it remains below the proof-claim line.
- Verification commands:
  - `shasum -a 256 -c EXP-MATH-EHP114-RADIAL-HYPERGEOMETRIC-20260502-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-RADIAL-HYPERGEOMETRIC-CALIBRATION-20260502-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-RADIAL-PUISEUX-INTERVAL-TARGET-20260505-01_RESULTS.sha256`
- SHA status: PASS on 2026-05-06 for the three listed checksum files.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L02 - Hessian/Shape Cone/Local Taylor Stability

- Date: 2026-05-05.
- Method or route: Hessian fantasy, interval Taylor M14 packet, spectral diagnosis, and shape interval matrix.
- Artifact IDs:
  - `EXP-MATH-EHP114-HESSIAN-FANTASY-20260505-01`
  - `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01`
  - `EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01`
  - `EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01`
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-HESSIAN-FANTASY-20260505-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01_RESULTS.json`
- Status: `INTERVAL_TAYLOR_MATRIX_FAIL`, `RADIAL_BASE_SHAPE_SOFTENING_DETECTED`, and `SHAPE_INTERVAL_MATRIX_CERTIFIED` across the cited diagnostics.
- What we tried: test whether local shape modes preserve uniform positivity after radial contraction and whether the finite-difference shape matrix could be certified.
- Why it was reasonable: a quadratic/Hessian-style local stability theorem would be the cleanest closure if it were true.
- Observed result: the naive uniform positivity theorem failed, while the shape matrix and spectral diagnostics identified radial-base shape softening as the real local phenomenon.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the simple theorem "the transported local shape cone is uniformly positive" as currently framed.
- What it does not rule out: it does not rule out a mixed theorem where radial reserve absorbs negative shape curvature on a smaller or scaled cone.
- Next dependency: formulate and certify a mixed radial/shape remainder theorem.
- Claim ceiling: diagnostic for interval Taylor packet; not proof or disproof of EHP #114.
- Orthodox reader note: this is exactly the kind of failed-simple-theorem record that prevents overfitting a pretty story.
- Verification commands:
  - `shasum -a 256 -c EXP-MATH-EHP114-HESSIAN-FANTASY-20260505-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-M14-PACKET-20260505-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-INTERVAL-TAYLOR-SPECTRAL-DIAGNOSIS-20260505-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-SHAPE-INTERVAL-MATRIX-20260505-01_RESULTS.sha256`
- SHA status: PASS on 2026-05-06 for the four listed checksum files.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L03 - Low-Dimensional Cone Probes

- Date: 2026-05-05.
- Method or route: low-dimensional cone, coefficient-box, root-affine, and variation probes.
- Artifact ID: `EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/low_dim_cone/EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_RESULTS.json`
- Status: `LOW_DIM_CONE_GRID_PASS`.
- What we tried: reduce the hard local cone to lower-dimensional search boxes and look for scalar-theorem evidence.
- Why it was reasonable: low-dimensional probes are cheap falsification tests before committing to high-dimensional certificates.
- Observed result: sampled grid evidence was favorable, but the artifact itself states it is not a continuous cone certificate.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out some obvious low-dimensional counter-behavior on the sampled grid.
- What it does not rule out: it does not certify the continuous cone, admissibility, or scalar deficit lower bound.
- Next dependency: promote sampled grid evidence to interval boxes bound to admissibility and scalar-deficit lower bounds.
- Claim ceiling: finite 2D coefficient-disk grid evidence only.
- Orthodox reader note: a sampled pass is useful for choosing a theorem target, not for asserting the theorem.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-LOW-DIM-CONE-BOX-PROBE-20260505-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L04 - Exact-Length Budget Extraction

- Date: 2026-05-05.
- Method or route: budget extraction from marching-squares certificate.
- Artifact ID: `EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/exact_length_lift/EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.json`
- Status: `BUDGET_ONLY_NOT_EXACT_LENGTH_CERTIFICATE`.
- What we tried: extract and formalize the local hard-cell length cap and budget bookkeeping.
- Why it was reasonable: every later certificate needs a fixed numerical cap and margin semantics.
- Observed result: the artifact produced a budget framing but explicitly did not compute exact lemniscate length.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out moving forward without a named cap and budget target.
- What it does not rule out: it does not validate marching squares against exact length.
- Next dependency: replace budget extraction with direct validated implicit-curve length.
- Claim ceiling: budget bookkeeping only.
- Orthodox reader note: this is an accounting artifact, not a geometric certificate.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-EXACT-LENGTH-LIFT-BUDGET-20260505-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L05 - Direct Validated-Length Patches

- Date: 2026-05-05.
- Method or route: direct patchwise implicit-curve length enclosure.
- Artifact ID: `EXP-MATH-EHP114-N14-VALIDATED-LENGTH-HARD-CELL-20260505-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-VALIDATED-LENGTH-HARD-CELL-20260505-01/EXP-MATH-EHP114-N14-VALIDATED-LENGTH-HARD-CELL-20260505-01_RESULTS.json`
- Status: `VALIDATED_LENGTH_FAIL_BUDGET`.
- What we tried: certify active implicit-curve boxes directly using patchwise length upper bounds.
- Why it was reasonable: a direct implicit-curve length certificate is closer to the desired theorem than a marching-squares bridge.
- Observed result: `unresolved_box_count = 0` and `excluded_box_count = 12340542`, but `total_validated_length_upper = 119.57917312364832`, far above the hard-cell cap. The method counted a thick interval tube rather than one owned branch.
- Discard decision: `retired`.
- What it rules out: it rules out thick-tube patch summation as the proof anchor.
- What it does not rule out: it does not rule out the true curve length being below cap.
- Next dependency: isolate owned branches before summing length.
- Claim ceiling: local diagnostic failure only.
- Orthodox reader note: the failure is structural and informative; it prevented a false certificate from being promoted.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-VALIDATED-LENGTH-HARD-CELL-20260505-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L06 - Branch Census And Slab Validator

- Date: 2026-05-06.
- Method or route: branch-count census and slab root-isolation validated length.
- Artifact IDs:
  - `EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01`
  - `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03`
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03_RESULTS.json`
- Status: `BRANCH_CENSUS_FAIL_AMBIGUOUS_ROOT_ISOLATION` and `SLAB_FAIL_ROOT_ISOLATION`.
- What we tried: replace thick-tube counting with slab-based branch isolation and branch-length summation.
- Why it was reasonable: once thick-tube summation failed, the natural repair was one graph segment per isolated branch per slab.
- Observed result: the z32 slab artifact achieved local accepted length `20.316451752723314`, below cap `20.672796062619668`, but left `4398` unresolved branch tubes.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the simple one-branch-per-slab assumption.
- What it does not rule out: it does not rule out a branch-atlas certificate; unresolved tubes remain the blocker.
- Next dependency: resolve remaining tubes by stronger isolation/certification, not by counting all active boxes.
- Claim ceiling: local hard-cell diagnostic; not a global proof.
- Orthodox reader note: this is the first route that separated length budget from isolation budget.
- Verification commands:
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01_RESULTS.sha256`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03_RESULTS.sha256`
- SHA status: PASS on 2026-05-06 for both listed checksum files.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L07 - Targeted Tube Resolver

- Date: 2026-05-06.
- Method or route: recursive targeted resolver over unresolved branch tubes.
- Artifact ID: `EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01_RESULTS.json`
- Status: `TUBE_RESOLVER_FAIL_DEPTH_LIMIT`.
- What we tried: recursively subdivide unresolved tubes, first excluding empty pieces, then certifying chart directions.
- Why it was reasonable: it attacked only the unresolved set from the slab artifact, preserving the accepted length baseline.
- Observed result: accepted length stayed `20.316451752723314` with margin `0.3563443098963539`, but `73080` unresolved pieces remained after the resolver limits; `93132` pieces were excluded.
- Discard decision: `retired`.
- What it rules out: it rules out more interval-subdivision bookkeeping as sufficient in this representation.
- What it does not rule out: it does not rule out stronger analytic or dependency-preserving certificates on the same pieces.
- Next dependency: use dependency-aware arithmetic or a root-collar theorem.
- Claim ceiling: local diagnostic failure only.
- Orthodox reader note: this route failed after narrowing the target set, so the failure is not from global overbreadth.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-TUBE-RESOLVER-HARD-CELL-20260506-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L08 - Bernstein Tube Certificate

- Date: 2026-05-06.
- Method or route: Bernstein-form tube certificate on remaining pieces.
- Artifact ID: `EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.json`
- Status: `BERNSTEIN_TUBE_FAIL_DEPTH_LIMIT`.
- What we tried: use Bernstein-form enclosure to exclude/certify tube pieces.
- Why it was reasonable: Bernstein coefficients are a standard way to get polynomial range enclosures and exclusion certificates.
- Observed result: accepted length stayed under cap, `87216` pieces were excluded, but `165782` pieces remained unresolved and no Krawczyk-style piece was certified.
- Discard decision: `retired`.
- What it rules out: it rules out this Bernstein tube certificate as the next proof-facing closure route.
- What it does not rule out: it does not rule out Taylor-model or affine arithmetic, which preserve dependency differently.
- Next dependency: move away from plain Bernstein enclosure on these pieces.
- Claim ceiling: local certificate route failure only.
- Orthodox reader note: this is a normal validated-numerics outcome: a basis can be rigorous but too wide to certify.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-BERNSTEIN-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L09 - One-Dimensional Krawczyk Tube Certificate

- Date: 2026-05-06.
- Method or route: one-dimensional Krawczyk chart certification on unresolved tube pieces.
- Artifact ID: `EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.json`
- Status: `KRAWCZYK_TUBE_FAIL_DEPTH_LIMIT`.
- What we tried: certify isolated chart roots with Krawczyk-style inclusion while carrying the root-affine uncertainty.
- Why it was reasonable: the implicit function theorem suggests local one-dimensional charts should be certifiable if derivative signs and contraction are tight enough.
- Observed result: accepted length remained below cap, `60796` pieces were excluded, but `270768` pieces remained unresolved and no Krawczyk-certified pieces were produced.
- Discard decision: `retired`.
- What it rules out: it rules out one-dimensional Krawczyk under the whole root-affine cell as sufficient.
- What it does not rule out: it does not rule out Krawczyk after a better dependency representation or an analytic collar.
- Next dependency: reduce parameter dependency or switch arithmetic model.
- Claim ceiling: local diagnostic failure only.
- Orthodox reader note: the failure identifies parameter dependency as load-bearing, not length budget.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-KRAWCZYK-TUBE-CERT-HARD-CELL-20260506-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L10 - Parameter-Sliced Krawczyk

- Date: 2026-05-06.
- Method or route: root-parameter-sliced Krawczyk resolver.
- Artifact ID: `EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01`.
- Artifact path:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01_RESULTS.json`
- Status: `PARAM_SLICED_KRAWCZYK_PILOT_FAILS`.
- What we tried: split the root-affine parameter cell and rerun Krawczyk on tighter parameter tiles.
- Why it was reasonable: if full-cell dependency was the issue, parameter slicing was the cheapest targeted repair.
- Observed result: the pilot still certified zero pieces.
- Discard decision: `retired`.
- What it rules out: it rules out cheap parameter slicing as the missing ingredient.
- What it does not rule out: it does not rule out full affine/Taylor models or a hand-proved analytic collar.
- Next dependency: stop adding slices; switch representation.
- Claim ceiling: pilot failure only.
- Orthodox reader note: this is a clean example of retiring the cheap repair before investing in heavier machinery.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-PARAM-SLICED-KRAWCZYK-PILOT-HARD-CELL-20260506-01_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L11 - Bivariate Bernstein/Krawczyk

- Date: 2026-05-06.
- Method or route: full bivariate Bernstein/Krawczyk pilot.
- Artifact IDs:
  - `EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-01`
  - `EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02`
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02/EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02_RESULTS.json`
- Status: `BIVARIATE_BERNSTEIN_KRAWCZYK_PILOT_FAILS`.
- What we tried: use bivariate Bernstein control intervals for exclusion, derivative sign/lower-bound certification, chart eligibility, and Krawczyk-style inclusion on `4096` unresolved pieces.
- Why it was reasonable: this was the natural stronger polynomial-certificate route after one-dimensional chart methods failed.
- Observed result: `processed_piece_count = 4096`, `bivariate_excluded_piece_count = 0`, `krawczyk_certified_piece_count = 0`, `remaining_unresolved_piece_count = 4096`, with accepted length still below cap at `20.316451752723314`. The later review records ordinary interval signs as tight while Bernstein hulls are extremely wide, indicating dependency/coefficient swelling.
- Discard decision: `retired`; `20260506-01` is superseded by schema-correct `20260506-02`.
- What it rules out: it rules out bivariate Bernstein/Krawczyk in this coordinate representation as the next closure route.
- What it does not rule out: it does not rule out the curve geometry, the accepted length budget, affine/Taylor arithmetic, or an analytic root-collar theorem.
- Next dependency: prototype affine/Taylor model arithmetic or prove a root-collar theorem on the hardest unresolved pieces.
- Claim ceiling: local method retirement only.
- Orthodox reader note: the correct conclusion is narrow: representation failure, not theorem failure.
- Verification command: `shasum -a 256 -c EXP-MATH-EHP114-N14-BIVARIATE-BERNSTEIN-KRAWCZYK-PILOT-HARD-CELL-20260506-02_RESULTS.sha256`.
- SHA status: PASS on 2026-05-06 for `20260506-02`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L12 - Third-Party Review

- Date: 2026-05-06.
- Method or route: external/model methodological review of bivariate Bernstein/Krawczyk failure.
- Artifact IDs:
  - `Analysis & Opinion_ EHP114 Bivariate Bernstein_Kra.md`
  - `THIRD_PARTY_REVIEW_EHP114_BIVARIATE_BERNSTEIN_KRAWCZYK_FAILURE_2026-05-06.md`
- Artifact paths:
  - `/Users/kenbengoetxea/Downloads/Analysis & Opinion_ EHP114 Bivariate Bernstein_Kra.md`
  - `/Users/kenbengoetxea/Downloads/THIRD_PARTY_REVIEW_EHP114_BIVARIATE_BERNSTEIN_KRAWCZYK_FAILURE_2026-05-06.md`
- Status: review materialized locally.
- What we tried: ask whether retiring bivariate Bernstein/Krawczyk is methodologically justified.
- Why it was reasonable: before abandoning a standard validated-numerics route, the retirement decision should be challenged externally.
- Observed result: the review agrees that the route should be retired and identifies dependency/coefficient swelling as the likely cause; it recommends affine/Taylor model arithmetic or an analytic root-collar theorem.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out pretending the bivariate Bernstein/Krawczyk failure is just a routine tuning problem.
- What it does not rule out: it does not constitute a literature survey or proof certificate.
- Next dependency: implement root-collar/affine/Taylor prototype with predeclared gates.
- Claim ceiling: reasoning check only.
- Orthodox reader note: the review is useful as methodological corroboration, but the local artifacts remain primary.
- Verification commands: no checksum for Downloads review documents; local file presence verified on 2026-05-06.
- SHA status: `NO_SHA_SOURCE`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L13 - Root-Collar Taylor/Affine Hard-Cell Pilot

- Date: 2026-05-06.
- Method or route: root-collar theorem target plus affine/Taylor model arithmetic on the hardest unresolved pieces.
- Artifact ID: `EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01_REPORT.md`
- Status: `ROOT_COLLAR_TAYLOR_AFFINE_PILOT_DEPENDENCY_RECOVERED_NO_COLLAR`.
- What we tried: preserve first-order dependency through a centered Taylor/affine local model on the first 10 unresolved pieces from the bivariate Bernstein/Krawczyk failure artifact.
- Why it is reasonable: ordinary interval diagnostics saw tight derivative signs where Bernstein hulls widened dramatically; affine/Taylor methods are designed for this failure mode.
- Observed result: processed 10 pieces; source derivative signs were stable on 10/10; Taylor derivative signs were stable on 10/10; Taylor intervals were materially tighter than Bernstein on 10/10; strict root-collar candidates were 0/10.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating the Bernstein failure as pure geometric difficulty; dependency-aware Taylor coordinates recover the local derivative structure.
- What it does not rule out: it does not produce a proof-grade collar, does not close unresolved pieces, and does not include full root-affine Taylor validation.
- Next dependency: promote the centered Taylor/affine diagnostic into a validated Taylor model over root-affine uncertainty, or prove a hand analytic root-collar lemma.
- Claim ceiling: local root-collar Taylor/affine pilot only.
- Orthodox reader note: this is a useful split result: representation improves, but certification still fails.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_root_collar_taylor_affine_pilot`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-ROOT-COLLAR-TAYLOR-AFFINE-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L14 - Validated Taylor Root-Collar Pilot

- Date: 2026-05-06.
- Method or route: validated Taylor model over root-affine uncertainty, with interval Hessian bounds and interval Krawczyk collar test.
- Artifact ID: `EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01_REPORT.md`
- Status: `VALIDATED_TAYLOR_ROOT_COLLAR_DERIVATIVE_SIGNED_NO_COLLAR`.
- What we tried: convert L13 dependency recovery into a rigorous enclosure that carries root-affine uncertainty and tests a strict interval Newton/Krawczyk collar on the first 10 worst pieces.
- Why it is reasonable: L13 showed the Taylor/affine coordinate choice fixes the Bernstein dependency blow-up on the first 10 pieces, but the midpoint diagnostic is not proof-grade.
- Observed result: processed 10 pieces; validated Taylor derivative signs remained stable on 10/10; strict root-collar candidates were 0/10. The first piece had validated `Fx` interval approximately `[-52.21, -50.22]` and `Fy` interval approximately `[-8.06, -6.71]`, but Krawczyk still was not a strict subset.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the idea that root-affine uncertainty immediately destroys derivative sign stability on the sampled hard pieces.
- What it does not rule out: it does not certify any branch collar, close unresolved pieces, or prove the local hard-cell length certificate.
- Next dependency: sharpen to an analytic root-collar lemma or move to a higher-order Taylor model with explicit third-derivative remainder.
- Claim ceiling: local validated Taylor root-collar pilot only.
- Orthodox reader note: this is stronger than L13 because root uncertainty is carried explicitly; it still fails exactly at the strict-collar step.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_validated_taylor_root_collar`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-VALIDATED-TAYLOR-ROOT-COLLAR-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L15 - Superseded Planning Route: Analytic Collar Or Higher-Order Taylor Model

- Date: proposed 2026-05-06.
- Method or route: sharper analytic root-collar lemma, or higher-order Taylor model with explicit third-derivative remainder.
- Proposed artifact ID: `EXP-MATH-EHP114-N14-ANALYTIC-ROOT-COLLAR-LEMMA-TARGET-20260506-01`.
- Artifact paths: not yet emitted.
- Status: `superseded_by_L18_third_order_pilot`.
- What we tried: planning lane for either a sharper analytic root-collar lemma or a higher-order Taylor model with explicit third-derivative remainder.
- Why it is reasonable: L14 shows derivative sign stability survives root-affine uncertainty; the missing piece is contraction/sharpness, not sign regularity.
- Observed result: superseded on 2026-05-06 by L18, which implemented the higher-order Taylor remainder branch of this planned route.
- Discard decision: `superseded`.
- What it rules out: nothing directly; this was a planning lane and is now replaced by the executed L18 artifact plus the narrower active analytic theorem lane.
- What it does not rule out: it does not rule out analytic critical-point exclusion, a sharper root-collar theorem, or future affine/Taylor model arithmetic.
- Next dependency: see EHP114-L18 for the executed third-order pilot and EHP114-L19 for the active analytic critical-point theorem target.
- Claim ceiling: planned local hard-cell theorem-target route only.
- Orthodox reader note: planning rows should not remain active after the corresponding diagnostic has run; the active lane is now the narrower analytic theorem target.
- Verification commands: superseded by L18 artifact verification.
- SHA status: `NOT_YET_RUN`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L16 - Branch Isolation Collar-Atlas Pilot

- Date: 2026-05-06.
- Method or route: branch-atlas pilot using interval exclusion plus monotone collar wall certification.
- Artifact ID: `EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01_REPORT.md`
- Status: `BRANCH_ISOLATION_FAIL_COLLAR`.
- What we tried: consume the z32 slab unresolved branches, preserve source ownership keys, try interval exclusion first, then certify monotone graph collars from fixed nonzero derivative plus opposite signed collar walls on a 16-branch stratified pilot.
- Why it was reasonable: the z32 slab artifact had accepted length under cap and zero duplicate ownership, but 4398 unresolved branch tubes. A branch atlas is the right layer between topology and length.
- Observed result: processed 16 source unresolved branches. The recursion produced 64 unresolved subregions, 0 certified branches, 0 excluded regions, and 0 ownership duplicates. The first failed condition was `collar_wall_sign_separation_failed`; 56 regions also hit `critical_collar_derivatives_not_sign_stable` at depth 2.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out first-order monotone wall-collar certification at depth 2 as an immediately closing branch-isolation route on the pilot sample.
- What it does not rule out: it does not rule out branch-atlas certification, rotated charts, analytic critical-point exclusion, higher-order Taylor collars, or the favorable local length budget.
- Next dependency: move to a sharper analytic root-collar theorem, rotated-chart critical-point exclusion, or higher-order Taylor model with explicit remainder before increasing raw subdivision depth.
- Claim ceiling: local branch-isolation collar-atlas pilot only.
- Orthodox reader note: topology ownership bookkeeping is now explicit; the hard inequality is wall sign separation under root-affine uncertainty.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_branch_isolation_collar_atlas`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-BRANCH-ISOLATION-COLLAR-ATLAS-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L17 - Normal-Collar Critical-Point Exclusion Pilot

- Date: 2026-05-06.
- Method or route: normal/tangent collar pilot using midpoint-gradient coordinates.
- Artifact ID: `EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01_REPORT.md`
- Status: `NORMAL_COLLAR_FAIL_CRITICAL_POINT_EXCLUSION`.
- What we tried: rotate each of the 64 unresolved branch-atlas subregions into midpoint-gradient normal/tangent coordinates, then test whether fixed nonzero normal derivative plus opposite signed normal walls certifies a graph branch over the tangent direction.
- Why it was reasonable: the coordinate-axis collar failed partly because many boxes lost a stable x/y chart. A gradient-normal chart is the natural coordinate-free repair before abandoning first-order collar geometry.
- Observed result: processed 64 regions; normal-collar certified 0, excluded 0, remaining unresolved 64. There were 56 critical-point-exclusion failures and 8 wall-separation failures. Length remained under cap because no new branch length was promoted.
- Discard decision: `retired`.
- What it rules out: it rules out first-order normal/tangent collar geometry as an immediately useful proof-facing closure route on the current branch-atlas pilot residuals.
- What it does not rule out: it does not rule out analytic critical-point exclusion, higher-order Taylor models with explicit remainder, smaller theorem-shaped collars, or the favorable local length budget.
- Next dependency: move to analytic critical-point exclusion or explicit third-derivative Taylor remainder control; do not keep escalating first-order collar geometry.
- Claim ceiling: local normal-collar critical-point exclusion pilot only.
- Orthodox reader note: the failed normal chart is valuable because it shows the obstacle is not merely a poor axis choice; the next inequality must control critical-point exclusion or higher-order remainders.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_normal_collar_critical_exclusion_pilot`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-NORMAL-COLLAR-CRITICAL-EXCLUSION-PILOT-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L18 - Third-Order Root-Collar Remainder Pilot

- Date: 2026-05-06.
- Method or route: third-order Taylor collar remainder pilot in midpoint-gradient normal/tangent coordinates.
- Artifact ID: `EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_REPORT.md`
- Status: `THIRD_ORDER_FAIL_CRITICAL_EXCLUSION`.
- What we tried: consume the 64 normal-collar unresolved regions, compute `p`, `p'`, `p''`, and `p'''` by recurrence, and test explicit Taylor inequalities for critical-point exclusion and normal-wall separation.
- Why it was reasonable: L17 showed first-order normal collars were too weak; third-order remainder bounds are the cheapest theorem-shaped upgrade before abandoning sampled local collars entirely.
- Observed result: processed 64 regions. Third-order certified 0, excluded 0, and left 64 unresolved. Critical exclusion passed on 19 regions but wall separation passed on 0. The first failing inequality was `sup_r |F(0,r)| + normal_remainder = 0.8371883316743337` not below `S * lower(|F_n|) = 0.03309058624726584`. The total length stayed at `20.316451752723314` under the cap `20.672796062619668` because no new branch length was promoted.
- Discard decision: `retired`.
- What it rules out: it rules out the current sampled third-order local collar remainder inequality as a closure route for the normal-collar residual set.
- What it does not rule out: it does not rule out a sharper analytic critical-point theorem, a different collar scale, a root-collar theorem with better constants, or the favorable local length budget.
- Next dependency: write the analytic critical-point/root-collar theorem target with constants for wall separation and derivative lower bounds before running more subdivision bookkeeping.
- Claim ceiling: local third-order root-collar remainder pilot only.
- Orthodox reader note: the failure is specific: higher-order remainder control improved diagnosis but did not bridge center-strip wall separation. It is not evidence against the mathematical route, only against this sampled inequality at this scale.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_third_order_collar_remainder_pilot`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L19 - Superseded Planning Route: Analytic Critical-Point/Root-Collar Theorem Target

- Date: proposed 2026-05-06.
- Method or route: analytic critical-point/root-collar theorem target.
- Proposed artifact ID: `EXP-MATH-EHP114-N14-ANALYTIC-CRITICAL-POINT-COLLAR-THEOREM-TARGET-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_ANALYTIC_CRITICAL_POINT_COLLAR_THEOREM_TARGET_2026-05-06.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md`
- Status: `superseded_by_L20_moving_frame_pilot`.
- What we tried: wrote the theorem-target lane document for the next analytic critical-point/root-collar route; no proof artifact has run yet.
- Why it is reasonable: L16 through L18 show that raw local collars fail at wall separation or critical exclusion even when length budget is favorable. The next proof-facing object must be a theorem-shaped critical-point and wall-separation inequality, not another bookkeeping refinement.
- Observed result: superseded by L20, which implemented the branch-centered moving-frame collar pilot on the Class A wall failures.
- Discard decision: `superseded`.
- What it rules out: nothing directly; this planning row is replaced by the executed L20 artifact.
- What it does not rule out: it does not rule out global analytic critical-point exclusion or a different residual-domain decomposition.
- Next dependency: see EHP114-L20 for the executed pilot and EHP114-L21 for the active global analytic critical-point exclusion target.
- Claim ceiling: planned local hard-cell theorem-target route only.
- Orthodox reader note: planning rows should not remain active after execution; the active lane is now the global critical-point theorem target.
- Verification commands: overclaim scan over `EHP114_ANALYTIC_CRITICAL_POINT_COLLAR_THEOREM_TARGET_2026-05-06.md` and `EHP114_CENTER_STRIP_ROOT_COLLAR_REDUCTION_PACKET_2026-05-06.md`.
- SHA status: `NO_SHA_SOURCE`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L20 - Branch-Centered Moving-Frame Collar Pilot

- Date: 2026-05-06.
- Method or route: branch-centered moving-frame collar pilot on Class A wall failures.
- Artifact ID: `EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01_REPORT.md`
- Status: `MOVING_FRAME_COLLAR_PILOT_FAIL_BRANCH_CENTERING`.
- What we tried: processed the 14 Class A wall failures from L18, prioritized key `2286:0`, and tested whether branch centering on the normal axis plus a moving-frame tangent-zero hypothesis makes the quadratic center-strip bound fit the wall budget.
- Why it was reasonable: L18 showed the first wall failures were dominated by raw `C0` rather than normal remainder. A branch-centered moving frame was the cleanest way to kill the tangent-linear center-strip term before escalating to a global theorem.
- Observed result: processed 14 Class A candidates. Validated center points 1/14, moving-frame tangent-zero count 1/14, quadratic center-strip passes 0/14, remaining Class A unresolved 14/14. The easiest family `2286:0` failed branch centering on the normal axis. Length remained `20.316451752723314` under cap `20.672796062619668`; no branch length was promoted.
- Discard decision: `retired`.
- What it rules out: it rules out the current branch-centered moving-frame local-collar pilot as a proof-facing closure route on the Class A residuals, including the easiest `2286:0` family.
- What it does not rule out: it does not rule out a global analytic critical-point exclusion theorem, a different residual-domain decomposition, or the favorable local length budget.
- Next dependency: move to a global analytic critical-point exclusion target for the 64 L18 residual regions; do not keep iterating local normal-axis collar variants.
- Claim ceiling: local branch-centered moving-frame collar pilot only.
- Orthodox reader note: this is a strong negative control: even after attacking the exact `C0` weakness, the inherited local collar geometry failed at branch centering or quadratic budget. The next step must change the theorem, not the sampling density.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_branch_centered_moving_frame_collar_pilot`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-BRANCH-CENTERED-MOVING-FRAME-COLLAR-PILOT-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L21 - Global Critical-Point Exclusion Target Diagnostic

- Date: 2026-05-06.
- Method or route: global critical-point exclusion target diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_GLOBAL_ANALYTIC_CRITICAL_POINT_EXCLUSION_TARGET_2026-05-06.md`
- Status: `GLOBAL_CRITICAL_POINT_TARGET_REGULAR_REGIONS_FOUND`.
- What we tried: built and ran a Rust diagnostic over the 64 L18 residual regions, classifying each region as excluded, regular, or still a critical candidate under the full root-affine uncertainty.
- Why it was reasonable: the local-collar ladder had failed in fixed axes, gradient-normal axes, third-order collars, and branch-centered moving frames. The next proof object needed to identify whether the residual domain actually contains regular subregions before asking a collar theorem to carry every piece.
- Observed result: processed 64 residual regions. The diagnostic found 8 regular regions, all x-dominant, 0 excluded regions, and 56 critical candidates. The minimum regular gradient lower-bound candidate was 33.07361176255703. Length remained `20.316451752723314` under cap `20.672796062619668` with margin `0.3563443098963539`; no new branch length was promoted.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating all 64 L18 residuals as equally critical under the current interval model; eight have a robust x-direction derivative certificate target.
- What it does not rule out: it does not close the hard-cell certificate, does not resolve the 56 critical candidates, and does not prove global n=14 length or Erdos #114.
- Next dependency: prove the eight regular regions as a separate regular-slice branch atlas, then attack the 56 critical candidates with a distinct critical-candidate theorem.
- Claim ceiling: local hard-cell residual classification diagnostic only; not a claim upgrade.
- Orthodox reader note: this is the first post-collar narrowing that separates a certifiable regular slice from the genuinely hard critical-candidate slice. A skeptical reader can now audit two smaller obligations rather than one opaque failure bucket.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_global_critical_point_exclusion_target`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-GLOBAL-CRITICAL-POINT-EXCLUSION-TARGET-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L22 - Regular Residual Decomposition Diagnostic

- Date: 2026-05-06.
- Method or route: regular residual decomposition diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_REGULAR_RESIDUAL_DECOMPOSITION_TARGET_2026-05-06.md`
- Status: `REGULAR_RESIDUAL_FAIL_WALL_SEPARATION`.
- What we tried: built and ran a Rust diagnostic over exactly the eight L21 x-dominant regular regions, checking `Fx` sign stability, x-wall sign separation, ownership uniqueness, and graph-length budget while leaving the 56 critical candidates untouched.
- Why it was reasonable: L21 proved these eight residuals have robust x-direction derivative lower bounds. The orthodox next move was to test whether they already form owned `x=f(y)` collars before attempting the harder critical-candidate theorem.
- Observed result: processed 8 regular regions. `Fx` was sign-stable in all 8, derivative sign failures were 0, ownership duplicates were 0, but wall separation failures were 8 and certified regular regions were 0. Length stayed `20.316451752723314` under cap `20.672796062619668`; no new branch length was promoted.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out raw inherited x-wall interval evaluation as a sufficient certificate for the eight regular regions.
- What it does not rule out: it does not rule out regularity of the eight regions, a dependency-aware Taylor or affine wall certificate, the separate 56-region critical-candidate theorem, or the favorable local length budget.
- Next dependency: replace raw wall intervals with a dependency-aware Taylor/affine wall separation certificate for the eight regular regions.
- Claim ceiling: local hard-cell regular-slice diagnostic only; not a claim upgrade.
- Orthodox reader note: this is a useful failure: the derivative theorem is alive, but the wall theorem is too crude. A skeptical reader should see this as separation of proof obligations rather than retreat.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_regular_residual_decomposition`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-REGULAR-RESIDUAL-DECOMPOSITION-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L23 - Regular-Slice Taylor-y Wall Separation Diagnostic

- Date: 2026-05-06.
- Method or route: regular-slice Taylor-y wall separation diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_REGULAR_SLICE_WALL_TAYLOR_TARGET_2026-05-06.md`
- Status: `REGULAR_RESIDUAL_FAIL_WALL_SEPARATION`.
- What we tried: ran dependency-aware Taylor-y wall evaluation on the eight L22 regular regions. The diagnostic evaluated each x-wall by centering in y and using `Fy` over the wall strip, instead of direct raw interval evaluation over the full wall.
- Why it was reasonable: L22 showed derivative regularity survived but raw wall intervals were too wide. Taylor-y was the cheapest rigorous dependency repair before moving to full affine arithmetic.
- Observed result: processed 8 regular regions. Taylor-y reduced aggregate wall width by `14.679106803041705`, but certified 0 wall-separated branch collars and left 8 wall failures under the original branch-collar pass rule.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out simple Taylor-y wall separation as a branch-collar certificate for the eight regular boxes.
- What it does not rule out: it does not rule out monotone exclusion, full affine wall arithmetic, critical-candidate exclusion, or the favorable length budget.
- Next dependency: use the Taylor-y sign information with the stable `Fx` monotonicity rule; if both walls are same-signed, classify the region as excluded rather than failed branch certification.
- Claim ceiling: local hard-cell Taylor-y wall diagnostic only; not a claim upgrade.
- Orthodox reader note: this row is a methodological pivot point: Taylor-y did not make collars, but it exposed same-sign walls that become useful after monotonicity is accounted for.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_regular_residual_decomposition`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-REGULAR-SLICE-WALL-TAYLOR-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L24 - Regular-Slice Monotone Taylor Exclusion Diagnostic

- Date: 2026-05-06.
- Method or route: regular-slice monotone Taylor exclusion diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01_REPORT.md`
- Status: `REGULAR_RESIDUAL_PASS_NOT_GLOBAL_PROOF`.
- What we tried: re-ran the eight regular regions with Taylor-y wall bounds plus the stable `Fx` monotonicity rule. A region was closed if the chosen wall signs and monotone `Fx` excluded level-set crossing.
- Why it was reasonable: the L23 Taylor-y wall intervals were strictly positive on both proof walls while `Fx` was strictly negative. Monotonicity then implies `F` stays positive across each regular box, so these boxes are exclusions rather than branch collars.
- Observed result: processed 8 regular regions. Excluded 8 by monotone x-chart Taylor walls, with 0 ownership duplicates and 0 remaining Taylor wall failures. No branch length was promoted; length stayed `20.316451752723314` under cap `20.672796062619668`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it closes the eight L21 regular regions as non-crossing residual buckets under Taylor-y monotone exclusion.
- What it does not rule out: it does not address the 56 L21 critical candidates and does not prove global n=14 or Erdos #114.
- Next dependency: move to the separate 56-region critical-candidate theorem target; do not revisit the eight regular regions unless auditing constants.
- Claim ceiling: local hard-cell regular-slice exclusion diagnostic only; not a claim upgrade.
- Orthodox reader note: this is the first clean closure of a residual sub-obligation: the regular slice is not a missing length contribution, it is safely excluded by monotonicity and dependency-aware wall signs.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_regular_residual_decomposition`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-REGULAR-SLICE-MONOTONE-TAYLOR-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L25 - Critical-Candidate Affine Gradient Diagnostic

- Date: 2026-05-06.
- Method or route: critical-candidate affine gradient diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_CRITICAL_CANDIDATE_THEOREM_TARGET_2026-05-06.md`
- Status: `CRITICAL_CANDIDATE_AFFINE_GRADIENT_FAILS`.
- What we tried: ran a two-root-parameter affine-style gradient diagnostic over exactly the 56 L21 critical candidates, with parameter subdivision 4 and no length promotion.
- Why it was reasonable: the 56 candidates were critical only under full root-affine interval dependency. A first-order affine parameter model was the cheapest way to test whether the zero-gradient obstruction was mostly dependency swelling.
- Observed result: processed 56 candidates. It closed 0 by F exclusion, 0 by affine-gradient regularity, and left 56 still critical. Length stayed `20.316451752723314` under cap `20.672796062619668`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out first-order affine tracking of the full gradient product `p_prime * conjugate(p)` as sufficient for the critical candidates.
- What it does not rule out: it does not rule out direct p-prime exclusion, derivative-root location, stronger affine/Taylor models, or the favorable length budget.
- Next dependency: use the analytic invariant that on `|p|=1`, `grad F` can vanish only when `p_prime` vanishes; test `p_prime` directly rather than the wider product.
- Claim ceiling: local hard-cell affine-gradient diagnostic only; not a claim upgrade.
- Orthodox reader note: this is an honest negative control: preserving root-parameter dependency in the gradient product still leaves all candidates critical, which motivates the sharper p-prime invariant.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_critical_candidate_affine_gradient`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-AFFINE-GRADIENT-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L26 - Direct p-Prime Exclusion Diagnostic, Parameter Subdivision 4

- Date: 2026-05-06.
- Method or route: direct p-prime exclusion diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01_REPORT.md`
- Status: `CRITICAL_CANDIDATE_AFFINE_GRADIENT_PARTIAL` with `test_mode = pprime`.
- What we tried: reused the affine parameter model but tested `p_prime` directly on the 56 critical candidates. A tile is regular if the complex `p_prime` interval avoids zero.
- Why it was reasonable: on the lemniscate `|p|=1`, `p` is nonzero, so `grad(|p|^2-1)=0` implies `p_prime=0`. Direct p-prime exclusion avoids multiplying by `conjugate(p)`, which was widening the gradient test.
- Observed result: at parameter subdivision 4, direct p-prime exclusion closed 8 of 56 candidates and left 48 still critical. No branch length was promoted.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the idea that all 56 candidates are equally hard under direct p-prime analysis; two ownership families close cleanly.
- What it does not rule out: it does not close the remaining 48 candidates and does not prove the hard-cell certificate.
- Next dependency: repeat with parameter subdivision 8 to see whether additional p-prime closures appear before changing theorem class.
- Claim ceiling: local hard-cell p-prime diagnostic only; not a claim upgrade.
- Orthodox reader note: this is the first positive critical-candidate signal: the right invariant is p-prime location, not the full gradient product.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_critical_candidate_affine_gradient`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L27 - Direct p-Prime Exclusion Diagnostic, Parameter Subdivision 8

- Date: 2026-05-06.
- Method or route: direct p-prime exclusion diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01_REPORT.md`
- Status: `CRITICAL_CANDIDATE_PPRIME_PARTIAL`.
- What we tried: ran the direct p-prime exclusion diagnostic again with parameter subdivision 8 under a new immutable artifact.
- Why it was reasonable: the subdivision-4 p-prime run had real signal. A single stronger parameter split was justified to distinguish parameter dependency from a deeper root-location obstruction.
- Observed result: subdivision 8 closed the same 8 of 56 candidates and left 48 still critical. The closed families are `2484:1/root/y*` and `3707:1/root/y*`. The first remaining obstruction is `2587:3/root/y0/y0`. No length was promoted.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out simple parameter tiling beyond p=4 as the next route; p=8 does not increase closure count.
- What it does not rule out: it does not rule out a derivative-root location theorem, Rouche-style p-prime exclusion, or spatial/root-location decomposition for the remaining 48.
- Next dependency: move from p-prime interval evaluation to p-prime root-location: locate or exclude derivative roots relative to the remaining 48 boxes.
- Claim ceiling: local hard-cell p-prime diagnostic only; not a claim upgrade.
- Orthodox reader note: the useful result is not closure count alone; it identifies the invariant and shows further tiling is not the proof idea.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_critical_candidate_affine_gradient`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CRITICAL-CANDIDATE-PPRIME-EXCLUSION-P8-HARD-CELL-20260506-01_RESULTS.sha256`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L28 - p-Prime Root-Location Diagnostic

- Date: 2026-05-06.
- Method or route: p-prime root-location diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_PPRIME_ROOT_LOCATION_TARGET_2026-05-06.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_pprime_root_location.rs`
- Status: `PPRIME_ROOT_LOCATION_PARTIAL`.
- What we tried: built and ran the p-prime root-location diagnostic over exactly the 48 L27 remaining critical candidates, using a root-free disk test based on lower_abs(p_prime at the box center) versus spatial_radius times upper_abs(p_second over the box).
- Why it was reasonable: on the lemniscate `|p| = 1`, critical points require `p_prime(z) = 0`. This compresses the blocker from branch-collar bookkeeping into derivative-root exclusion.
- Observed result: the diagnostic processed 48 candidates, excluded p-prime roots on 31 candidates, and left 17 still near possible p-prime roots. The accepted length upper stayed `20.316451752723314` below cap `20.672796062619668`, with margin `0.3563443098963539`. The first failed condition is `3102:1/root/y0/y0`. The worst root-box distance lower is `-0.02296883128931693` at `4488:1/root/y1/y1`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out 31 of the remaining 48 L27 candidates as possible p-prime root locations under the root-free disk inequality.
- What it does not rule out: it does not close the hard-cell certificate, does not resolve the 17 still-near candidates, and does not rule out sharper analytic p-prime root-location or Taylor/Rouche-style arguments.
- Next dependency: build a sharper p-prime root-location theorem for the 17 still-near candidates, starting from the failed margins recorded in the L28 artifact.
- Claim ceiling: local n=14 hard-cell p-prime root-location diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
- Orthodox reader note: this is the first large reduction of the critical-candidate pile by the invariant a conventional analyst would recognize: critical points on the unit lemniscate must be zeros of `p_prime`.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_pprime_root_location`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.sha256`
  - `jq '{status, processed_remaining_candidate_count, pprime_root_excluded_count, still_unresolved_count, margin_to_cap, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L28 artifact directory and `ehp114_n14_pprime_root_location.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L29 - Sharp p-Prime Root-Location Diagnostic

- Date: 2026-05-06.
- Method or route: sharp p-prime root-location diagnostic for the 17 L28 still-near candidates.
- Artifact ID: `EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_sharp_pprime_root_location.rs`
- Status: `SHARP_PPRIME_ROOT_LOCATION_PASS_NOT_GLOBAL_PROOF`.
- What we tried: built and ran a sharper p-prime root-location diagnostic over exactly the 17 L28 still-near candidates, preserving the L28 disk test and adding a Taylor/Rouche variation bound with local splitting to max depth 4.
- Why it was reasonable: L28 showed the remaining obstruction was p-second variation overwhelming the p-prime center lower bound. A Taylor/Rouche test with p-third remainder and local radius reduction attacks that exact failure mode.
- Observed result: processed 17 candidates and 2432 local leaves. All 17 candidates were separated from p-prime roots under the sharp test; unresolved leaf count was 0. The accepted length upper remained `20.316451752723314` below cap `20.672796062619668`, with margin `0.3563443098963539`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the 17 L28 still-near regions as p-prime root locations under the recorded sharp Taylor/Rouche local split test.
- What it does not rule out: it does not by itself prove Erdos #114, does not produce a global n=14 certificate, and does not replace the need for an integration check across the residual chain.
- Next dependency: run an integration certificate verifying that L24 regular exclusions, L27/L28 p-prime exclusions, and L29 sharp p-prime exclusions cover the intended residual chain without ownership gaps.
- Claim ceiling: local n=14 hard-cell sharp p-prime root-location diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
- Orthodox reader note: this is the first full closure of the derivative-root-near residual pile, but the orthodox next step is coverage integration rather than public claim promotion.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_sharp_pprime_root_location`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-SHARP-PPRIME-ROOT-LOCATION-HARD-CELL-20260506-01_RESULTS.sha256`
  - `jq '{status, processed_remaining_candidate_count, sharp_pprime_root_excluded_count, still_unresolved_count, processed_leaf_count, unresolved_leaf_count, max_depth, margin_to_cap, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L29 artifact directory and `ehp114_n14_sharp_pprime_root_location.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L30 - RH/MDL Beurling-Nyman Quantization-Stability Finite Audit

- Date: 2026-05-06.
- Method or route: RH/MDL Beurling-Nyman quantization-stability finite audit.
- Artifact ID: `EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/RH_BEURLING_NYMAN_MDL_QUANTIZATION_STABILITY_AUDIT_2026-05-06.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/rh_beurling_nyman_mdl_quantization_stability_check.py`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/EHP114_RH_MDL_NEXT_PHASE_SYNTHESIS_2026-05-06.md`
- Status: `BN_MDL_QUANTIZATION_STABILITY_CHECK_PASS_NOT_RH_EVIDENCE`.
- What we tried: created a finite checker for the frozen Beurling-Nyman MDL probe and sensitivity artifacts, requiring quantization step, quantization penalty, certified upper residual, condition number, and claim ceiling fields.
- Why it was reasonable: the Baez-facing critique identified invariance and bit-cost stability as the load-bearing gap. A finite quantization audit is the modest theorem-facing move before any asymptotic or RH-shaped language.
- Observed result: the checker verified source SHA status PASS, audited 98 rows and 312 quantized entries, and found 0 field failures. The 3.2-bit bend remains named only as the first observed finite-N MDL conditioning crossover.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the immediate objection that the current finite MDL rows hide missing quantization or condition-number fields.
- What it does not rule out: it does not establish dictionary invariance, asymptotic residual decay, theorem progress toward RH, or RH evidence.
- Next dependency: if this lane continues, formalize the finite quantization-stability inequality and stress-test quadrature/dictionary invariance under new immutable artifacts.
- Claim ceiling: finite Beurling-Nyman MDL quantization-stability audit only; not RH evidence, not theorem progress, and not an asymptotic statement.
- Orthodox reader note: this is the Hilbert-space MDL lane, not the #114 proof lane; it makes the finite object auditable without claiming RH significance.
- Verification commands:
  - `python3 rh_beurling_nyman_mdl_quantization_stability_check.py`
  - `shasum -a 256 -c EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_RESULTS.sha256`
  - `jq '{status, total_row_count, total_quantized_entry_count, total_failure_count, source_sha_status, bend_name, claim_ceiling}' EXP-MATH-RH-BEURLING-NYMAN-MDL-QUANTIZATION-STABILITY-CHECK-20260506-01_RESULTS.json`
  - overclaim scan over RH/MDL checker artifacts and synthesis packet.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L31 - Residual-Chain Integration Certificate

- Date: 2026-05-06.
- Method or route: residual-chain integration certificate over L24/L27/L28/L29.
- Artifact ID: `EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_residual_chain_integration_cert.rs`
- Status: `RESIDUAL_CHAIN_INTEGRATION_PASS_NOT_GLOBAL_PROOF`.
- What we tried: built and ran a bookkeeping integration certificate over L24 regular monotone exclusions, L27 p-prime P8 partition closures, L28 p-prime root-location exclusions, and L29 sharp p-prime root-location exclusions.
- Why it was reasonable: after L29 closed the derivative-root-near pile, the proof-facing risk shifted from local inequalities to coverage integrity: source checksums, source filters, candidate counts, ownership uniqueness, and length-cap consistency.
- Observed result: source SHA verification passed for L24, L27, L28, and L29. The chain closed 8 regular regions, 8 L27 partition-closed critical regions, 31 L28 p-prime root-location regions, and 17 L29 sharp root-location regions. Final critical closed count was 56/56, residual chain total was 64/64, with 0 ownership duplicates, 0 source-filter mismatches, 0 count drift, and margin to cap `0.3563443098963539`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the immediate bookkeeping objections that the L24/L27/L28/L29 local closures fail to compose because of source-filter drift, candidate-count drift, duplicate ownership, checksum failure, or length-budget inconsistency.
- What it does not rule out: it does not prove Erdos #114 globally, does not certify all n=14 subcells, and does not replace the need for an orthodox local hard-cell theorem packet and independent audit before public claim promotion.
- Next dependency: assemble a local hard-cell certificate packet that states the theorem shape, cites the accepted length source and L31 integration certificate, and then run independent review before scaling beyond subcell `(6,4)`.
- Claim ceiling: local n=14 hard-cell residual-chain integration certificate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate. shadow signature, not universal law.
- Orthodox reader note: this is the expected audit step after local exclusions: it proves the pieces compose as a residual-chain certificate, while keeping the claim local and below global n=14 or Erdos #114 language.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_residual_chain_integration_cert`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-HARD-CELL-20260506-01_RESULTS.sha256`
  - `jq '{status, regular_regions_closed_by_l24, critical_regions_closed_by_l27_partition, critical_regions_closed_by_l28_root_location, critical_regions_closed_by_l29_sharp_root_location, final_critical_closed_count, residual_chain_total_closed_count, ownership_duplicate_count, source_filter_mismatch_count, candidate_count_drift_count, margin_to_cap, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over sidecar notes, the L31 artifact directory, and `ehp114_n14_residual_chain_integration_cert.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L32 - Local Hard-Cell Certificate Packet

- Date: 2026-05-06.
- Method or route: local hard-cell certificate packet over L31 residual-chain integration.
- Artifact ID: `EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01/EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_local_hard_cell_certificate_packet.rs`
- Status: `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF`.
- What we tried: built and ran a theorem-facing local hard-cell packet that verifies the L31 residual-chain source, recomputes source checksums for L31 and its child L24/L27/L28/L29 artifacts, and binds the accepted length upper to the local cap.
- Why it was reasonable: after L31 passed, the next skeptical-reader objection was not another local inequality but whether the local hard-cell result had a single theorem-shaped packet with hashes, coverage conditions, length budget, and claim ceiling in one place.
- Observed result: the packet passed with source SHA fail count `0`, L31 failure count `0`, residual-chain total closed count `64`, total validated length upper `20.316451752723314`, exact length cap `20.672796062619668`, margin to cap `0.3563443098963539`, and first failed condition `none`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the immediate local-hard-cell packaging objection: the L31 closure and child artifact hashes are not merely scattered diagnostics but are bound into a theorem-facing local certificate packet.
- What it does not rule out: it does not certify all n=14 root-affine cells, does not make Tao high-degree bounds operational for finite completion, and does not prove the full EHP114 statement.
- Next dependency: run L33 global n=14 atlas coverage over all root-affine cells, using L32 as the hard-cell theorem packet and failing on any missing cell, duplicate ownership, source-filter mismatch, count drift, or budget violation.
- Claim ceiling: local n=14 hard-cell certificate packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the expected local theorem-wrapper step: the local hard cell is now packetized and auditable, while the global theorem remains blocked on all-cell n=14 coverage and a high-degree bridge.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_local_hard_cell_certificate_packet`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-LOCAL-HARD-CELL-CERTIFICATE-PACKET-20260506-01_RESULTS.sha256`
  - `jq '{status, local_hard_cell_certificate_pass, source_sha_fail_count, l31_failure_count, residual_chain_total_closed_count, total_validated_length_upper, exact_length_cap, margin_to_cap, first_failed_condition, next_dependency}' *_RESULTS.json`
  - overclaim scan over the L32 artifact directory and `ehp114_n14_local_hard_cell_certificate_packet.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33 - Global n=14 Atlas Certificate Gate

- Date: 2026-05-06.
- Method or route: global n=14 atlas certificate gate over all root-affine cells.
- Artifact ID: `EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01/EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_global_atlas_certificate.rs`
- Status: `GLOBAL_ATLAS_FAIL_MISSING_CELL_CERTIFICATES`.
- What we tried: built and ran an atlas gate that verifies the 64-cell root-affine coverage skeleton, scans for theorem-grade local cell packets, checks their source hashes and local acceptance conditions, and fails on any missing, duplicate, rejected, or over-budget cell.
- Why it was reasonable: after L32 packetized the hard cell, the next orthodox objection was local-to-global coverage: a single hard-cell theorem packet cannot imply a global n=14 atlas certificate without one accepted local packet for every root-affine cell.
- Observed result: the root-affine skeleton enumerated `64/64` cells with source SHA failures `0`. The gate accepted `1` theorem-grade cell packet, subcell `(6,4)`, and found `63` missing cell certificates. Coverage failures, duplicate cell certificates, rejected discovered packets, and budget failures were all `0`. First failed condition: missing theorem-grade local cell certificate for subcell `(0,0)`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating the L32 hard-cell packet as a global n=14 certificate. It also checks that the 64-cell skeleton itself is present and not the immediate source of coverage drift.
- What it does not rule out: it does not rule out certifying the missing 63 cells by generating theorem-grade residual-chain packets for them; this is a missing-certificate integration failure, not a geometric counterexample.
- Next dependency: generate or import theorem-shaped local residual-chain packets for the missing 63 n=14 root-affine cells, then rerun the L33 atlas gate.
- Claim ceiling: L33 global atlas integration gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the correct failure mode: the global proof ladder now has a machine-readable list of missing local theorem packets rather than a vague scaling TODO.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_global_atlas_certificate`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-GLOBAL-ATLAS-CERTIFICATE-20260506-01_RESULTS.sha256`
  - `jq '{status, global_atlas_certificate_pass, root_affine_skeleton_cell_count, certified_cell_count, missing_cell_certificate_count, source_sha_fail_count, coverage_failure_count, duplicate_cell_certificate_count, rejected_cell_certificate_count, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33 artifact directory and `ehp114_n14_global_atlas_certificate.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33A - Missing-Cell Local Packet Generation Gate

- Date: 2026-05-06.
- Method or route: missing-cell local packet generation gate and pipeline parameterization audit.
- Artifact ID: `EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02/EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_missing_cell_local_packets.rs`
- Status: `MISSING_CELL_LOCAL_PACKETS_BLOCKED_PIPELINE_PARAMETERIZATION`.
- What we tried: built and ran a missing-cell packet gate that consumes L33, verifies its checksum, emits the 63-cell worklist, and audits whether the local residual-chain binaries can honestly run on cells other than `(6,4)`. The initial `-01` audit used an overly loose CLI detector and was versioned forward; `-02` is the canonical corrected artifact.
- Why it was reasonable: after L33 found 63 missing local packets, the immediate risk was accidentally renaming the hard-cell pipeline output for other cells. The correct next gate is to prove the pipeline is parameterized before launching a batch run.
- Observed result: L33 checksum passed and the missing-cell worklist contains `63` cells. The audit found `7/7` required binaries blocked by hard-coded `(6,4)` constants: `ehp114_n14_global_critical_point_exclusion_target`, `ehp114_n14_regular_residual_decomposition`, `ehp114_n14_critical_candidate_affine_gradient`, `ehp114_n14_pprime_root_location`, `ehp114_n14_sharp_pprime_root_location`, `ehp114_n14_residual_chain_integration_cert`, and `ehp114_n14_local_hard_cell_certificate_packet`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out starting L33 cell generation by simply replaying the current residual-chain binaries, because those binaries are not yet parameterized by subcell.
- What it does not rule out: it does not rule out certifying the 63 missing cells; it only shows that the software interface must be repaired before proof-facing batch generation can begin.
- Next dependency: parameterize the seven required L21/L24/L27/L28/L29/L31/L32 binaries with `--sub-i`, `--sub-j`, per-cell experiment IDs, per-cell output directories, and per-cell source-path wiring; then rerun L33A until it returns ready to run.
- Claim ceiling: missing-cell local-packet generation gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is a useful software-proof boundary: the audit prevents a coverage claim from being manufactured by path changes while constants still point at the hard cell.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_missing_cell_local_packets`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-MISSING-CELL-LOCAL-PACKETS-20260506-02_RESULTS.sha256`
  - `jq '{status, missing_cell_count, source_sha_fail_count, required_binary_count, blocking_hardcoded_binary_count, blockers:[.pipeline_audit[] | select(.blocking) | {binary_name, role, evidence}], first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33A artifact directory and `ehp114_n14_missing_cell_local_packets.rs`.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33B - Parameterized Local-Pipeline Smoke

- Date: 2026-05-06.
- Method or route: parameterized local residual-chain pipeline smoke gate.
- Artifact ID: `EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01/EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01/EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/ehp114_n14_cell.rs`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_parameterized_local_pipeline_smoke.rs`
- Status: `PARAMETERIZED_LOCAL_PIPELINE_BLOCKED_SOURCE_ARTIFACT_MISSING`.
- What we tried: added a shared n=14 `CellSpec` helper, parameterized the seven L21/L24/L27/L28/L29/L31/L32 local residual-chain binaries with `--sub-i`, `--sub-j`, and `--experiment-id`, added source/subcell contract checks, built the target binaries, and ran a `CELL-00-00` smoke gate.
- Why it was reasonable: L33A showed the immediate proof risk was software reuse of the hard-cell constants. A skeptic does not need another local inequality until the pipeline can prove which cell its inputs describe.
- Observed result: the smoke artifact reports `hardcoded_blocker_count = 0`: all seven required binaries expose subcell CLI, have no fixed `SUB_I` / `SUB_J` / `SUBCELL` constants, and use source-subcell checks. The deliberate `CELL-00-00` probe against the hard-cell source reports `SOURCE_SUBCELL_MISMATCH`, and the overall smoke status is blocked because the per-cell upstream source artifact for `CELL-00-00` does not yet exist.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the L33A software-proof blocker: the local residual-chain binaries can no longer silently ignore the requested subcell or reuse hard-cell data without a source/subcell mismatch.
- What it does not rule out: it does not certify any of the 63 missing cells, does not generate per-cell upstream source artifacts, and does not make L33 global n=14 coverage pass.
- Next dependency: generate or parameterize the upstream per-cell source artifacts, starting with the L18/L21 source chain for `CELL-00-00`, then run the parameterized local pipeline on that cell.
- Claim ceiling: parameterized local-pipeline smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the kind of software-proof guard a skeptical reader expects: the code now fails on source mismatch instead of letting a path rename masquerade as new cell coverage.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_parameterized_local_pipeline_smoke --bin ehp114_n14_global_critical_point_exclusion_target --bin ehp114_n14_regular_residual_decomposition --bin ehp114_n14_critical_candidate_affine_gradient --bin ehp114_n14_pprime_root_location --bin ehp114_n14_sharp_pprime_root_location --bin ehp114_n14_residual_chain_integration_cert --bin ehp114_n14_local_hard_cell_certificate_packet`
  - `./target/release/ehp114_n14_global_critical_point_exclusion_target --sub-i 8 --outdir /tmp/ehp114-invalid-subcell-smoke`
  - `./target/release/ehp114_n14_global_critical_point_exclusion_target --sub-i 0 --sub-j 0 --experiment-id EXP-MATH-EHP114-N14-SOURCE-CONTRACT-SMOKE-CELL-00-00-20260506-01 --outdir /tmp/ehp114-source-contract-smoke`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01_RESULTS.sha256`
  - `jq '{status, smoke_subcell, hardcoded_blocker_count, source_contract_status, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33B artifact directory and patched L21/L24/L27/L28/L29/L31/L32 Rust binaries.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33C - Per-Cell Source Generation Smoke

- Date: 2026-05-06.
- Method or route: per-cell source-generation smoke gate for `CELL-00-00`.
- Artifact ID: `EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-02`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-01/EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-02/EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-02_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_per_cell_source_generation_smoke.rs`
- Status: `PER_CELL_SOURCE_CHAIN_READY`.
- What we tried: parameterized the direct upstream source-chain binaries feeding L21, added the L33C smoke gate, and ran `CELL-00-00` as the first missing-cell source test without synthesizing substitute data.
- Why it is reasonable: the parameterized L21/L24/L27/L28/L29/L31/L32 pipeline still needs honest source artifacts for each requested cell. Starting with `CELL-00-00` tests the first missing cell without launching a 63-cell batch.
- Observed result: the initial `-01` smoke correctly stopped on the missing `CELL-00-00` slab source. After L33D emitted that source, the `-02` smoke reports `hardcoded_blocker_count = 0`, `source_contract_status = SOURCE_SUBCELL_MATCH`, `first_missing_source_artifact = null`, and `first_failed_condition = none`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out both prior software blockers for the first missing cell: the direct source-chain binaries are parameterized and the first slab source artifact now exists with matching `CELL-00-00` metadata.
- What it does not rule out: it does not certify `CELL-00-00`; it only says the source chain can now start from a real cell-local source object.
- Next dependency: evaluate whether the `CELL-00-00` slab source is proof-usable. If it is over budget or unresolved, sharpen the source before downstream L21-L32 work.
- Claim ceiling: per-cell source-generation smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the right kind of failure: the proof pipeline refuses to move from hard-cell data to another cell until the source artifact itself carries that cell identity.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_per_cell_source_generation_smoke --bin ehp114_n14_branch_isolation_collar_atlas --bin ehp114_n14_normal_collar_critical_exclusion_pilot --bin ehp114_n14_third_order_collar_remainder_pilot --bin ehp114_n14_global_critical_point_exclusion_target`
  - `./target/release/ehp114_n14_branch_isolation_collar_atlas --sub-i 8 --outdir /tmp/ehp114-l33c-invalid-subcell`
  - `./target/release/ehp114_n14_branch_isolation_collar_atlas --sub-i 0 --sub-j 0 --experiment-id EXP-MATH-EHP114-N14-L33C-SOURCE-CONTRACT-SMOKE-CELL-00-00-20260506-01 --outdir /tmp/ehp114-l33c-source-contract-smoke`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-PER-CELL-SOURCE-GENERATION-SMOKE-20260506-02_RESULTS.sha256`
  - `jq '{status, smoke_subcell, hardcoded_blocker_count, source_contract_status, first_missing_source_artifact, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33C artifact directory and touched upstream Rust binaries.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33D - CELL-00-00 Slab Source Generation

- Date: 2026-05-06.
- Method or route: `CELL-00-00` slab validated-length source generation.
- Artifact ID: `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_validated_length_slab.rs`
- Status: `SLAB_FAIL_ROOT_ISOLATION`.
- What we tried: parameterized the slab validated-length generator with explicit subcell and experiment-id arguments, then ran z32/local-bound-6 for `CELL-00-00` to create the first real non-hard-cell source artifact.
- Why it is reasonable: the local residual-chain pipeline now has cell-aware contracts, but the first source consumed by branch isolation is still the slab validated-length artifact. A real `CELL-00-00` slab source is the shortest invariant-preserving next step.
- Observed result: the artifact has matching subcell metadata `CELL-00-00` and SHA verification passes. It reports `slab_branch_count = 2427`, `unresolved_branch_count = 4412`, `ownership_duplicate_count = 0`, `total_validated_length_upper = 30.658351028561448`, exact length cap `20.672796062619668`, and `margin_to_cap = -9.98555496594178`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out source absence as the blocker for `CELL-00-00` and shows that the current z32 slab source is not proof-usable for this cell because the certified branch-length upper already exceeds the cap and root isolation remains unresolved.
- What it does not rule out: it does not rule out the true `CELL-00-00` length being below cap; it only retires this coarse z32 slab source as a sufficient proof anchor for that cell.
- Next dependency: sharpen the `CELL-00-00` source bound before downstream L21-L32 work: either tighter slab/vertical isolation, rotated/normal chart source generation, or a source-quality diagnostic that explains the `30.658` overage.
- Claim ceiling: `CELL-00-00` slab source diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is now a mathematical/computational sharpness blocker, not a metadata blocker: the cell source exists, but its first length bound is too coarse to feed a theorem packet.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_validated_length_slab`
  - `./target/release/ehp114_n14_validated_length_slab --sub-i 0 --sub-j 0 --experiment-id EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01 --outdir <artifact-dir> --z-subdivision 32 --local-bound-subdivision 6`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-00-20260506-01_RESULTS.sha256`
  - `jq '{experiment_id,status,subcell,cell_tag,total_validated_length_upper,exact_length_cap:.length_budget.exact_length_cap,margin_to_cap,slab_branch_count,unresolved_branch_count,ownership_duplicate_count}' *_RESULTS.json`
  - overclaim scan over the L33D artifact directory and slab generator binary.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33E - CELL-00-00 Source Sharpening Diagnostic

- Date: 2026-05-06.
- Method or route: `CELL-00-00` slab source sharpening diagnostic.
- Artifact ID: `EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01_REPORT.md`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/scripts/erdos-114/src/bin/ehp114_n14_cell00_slab_source_sharpening.rs`
- Status: `CELL00_SLAB_SOURCE_SHARPENING_DIAGNOSTIC_FAIL_OVER_BUDGET`.
- What we tried: built and ran a source-quality diagnostic comparing the `CELL-00-00` slab source against the hard-cell slab baseline. The diagnostic computes length, slope, denominator, y-run, top-branch, top-slab, and unresolved-branch concentration statistics.
- Why it is reasonable: the source-chain metadata problem is fixed, but the source length upper is already above cap. Downstream residual closure only adds length in the current architecture, so the proof-facing bottleneck must be attacked at the source-bound level first.
- Observed result: the source remains over budget: `total_validated_length_upper = 30.658351028561448` against cap `20.672796062619668`. The overage is concentrated: `top_10_branch_sum = 17.77377731502446`, which exceeds the `9.98555496594178` cap excess. The top three branches have slopes `8467.56`, `5667.88`, and `3759.27` with denominator lower bounds `0.007006`, `0.001896`, and `0.006212`. Recommended next route is `targeted_high_slope_slab_repair_before_downstream_pipeline`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out treating the `CELL-00-00` overage as a diffuse global failure at the current diagnostic level. The first repair target is a small set of high-slope, low-denominator slabs rather than a broad 63-cell batch.
- What it does not rule out: it does not certify `CELL-00-00`, does not prove that targeted repair will close the source, and does not rule out a deeper geometric theorem being needed if the high-slope slabs cannot be repaired.
- Next dependency: run targeted high-slope slab repair on the top slabs, starting with `ix = 3713`, `4063`, and `2584`, before downstream L21-L32 promotion.
- Claim ceiling: `CELL-00-00` slab source-sharpening diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the correct next narrowing: the source bound is over cap, but the excess is concentrated enough that a targeted local repair has a falsifiable path.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_cell00_slab_source_sharpening`
  - `./target/release/ehp114_n14_cell00_slab_source_sharpening`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CELL-00-00-SLAB-SOURCE-SHARPENING-20260506-01_RESULTS.sha256`
  - `jq '{status, source_cell:.source_stats.cell_tag, source_total:.source_stats.total_validated_length_upper, cap:.source_stats.exact_length_cap, margin:.source_stats.margin_to_cap, top_10_branch_sum:.overage_analysis.top_10_branch_sum, recommended_next_route, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33E artifact directory and source-sharpening binary.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33F - CELL-00-00 High-Slope Slab Repair

- Date: 2026-05-06.
- Method or route: targeted high-slope slab repair for `CELL-00-00`.
- Artifact ID: `EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01_REPORT.md`
- Status: `HIGH_SLOPE_SLAB_REPAIR_ADJUSTED_UNDER_CAP_DIAGNOSTIC`.
- What we tried: built and ran a branch-local repair over the three largest high-slope source branches: `2584:3`, `3713:1`, and `4063:0`. Each target x-slab was split into 32 subslabs; each original y-tube was rescanned into 256 bins; candidate bins were regrouped into mini-collars and recertified by endpoint sign separation plus a stable derivative chart.
- Why it is reasonable: L33E shows the overage is concentrated in a small set of high-slope/low-denominator slab branches. Repairing those slabs is the shortest invariant-preserving route before downstream local-pipeline promotion.
- Observed result: all three target branches certify with zero unresolved repair groups. Their original length sum was `15.251173022324249`; repaired target length upper is `0.003308273015128627`; replacement accounting gives an adjusted diagnostic total `15.410486279252329` against cap `20.672796062619668`, with margin `5.262309783367339`. The source still has `4412` unresolved branches, so this is not a cell certificate.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out the largest `CELL-00-00` overage terms as intrinsic branch length. They were artifacts of coarse high-slope slab tubes with tiny derivative denominators.
- What it does not rule out: it does not certify `CELL-00-00`, close the remaining unresolved branches, produce a theorem-grade source overlay, cover the other 63 cells, or complete global n=14 coverage.
- Next dependency: materialize a repaired source-overlay artifact that records branch replacement accounting and can be consumed by downstream L21-L32 without silently reusing the stale over-budget source.
- Claim ceiling: `CELL-00-00` high-slope slab repair diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is a useful correction, not a proof step by itself. It shows the overage was mostly a coordinate/tube artifact, but the repaired accounting must be bound into a source contract before downstream certificates may use it.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_cell00_high_slope_slab_repair`
  - `./target/release/ehp114_n14_cell00_high_slope_slab_repair`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CELL-00-00-HIGH-SLOPE-SLAB-REPAIR-20260506-01_RESULTS.sha256`
  - `jq '{status, target_count, target_original_length_sum, repaired_target_length_upper, target_length_delta, adjusted_total_if_replaced, exact_length_cap, adjusted_margin_to_cap, target_unresolved_group_count, remaining_source_unresolved_branch_count, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33F artifact directory and high-slope slab repair binary.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33G - CELL-00-00 Repaired Source Overlay

- Date: 2026-05-06.
- Method or route: repaired source-overlay materialization for `CELL-00-00`.
- Artifact ID: `EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01/EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01_REPORT.md`
- Status: `REPAIRED_SLAB_SOURCE_OVERLAY_READY_FOR_DOWNSTREAM_DIAGNOSTIC`.
- What we tried: materialized L33F replacement accounting into a source-like overlay, replacing the three repaired branch lengths while preserving source subcell metadata, unresolved branch rows, ownership metadata, and SHA audit.
- Why it is reasonable: L33F gives verified branch-local replacement accounting, but downstream proof-facing binaries must not infer replacement from prose. They need a source artifact or overlay with explicit repaired branch keys, replacement length, remaining unresolved list, ownership policy, and SHA.
- Observed result: overlay status is ready for downstream diagnostics. It records `replacement_count = 3`, source total `30.658351028561448`, repaired total `15.410486279252329`, cap `20.672796062619668`, margin `5.262309783367339`, and `unresolved_branch_count = 4412`. Rerun source-chain smoke `20260506-03` passes with `SOURCE_SUBCELL_MATCH` and zero hardcoded blockers.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out narrative replacement accounting as the blocker. The repaired source is now a machine-readable input contract.
- What it does not rule out: it does not close `CELL-00-00`; it only gives downstream diagnostics a usable source object.
- Next dependency: run the downstream residual-chain pipeline against the repaired overlay.
- Claim ceiling: `CELL-00-00` repaired source-overlay route only; no claim upgrade.
- Orthodox reader note: replacement accounting must be explicit and machine-readable; otherwise this would be a narrative patch rather than a certificate input.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_cell00_repaired_slab_source_overlay`
  - `./target/release/ehp114_n14_cell00_repaired_slab_source_overlay`
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-CELL-00-00-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01_RESULTS.sha256`
  - `./target/release/ehp114_n14_per_cell_source_generation_smoke --source <overlay-results> --sub-i 0 --sub-j 0`
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33H - CELL-00-00 Local Certificate Packet

- Date: 2026-05-06.
- Method or route: first non-hard-cell residual-chain and local certificate packet for `CELL-00-00`.
- Artifact ID: `EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-00-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-00-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-00-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-00-20260506-01_REPORT.md`
- Status: `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF`.
- What we tried: ran the downstream local chain from the repaired overlay: branch atlas, normal collar, third-order collar, global critical-point target, regular monotone Taylor closure, direct p-prime closure, p-prime root-location, sharp p-prime root-location, residual-chain integration, and local packet assembly.
- Why it is reasonable: after L33G made the repaired source explicit and under cap, the honest test was whether the already parameterized residual-chain architecture could close the first non-hard cell without source/count/ownership drift.
- Observed result: local packet passes. Residual-chain integration closes `64/64`, with `8` regular regions closed by L24 and `56` critical regions closed across L27/L28/L29. Source SHA fail count is `0`, source-filter mismatch count is `0`, candidate-count drift is `0`, ownership duplicate count is `0`, total validated length upper is `15.410486279252329`, cap is `20.672796062619668`, and margin is `5.262309783367339`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out `CELL-00-00` as merely a software smoke success; this cell now has a local theorem-packet-shaped artifact.
- What it does not rule out: it does not certify any of the other 62 missing cells, does not prove global n=14 coverage, and does not complete EHP114.
- Next dependency: repeat the source repair plus residual-chain packet on the next missing cell under the same source/subcell contracts before broad batching.
- Claim ceiling: local `CELL-00-00` packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the first evidence that the hard-cell residual-chain architecture transfers to a non-hard cell, but one transferred cell is not a global atlas.
- Verification commands:
  - `shasum -a 256 -c *_RESULTS.sha256` for the L33G/L33H artifact family.
  - `jq '{status, local_hard_cell_certificate_pass, residual_chain_total_closed_count, source_sha_fail_count, total_validated_length_upper, exact_length_cap, margin_to_cap, first_failed_condition}' *_RESULTS.json`
  - overclaim scan over the L33G/L33H artifacts and patched binaries.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33I - CELL-00-01 Controlled Next-Cell Generation

- Date: 2026-05-06.
- Method or route: controlled next-cell generation for `CELL-00-01`.
- Artifact ID: `EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-01-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-01-20260506-01/EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-CELL-00-01-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-01-HIGH-SLOPE-SLAB-REPAIR-20260506-01/EXP-MATH-EHP114-N14-CELL-00-01-HIGH-SLOPE-SLAB-REPAIR-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-CELL-00-01-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01/EXP-MATH-EHP114-N14-CELL-00-01-REPAIRED-SLAB-SOURCE-OVERLAY-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-00-01-20260506-01/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-00-01-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-01-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-01-20260506-01_RESULTS.json`
- Status: `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF`.
- What we tried: generated the real z32 slab source for `CELL-00-01`, verified the source contract, diagnosed a small over-cap source bound, repaired the three dominant high-slope branch contributions, materialized a repaired source overlay, and ran the L21/L24/L27/L28/L29/L31/L32 local packet chain.
- Why it is reasonable: `CELL-00-00` showed that the pipeline can transfer after explicit source repair. A second controlled non-hard cell tests reproducibility without hiding failures inside a broad remaining-cell batch.
- Observed result: the raw slab source was over cap by `0.5066662394395678`, but the top branch lengths were concentrated. Repairing targets `2584:3`, `4063:0`, and `2412:1` reduced the adjusted total to `15.078165474813495`. The residual chain then closed `64/64` with no source SHA failures, no ownership duplicates, no source-filter mismatch, and no candidate-count drift.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out `CELL-00-01` as a one-off source-transfer blocker and shows that the high-slope repair pattern transfers to a second non-hard cell.
- What it does not rule out: it does not prove global n=14 coverage, does not certify the remaining missing cells, and does not replace the eventual global atlas certificate.
- Next dependency: move from single-cell control to a small representative-cell harness over 3-5 additional missing cells before attempting full remaining-cell generation.
- Claim ceiling: local `CELL-00-01` packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: the useful fact is reproducibility: two non-hard cells now pass through the same source-repair and residual-chain structure, but the burden is still coverage across the remaining cells.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_validated_length_slab --bin ehp114_n14_per_cell_source_generation_smoke --bin ehp114_n14_residual_chain_integration_cert --bin ehp114_n14_local_hard_cell_certificate_packet`
  - `find /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length -maxdepth 2 -path '*CELL-00-01*/*_RESULTS.sha256' -print0 | while IFS= read -r -d '' file; do dir=$(dirname "$file"); base=$(basename "$file"); (cd "$dir" && shasum -a 256 -c "$base"); done`
  - `jq '{status, subcell, residual_chain_total_closed_count, source_sha_fail_count, ownership_duplicate_count, source_filter_mismatch_count, candidate_count_drift_count, total_validated_length_upper, exact_length_cap, margin_to_cap, first_failed_condition}' *_RESULTS.json`
  - forbidden-claim scan over `CELL-00-01` artifacts and touched Rust binaries returned no matches.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33J - Representative-Cell Harness

- Date: 2026-05-06.
- Method or route: representative-cell harness over `CELL-00-02` and `CELL-03-03`.
- Artifact ID: `EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-20260506-01/EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-02-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-00-02-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-03-03-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-03-03-20260506-01_RESULTS.json`
- Status: `REPRESENTATIVE_CELL_HARNESS_BLOCKED_REGULAR_RESIDUAL_DRIFT`.
- What we tried: ran the local source and residual-chain pipeline on `CELL-00-02` and `CELL-03-03`, with planned `CELL-07-07` intentionally skipped after the first non-transferable blocker.
- Why it is reasonable: two non-hard cells already passed, but a small representative harness is the right gate before broad remaining-cell generation. `CELL-00-02` tests another edge-adjacent source; `CELL-03-03` tests an interior geometry.
- Observed result: `CELL-00-02` passed with raw source length `17.183471776637536` under cap and residual chain `64/64`. `CELL-03-03` had source length repaired from `22.71308313038541` to `18.01426357052187` under cap, and the derivative-root chain closed all critical candidates, but L24 saw `12` regular regions, excluded only `8`, and left `4` regular wall-separation failures. L31 correctly failed with count drift rather than promoting the packet.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out immediate broad remaining-cell generation under the current fixed eight-regular-region residual-chain assumption.
- What it does not rule out: it does not rule out `CELL-03-03` mathematically. The blocker is regular-region closure and integration generality, not length budget or p-prime critical-candidate closure.
- Next dependency: L33K should generalize L24/L31 to variable regular-region counts and add a proof-facing repair for the four `CELL-03-03` regular wall-separation failures before returning to representative-cell batching.
- Claim ceiling: representative-cell harness diagnostic only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the kind of failure a skeptical reader should want: the pipeline stopped when an interior cell changed the residual topology instead of silently forcing it into the old eight-region template.
- Verification commands:
  - `shasum -a 256 -c EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-20260506-01_RESULTS.sha256`
  - `find /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length -maxdepth 2 \( -path "*CELL-00-02*/*_RESULTS.sha256" -o -path "*CELL-03-03*/*_RESULTS.sha256" -o -path "*REPRESENTATIVE-CELL-HARNESS-20260506-01/*_RESULTS.sha256" \) -print0 | while IFS= read -r -d "" file; do dir=$(dirname "$file"); base=$(basename "$file"); (cd "$dir" && shasum -a 256 -c "$base"); done`
  - forbidden-claim scan over `CELL-00-02`, `CELL-03-03`, and L33J harness artifacts returned no matches.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33K - Variable Regular-Region Closure

- Date: 2026-05-06.
- Method or route: variable regular-region closure for `CELL-03-03`.
- Artifact ID: `EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-03-03-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-03-03-20260506-01/EXP-MATH-EHP114-N14-VARIABLE-REGULAR-RESIDUAL-CLOSURE-CELL-03-03-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-03-03-20260506-02/EXP-MATH-EHP114-N14-RESIDUAL-CHAIN-INTEGRATION-CERT-CELL-03-03-20260506-02_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-03-03-20260506-02/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-03-03-20260506-02_RESULTS.json`
- Status: `LOCAL_HARD_CELL_CERTIFICATE_PASS_NOT_GLOBAL_PROOF`.
- What we tried: added adaptive y-subdivision inside L24 while preserving original source-region ownership keys, then updated L31 to validate the actual L24 regular-region count instead of assuming exactly eight. Reran the `CELL-03-03` L24/L31/local-packet chain.
- Why it is reasonable: `CELL-03-03` was not failing length or p-prime exclusion. The only blocker was that one interior source branch produced four regular wall-separation failures and the integration check still assumed the hard-cell regular count. Adaptive y-subdivision is the smallest proof-facing repair because `Fx` was already sign-stable.
- Observed result: L24 closed all `12` regular regions as monotone exclusions using `37` adaptive leaf boxes, max depth used `5`, and zero unresolved adaptive leaves. L31 then closed `12` regular regions plus `52` critical candidates for `64/64` total with no source SHA failures, no ownership duplicates, no source-filter mismatch, and no candidate-count drift. The local packet passes with total validated length upper `18.01426357052187` below cap `20.672796062619668`.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out `CELL-03-03` as a remaining regular-region blocker and rules out the fixed-eight regular-region assumption as a necessary part of L31.
- What it does not rule out: it does not certify the remaining n=14 cells and does not prove global n=14 coverage. It also does not show that all future interior cells will close under the same adaptive depth.
- Next dependency: resume the representative-cell harness, starting with the previously skipped `CELL-07-07`, before launching a broad remaining-cell batch.
- Claim ceiling: local `CELL-03-03` packet only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the correct repair shape: keep the original residual ownership key, show a finite leaf partition closes the regular wall condition, and let the integration certificate count the geometry actually produced.
- Verification commands:
  - `cargo build --release --bin ehp114_n14_regular_residual_decomposition --bin ehp114_n14_residual_chain_integration_cert --bin ehp114_n14_local_hard_cell_certificate_packet`
  - `shasum -a 256 -c *_RESULTS.sha256` for L33K L24, L31, and local packet artifacts.
  - `jq` summaries confirmed `12/12` regular closures, `64/64` residual-chain closure, and local packet pass.
  - forbidden-claim scan over L33K artifacts and touched Rust binaries returned no matches.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33L - Representative-Cell Harness Continuation

- Date: 2026-05-06.
- Method or route: representative-cell harness continuation after variable regular-region repair.
- Artifact ID: `EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-CONTINUATION-20260506-01`.
- Artifact paths:
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-CONTINUATION-20260506-01/EXP-MATH-EHP114-N14-REPRESENTATIVE-CELL-HARNESS-CONTINUATION-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-07-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-07-20260506-01_RESULTS.json`
  - `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-00-20260506-01/EXP-MATH-EHP114-N14-LOCAL-CELL-CERTIFICATE-PACKET-CELL-07-00-20260506-01_RESULTS.json`
- Status: `REPRESENTATIVE_CELL_HARNESS_CONTINUATION_PASS_NOT_GLOBAL_PROOF`.
- What we tried: resumed representative sampling after L33K, running `CELL-07-07` and `CELL-07-00` through slab source generation, source sharpening, high-slope repair where needed, repaired source overlay, critical target splitting, regular residual closure, p-prime exclusions, residual-chain integration, and local packet assembly.
- Why it is reasonable: the first harness stopped at `CELL-03-03` because the regular-region count changed. L33K repaired that class of blocker, so the honest next test was whether the machinery transfers to additional edge/corner cells before broad remaining-cell generation.
- Observed result: both representative cells passed. `CELL-07-07` raw slab total was `21.39738442653137` and repaired to `18.186933338216782` by replacing target `3713:1`. `CELL-07-00` raw slab total was `26.681001315430695` and repaired to `17.13237274486412` by replacing target `2997:0`. Both residual chains closed `64/64` with source hashes passing and no ownership, source-filter, or count drift.
- Discard decision: `retained_as_diagnostic`.
- What it rules out: it rules out `CELL-07-07` and `CELL-07-00` as immediate representative blockers, and it shows the single-target high-slope repair plus variable regular-region integration pattern is not limited to the first non-hard cells.
- What it does not rule out: it does not certify the remaining n=14 cells and does not guarantee that all future cells have one repairable high-slope concentration or simple regular-region closure.
- Next dependency: start L33M, the controlled remaining-cell batch, with stop-on-first-blocker discipline and the same source repair plus variable regular residual closure machinery.
- Claim ceiling: representative-cell continuation only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
- Orthodox reader note: this is the right pre-batch gate: it tests transfer on additional boundary cells and records the repairs as explicit source overlays, not prose exceptions.
- Verification commands:
  - `find /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos114/validated_length -maxdepth 2 \( -path "*CELL-07-07*/*_RESULTS.sha256" -o -path "*CELL-07-00*/*_RESULTS.sha256" -o -path "*REPRESENTATIVE-CELL-HARNESS-CONTINUATION-20260506-01/*_RESULTS.sha256" \) -print0 | while IFS= read -r -d "" file; do dir=$(dirname "$file"); base=$(basename "$file"); (cd "$dir" && shasum -a 256 -c "$base"); done`
  - `jq` summaries confirmed both local packet passes and the harness continuation pass.
  - forbidden-claim scan over L33L, `CELL-07-07`, and `CELL-07-00` artifacts returned no matches.
- SHA status: `PASS_2026-05-06`.
- Overclaim scan status: `PASS_2026-05-06`.

### EHP114-L33M - Active Next Route: Controlled Remaining-Cell Batch Generation

- Date: 2026-05-06.
- Method or route: active next route, controlled remaining-cell batch generation.
- Proposed artifact ID: `EXP-MATH-EHP114-N14-CONTROLLED-REMAINING-CELL-BATCH-20260506-01`.
- Artifact paths: none yet.
- Status: `active_next`.
- What we tried: not yet run. This is the active route exposed after L33L passed on `CELL-07-07` and `CELL-07-00`.
- Why it is reasonable: the pipeline now has passing hard-cell, first non-hard-cell, interior-cell, and additional boundary-cell packets, and the previous software blockers around source contracts, variable regular-region counts, and source overlays have concrete artifacts.
- Observed result: target only. No artifact emitted yet.
- Discard decision: `active_next`.
- What it rules out: nothing yet.
- What it does not rule out: broad n=14 coverage remains open until the remaining-cell batch emits theorem-grade local packets or stops on a named blocker.
- Next dependency: generate local packets for the still-missing n=14 cells under the existing source/subcell contracts. Stop and emit the exact blocker on the first unrepairable source, regular closure, p-prime closure, residual-chain drift, or packet failure.
- Claim ceiling: planned remaining-cell batch only; no claim upgrade.
- Orthodox reader note: the batch must behave like a proof audit, not a production sweep. A failed cell is useful only if it names the exact mathematical or software contract that failed.
- Verification commands: no checksum source yet; run the full artifact contract after target artifact creation.
- SHA status: `NO_SHA_SOURCE`.
- Overclaim scan status: `PASS_2026-05-06`.

## Current Critical Interpretation

The live blocker is not the local length budget. The completed local cell packets now read:

```text
hard cell (6,4) accepted length upper = 20.316451752723314
cell (0,0) accepted length upper      = 15.410486279252329
cell (0,1) accepted length upper      = 15.078165474813495
cell (0,2) accepted length upper      = 17.183471776637536
cell (3,3) accepted length upper      = 18.014263570521870
cell (7,7) accepted length upper      = 18.186933338216782
cell (7,0) accepted length upper      = 17.132372744864120
local cap                             = 20.672796062619668
```

The live blocker is now representative coverage across the remaining missing cells, not local hard-cell length, derivative-root exclusion, theorem-packet assembly, the 64-cell skeleton itself, the existence of a missing-cell worklist, local-pipeline parameterization, direct source-chain parameterization, source absence, the first non-hard source-overage repairs, or the first interior regular-region drift. For the cells that have passed, the L21 regular x-chart regions are closed by monotone Taylor exclusion, direct p-prime exclusion handles part of the critical candidates, L28 root-free disk exclusion handles part of the remaining derivative-root candidates, L29 sharp Taylor/Rouche root-location closes the final derivative-root-near boxes, L31 verifies that those closures compose into a full residual-chain integration certificate, and the local packet binds the result with all source hashes passing.

L33 has now done the atlas audit: the root-affine skeleton enumerates all 64 intended cells, the hard-cell packet `(6,4)` is accepted, and the global gate fails because 63 theorem-grade local cell packets are missing. L33A then checked whether the missing-cell batch could start immediately and found the necessary software blocker: all seven required proof-facing local-packet binaries were hard-coded to `(6,4)`. L33B repaired that boundary: those seven binaries now accept explicit subcell inputs, reject source/subcell mismatches, and no longer contain the fixed-subcell constants flagged by L33A.

L33C then moved one step upstream and parameterized the direct source chain feeding L21: branch isolation, normal collar, third-order collar, and global critical-point exclusion now participate in the same cell/source contract. Its first smoke was intentionally blocked because the `CELL-00-00` slab source did not exist. L33D created that source object and the rerun L33C smoke now passes the source-chain contract.

The new blocker was sharper and more mathematical: the `CELL-00-00` slab source itself was too coarse, with total validated length upper `30.658351028561448` against cap `20.672796062619668` and `4412` unresolved branches. L33E showed this was not diffuse: the top 10 branch contributions sum to `17.77377731502446`, enough to cover the `9.98555496594178` cap excess if those branches could be repaired. L33F then repaired the three largest high-slope branches, reducing their combined length upper from `15.251173022324249` to `0.003308273015128627` and giving adjusted replacement accounting `15.410486279252329 < 20.672796062619668`. L33G materialized that replacement accounting as a source overlay, and L33H ran the downstream local packet to closure: `CELL-00-00` now has a theorem-packet-shaped local artifact.

L33I repeated the controlled process on `CELL-00-01`. The raw slab source was much closer but still over cap: `21.179462302059235` against `20.672796062619668`. The same high-slope repair pattern transferred: repairing targets `2584:3`, `4063:0`, and `2412:1` reduced the adjusted total to `15.078165474813495`, and the downstream residual chain again closed `64/64` with source hashes passing and no ownership/source/count drift.

L33J then ran the representative-cell harness. `CELL-00-02` passed cleanly: the raw slab source was already under cap at `17.183471776637536`, and the residual chain closed `64/64`. `CELL-03-03` exposed the first non-transferable blocker. Its length accounting was repairable, dropping from `22.71308313038541` to `18.01426357052187`, and the p-prime derivative-root chain closed all critical candidates. The failure was the regular lane: L24 saw `12` regular regions, excluded `8`, and left `4` regular wall-separation failures; L31 correctly reported count drift instead of promoting a packet.

L33K repaired that blocker. L24 now uses adaptive y-subdivision while preserving the original regular-region ownership key. On `CELL-03-03`, it closed all `12` regular regions using `37` adaptive leaf boxes with max depth used `5` and zero unresolved leaves. L31 was updated to validate the actual L24 regular-region count, and the chain now closes `12 + 52 = 64` residual buckets with no source SHA failures, ownership duplicates, source-filter mismatch, or candidate-count drift. The `CELL-03-03` local packet now passes under the same local cap.

L33L then resumed representative sampling. `CELL-07-07` and `CELL-07-00` both needed source repair, but each overage was concentrated enough for a single high-slope replacement to bring the local total below cap. Their residual chains both closed `64/64`, and their local packets passed with source hashes clean and no ownership/source/count drift. The active next route is L33M: controlled remaining-cell batch generation under the same contracts, stopping on the first unrepairable blocker rather than smoothing over it.

## Anti-Overclaim Rules

- Do not cite this ledger as a proof of Erdős #114.
- Do not cite this ledger as a global n=14 certificate.
- Do not treat a retired numerical route as a failed mathematical conjecture.
- Do not treat a favorable length budget as a certificate until unresolved branch pieces are closed.
- Do not promote sampled or diagnostic evidence above its artifact status.

## Verification Log

Checksum verification was run on 2026-05-06 for the artifact set cited in lanes EHP114-L01 through EHP114-L14 and EHP114-L16 through EHP114-L33L where SHA sources exist. Each listed checksum returned `OK`.

Overclaim scan: run the project-standard forbidden-claim regex from the implementation plan against both ledger copies, the sidecar notes, the L31, L32, L33A-L33L artifact directories, the patched upstream Rust binaries, and existing L29/RH-MDL artifacts. The literal regex is intentionally not repeated here, because embedding it inside the ledger would make the ledger match its own audit.

This ledger should be regenerated or versioned whenever a new proof-facing artifact is emitted. Do not overwrite this version after external use; create a new dated ledger instead.
