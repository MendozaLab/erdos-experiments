# EHP114 n=14 Parameterized Local Pipeline Smoke

Experiment: `EXP-MATH-EHP114-N14-PARAMETERIZED-LOCAL-PIPELINE-SMOKE-20260506-01`

## Verdict

- Status: `PARAMETERIZED_LOCAL_PIPELINE_BLOCKED_SOURCE_ARTIFACT_MISSING`
- Smoke subcell: `CELL-00-00`
- Required binaries: `7`
- Hardcoded blocker count: `0`
- Source contract status: `SOURCE_SUBCELL_MISMATCH`
- First failed condition: `per-cell upstream source artifact is missing for CELL-00-00; expected first source like EXP-MATH-EHP114-N14-THIRD-ORDER-COLLAR-REMAINDER-PILOT-CELL-00-00-20260506-01_RESULTS.json`

## Meaning

This is a software-proof rigor gate. It checks that the seven local residual-chain binaries now take explicit subcell arguments and contain no fixed hard-cell constants. It does not certify any missing cell. If the first non-hard-cell source artifact is absent, that is reported as the blocker rather than silently reusing hard-cell data.

## Claim Ceiling

Parameterized local-pipeline smoke gate only. Not a proof of Erdos #114, not a global n=14 certificate, and not an exact lemniscate-length certificate.
