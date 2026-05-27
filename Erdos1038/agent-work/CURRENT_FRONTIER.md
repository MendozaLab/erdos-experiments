# Erdős #1038 Current Frontier

Status: `PUBLIC_SAFE_FRONTIER_SUMMARY`

This is a public-safe summary for external agents. It intentionally omits
private paths, private repository links, Linear URLs, credentials, protected
scoring internals, and unpublished atlas internals.

## Problem

Erdős #1038 asks for the extremal measure of the real sublevel set
`{x : |f(x)| < 1}` for monic real-rooted polynomials with all roots in
`[-1, 1]`.

The problem remains open. This staging lane does not claim a solution.

## Current Route Shape

The active route is studying a fixed public-candidate cloud for a large-degree
root-factor projection. In that projection, the sublevel boundary has a
finite-component structure. Current local work has isolated a period/harmonic
measure obstruction that must be made legitimate before any residual audit can
mean anything.

Use the Everest route analogy from `EVEREST_ROUTE_FRAME.md`: the current
altitude is local route evidence, not public proof status. The current local
altitude remains `8525 m`, and Everest remains the full accepted proof.

The current sanitized invariants are:

- cloud component count: `25`;
- gap-period row count: `24`;
- normalization is separate and is not a period row;
- one boundary-near support atom is routed as endpoint-limit primary;
- the right-exterior route is not certified for that atom;
- weighted-QR basis conditioning is only an f64 target until interval-certified;
- the missing backend is a directed-interval basis certificate.

## Current Blocker

The useful next blocker is not a residual audit. The residual audit is blocked
until the period object, endpoint-source admissibility, and basis interval
certificate are all sound.

Primary current packet landed locally and is now mirrored here:

```text
EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-INTERVAL-BACKEND-IMPLEMENTATION-20260527-01
```

Status:

```text
GAP_PERIOD_BASIS_INTERVAL_BACKEND_IMPLEMENTATION_PASS__FAIL_CLOSED_FIXTURES_PASS__PRIVATE_NUMERIC_PAYLOADS_PENDING
```

Meaning: the Rust/Inari backend harness now exists and passes synthetic
fail-closed fixtures. It rejects the known bad cases before any real private
numeric payload is consumed:

- 27-vs-24 row mismatch;
- normalization row leakage;
- sign-convention mismatch;
- transform condition above threshold;
- missing recovered-direction independence witnesses;
- nonpositive smallest singular value.

It does not certify the real B1/B2/B3 payload. The next packet target is:

```text
EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-PRIVATE-PAYLOAD-INTEGRATION-20260527-01
```

Goal: feed real private numeric interval payloads into the backend for:

- B1: transform interval certificate;
- B2: independence witnesses for the recovered basis directions;
- B3: transformed period matrix interval entries and singular-value/condition
  bounds.

The weighted-QR route is structurally suspect: f64 transformed matrix condition
is about `264.93593439501007`, but the f64 transform condition is about
`4.726362562410606e16` and the seed evaluation rank is only `13`. Treat this as
high-risk until B1/B2 certify. Prepare the hyperelliptic-canonical basis
fallback in parallel.

Parallel theorem target:

```text
EXP-MATH-ERDOS1038-PHI-K-ENDPOINT-LIMIT-SOURCE-KERNEL-GATE-20260527-01
```

Goal: prove or falsify that the boundary-near endpoint-limit source has a
finite row-compatible period contribution and does not leak into normalization
or strict interior-gap source theorems.

## Claim Ceiling

Allowed from this directory:

- public-safe route critique;
- literature pointers;
- theorem scaffolds;
- patch proposals;
- numerical-backend designs;
- falsifiable next local gates.

Forbidden from this directory:

- claiming #1038 is solved;
- claiming public SOTA improvement;
- claiming altitude upgrade;
- claiming period legitimacy;
- claiming residual audit closure;
- claiming KKT closure;
- claiming global reduction;
- claiming Lean proof;
- claiming independent coefficient-box theorem.
