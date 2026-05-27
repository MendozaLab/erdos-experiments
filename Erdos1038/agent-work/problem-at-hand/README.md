# Problem at Hand: Gap-Period Basis Interval Backend

Status: `PUBLIC_SAFE_PROBLEM_AT_HAND_BUNDLE`

This bundle is intentionally narrow. It gives Perplexity Computer enough public
Git context to work on the current #1038 blocker without needing private repo
access.

## Target Packet

```text
EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-INTERVAL-BACKEND-IMPLEMENTATION-20260527-01
```

## Goal

Design the missing Rust/Inari backend for the Stage-B stable gap-period basis
certificate.

The existing local Stage-B packet says:

```text
STABLE_GAP_PERIOD_BASIS_INTERVAL_CERTIFICATE_PARTIAL__STRUCTURAL_PRE_AUDIT_PASS__STAGE_B_BLOCKED_BACKEND_REQUIRED
```

The required backend stages are:

- B1: interval certificate for the weighted-QR transform;
- B2: interval independence witnesses for the 11 rank-recovered directions;
- B3: transformed 24x24 period-matrix interval entries plus singular-value or
  condition-number bounds.

## Directory Layout

```text
private-route-artifacts/
  Sanitized packet receipts relevant to this blocker.

runner-source/
  Python runners for the Stage-A conditioning gate and Stage-B backend
  requirement packet.

backend-source/
  Cargo files and current Rust library surface.

reference-backends/
  Existing Rust/Inari backend patterns from related #1038 packets.
```

## What To Return

Return `WORK_PRODUCT_INTENT_ONLY`, not a claim that anything landed.

Required sections:

1. Backend architecture.
2. Patch proposal with files to add/edit.
3. Output JSON schema.
4. B1 transform interval certificate.
5. B2 recovered-direction independence certificate.
6. B3 transformed matrix certificate.
7. Verification commands.
8. Failure modes.
9. Orthogonal fallback if weighted QR cannot be certified.
10. Claim ceiling.

## Claim Ceiling

This bundle does not solve #1038, does not improve public SOTA, does not upgrade
altitude, does not prove period legitimacy, does not run a residual audit, does
not close KKT or global reduction, and does not prove a Lean theorem.
