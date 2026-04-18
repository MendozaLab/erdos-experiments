# Track A Author Approval Review
**Artifact:** Track A release packet for Erdős #30 spectral Sidon note  
**Status:** NOT APPROVED FOR DOI RELEASE AS CURRENTLY WRITTEN  
**Purpose:** Author-facing review memo for approval of revisions prior to any GitHub Release or Zenodo DOI minting.

## Files Reviewed

1. `Math/erdos-experiments/Erdos30/paper/sidon_choreography.tex`
2. `Math/Shadow_Numberverse/publication-package/ADVERSARIAL_CLAIM_MATRIX.md`
3. `Math/Math-Problems/Erdos-Standard/erdos-30/sidon_crystal_engine/Erdos30_SidonCrystal_Skeleton.lean`
4. `Math/Math-Problems/numberverse/Erdos-Standard/erdos-30/EXP-TDP-SIDON-003_RESULTS.json`
5. Supporting evidence checked:
   `sidon-choreography/SpectralSidon.lean`
   `sidon-choreography/EXP-TDP-SIDON-001_RESULTS.json`

## Executive Decision

Track A is potentially releasable as a **bounded computational note**.

Track A is **not** currently releasable as:

- a zero-sorry Lean theorem package
- a formal proof of the spectral Sidon claim
- an asymptotic advance on Erdős Problem #30
- a release claiming any resolution, barrier break, or proof-grade bridge from the spectral model to the classical conjecture

## Approval Condition

Author approval is recommended **only if** the packet is revised so that every public artifact uses the same honest scope:

- `Computationally verified finite-range result`
- `Proof architecture for the Lean bridge`
- `No claimed asymptotic improvement for Erdős #30`

## Required Revisions Before Approval

### Priority 1: Scope Correction

- Downgrade the main Track A claim from theorem/proof language to finite computational verification language.
- Remove all references to "Theorem-Grade," "formally verified," "compiler verified," and "zero sorries" for the spectral Sidon theorem itself.
- State explicitly that the formal Lean support currently covers the classical integer Sidon infrastructure, not a completed Lean proof of the spectral theorem.

### Priority 2: Remove Overclaiming

- Remove any statement that the work "bypasses the Erdős-Turán `O(N^{1/4})` limit."
- Remove any statement that the work "establishes full geometric resolution of Erdos 30."
- Remove any statement that the `N^{1/4}` barrier is broken.
- Remove any implication that the spectral computations constitute an asymptotic proof for the classical Sidon conjecture.

### Priority 3: Repair the Lean Narrative

- Reclassify `Erdos30_SidonCrystal_Skeleton.lean` as a proof skeleton only.
- Reclassify `sidon-choreography/SpectralSidon.lean` as proof architecture only.
- Do not cite either file as zero-sorry or proof-complete unless that state is actually achieved and re-verified.

### Priority 4: Repair the Bridge Claim

- Reframe the use of the compiled integer Sidon theorem as motivation or an analogue.
- Do not present the discretization/rescaling bridge from `D_n ⊂ R` to the integer Lean theorem as a proved corollary unless the bridge lemma is explicitly written and verified.

### Priority 5: Disclosure and Release Hygiene

- Add the required AI acknowledgment block to the manuscript before any public release.
- Ensure any public README or release description uses the same bounded scope.
- Exclude any artifact whose wording still overclaims the result.

## Claims That Can Safely Remain

The following claims are supportable by the currently reviewed evidence:

- The spectral distance set
  `D_n = {2 sin(pi m / n) : 1 <= m <= floor(n/2)}`
  is the object under study.
- Exhaustive finite verification supports that `D_n` is Sidon for `3 <= n <= 29`.
- There is an explicit collision at `n = 30`.
- In the tested range `3 <= n <= 200`, the observed collision values are
  `{30, 42, 60, 84, 90, 120, 126, 150, 168, 180}`.
- In that tested range, all observed collision values are divisible by `6`, while divisibility by `6` is not sufficient.

## Claims That Must Not Appear In Any Release

- "Zero-sorry Lean proof of the spectral Sidon theorem"
- "Formally verified Track A theorem"
- "Theorem-grade Lean proof" for the spectral claim
- "Barrier broken" language for the classical `N^(1/4)` term
- "Resolution of Erdős #30"
- Any suggestion that finite computation here proves the asymptotic conjecture

## Recommended Author-Facing Positioning

If approved after revision, Track A should be described as:

`A finite computational study of spectral distance sets associated with Z_n choreographies, together with a proof-architecture bridge to classical Sidon infrastructure in Lean 4.`

That framing is strong, honest, and citable without exposing the packet to an avoidable credibility failure.

## Suggested Approval Statement

Use one of the following:

### Approve for Revision and Release Preparation

`Approved to proceed only after the required scope corrections and disclosure fixes in this review are implemented. This is approved as a bounded computational Track A release, not as a theorem-grade Lean release.`

### Approve for Internal Circulation Only

`Approved for internal use and drafting only. Not approved for GitHub Release, Zenodo DOI minting, or public distribution in current form.`

### Reject Current Release Packet

`Not approved. Current packet materially overclaims the formal status and mathematical consequences of the Track A artifacts.`

## Author Sign-Off

**Author decision:**  
`[ ] Approve for revision and release preparation`  
`[ ] Approve for internal circulation only`  
`[ ] Reject current release packet`

**Scope accepted by author:**  
`[ ] Computational finite-range note`  
`[ ] Proof architecture only for spectral Lean layer`  
`[ ] No asymptotic claim on Erdős #30`

**Author name:** ____________________  
**Date:** ____________________  
**Notes:** ____________________
