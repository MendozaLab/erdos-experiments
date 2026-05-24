# EHP114 Reception-Gated Next Steps

Created: 2026-05-08

Purpose: prepare the next moves for Erdős Problem #114 without overclaiming. We went deep into a complete-proof direction, then deliberately pulled back to the strongest credible public claim: a finite certificate-backed statement for `1 <= n <= 14`, with the all-degree conjecture still open. This document keeps the stronger material organized, but gates any escalation on reception from erdosproblems.com, Google DeepMind formal-conjectures, or a serious mathematical reviewer.

## Current public state

- Zenodo corrected certificate packet is live:
  - version DOI: `10.5281/zenodo.20087919`
  - concept DOI: `10.5281/zenodo.19184467`
  - version: `v3.1.0`
  - key correction: `n = 13` now has `bb_total_evals = 197132288`, not the old zero-eval artifact.
- erdosproblems.com post:
  - URL: `https://www.erdosproblems.com/forum/thread/114#post-6330`
  - status checked 2026-05-08: visible to logged-in account and still `Awaiting moderator approval`.
  - canonical problem page still says there are no partial or complete solutions claimed in comments.
  - old March `n = 3..12` post was accepted as a forum comment, not incorporated into the canonical problem remarks.
- Google DeepMind formal-conjectures:
  - PR #3958: `https://github.com/google-deepmind/formal-conjectures/pull/3958`
  - state: open, review required.
  - head: `713d3fa`.
  - visible checks: `labeler` and `cla/google` green only. No full upstream Build check is visible, so do not say CI green.
  - Lean file keeps all-degree `erdos_114` open and finite sibling `erdos_114_finite_le_14` as externally certificate-backed with `sorry`.

## Public claim ladder

Use the lowest rung that matches the reception.

1. **Forum note / finite certificate**  
   Safe wording: "finite certificate-backed result for `1 <= n <= 14`; full conjecture remains open."  
   This is the current public posture.

2. **Finite-side preprint**  
   Trigger: erdosproblems moderator or a serious reviewer says the finite side should be submitted as a note/preprint.  
   Safe wording: "dependency-honest finite-degree verification architecture."  
   Not safe: "proof of Erdős #114" or "complete solution."

3. **Reviewer packet / computation audit**  
   Trigger: someone questions the computation, `n = 13`, reproducibility, or DOI contents.  
   Safe wording: "Here is the source release, version DOI, SHA sidecars, and rerun audit."  
   Do not defend by rhetorical confidence. Point to artifacts.

4. **Full-completion reserve lane**  
   Trigger: a credible reviewer asks whether the finite certificate can be combined with Tao's sufficiently-large-`n` theorem.  
   Safe wording: "There is a structured bridge program, but the Tao threshold remains opaque and no all-degree result is claimed."  
   This is not public-ready as a theorem.

## Reception decision tree

### A. erdosproblems approves post without comment

Action:
- Do not immediately add another comment.
- Record that the post is approved.
- Wait for any moderator/reviewer response before escalating.
- Then optionally prepare one short pointer post to the approved thread and Zenodo DOI.

Do not:
- Claim the canonical problem page accepted the result.
- Claim a partial solution was incorporated.
- Use "Erdos Problems approved the proof."

### B. erdosproblems asks for a preprint or shorter note

Action:
- Use `Math/preprints/ehp114-finite/` as the starting manuscript.
- Update DOI references from concept-only wording to:
  - version DOI `10.5281/zenodo.20087919`
  - concept DOI `10.5281/zenodo.19184467`
- Recompile the paper.
- Run publisher + crackpot-scrub + prepub-redteam on the exact PDF/post text.
- Deposit the preprint only after those gates pass.

Prepared response:

```text
Thanks. I will prepare this as a short finite-side note rather than as a claim on the all-degree problem. The intended scope is exactly `1 <= n <= 14`, with the all-degree conjecture left open and Tao's sufficiently-large-n result cited only as complementary context.
```

### C. erdosproblems challenges computational status

Action:
- Answer with the corrected packet, not a new claim.
- Lead with the `n = 13` correction and SHA.
- Offer the source release and JSON path.

Prepared response:

```text
The corrected version is pinned at version DOI 10.5281/zenodo.20087919. The `n = 13` row in pre-v3.1.0 material should be treated as superseded; the old symptom was `bb_total_evals = 0` with an empty `bb_levels` array. The corrected v3.1.0 canonical row has SHA-256 `c06c633b4053cdf2c4c6003327f30ee2ce683e6a9da70e047d5c5d43e685fd17` and `bb_total_evals = 197132288`.
```

### D. erdosproblems asks whether this is a solution

Action:
- Say no, not the all-degree problem.
- Keep the finite result useful but narrow.

Prepared response:

```text
No. I am not claiming a solution of the all-degree Erdős #114 conjecture. The claim is a finite, certificate-backed statement for `1 <= n <= 14`. Tao's theorem handles all sufficiently large `n`, but I am not extracting a numerical threshold from that theorem here, and I am not claiming anything for `n >= 15`.
```

### E. DeepMind formal-conjectures engages on PR #3958

Likely reviewer concerns:
- Does `@[category research solved]` fit a theorem with `sorry`?
- Should the finite sibling be included if the certificate logic is external?
- Should the DOI be concept DOI only or version DOI?
- Should MacLane be cited in the Lean docstring?

Response posture:
- Agree quickly on style changes.
- Preserve the proof-integrity boundary: no certificate axioms, no fake closure.
- If asked, downgrade category or wording rather than defend "solved" too hard.
- Keep all-degree `erdos_114` open.

Prepared response:

```text
Happy to adjust the category/wording to match repository convention. The important design choice is that the external certificate is not imported as an axiom or disguised as a Lean proof; the finite statement remains a `sorry` until the certificate logic is formalized in Lean.
```

### F. DeepMind remains silent

Action:
- Do not ping again immediately.
- The next meaningful update should be a maintainer-facing response, not a status nudge.
- If erdosproblems approves the forum post or asks for a preprint, then a short PR comment can mention the public forum/preprint state.

## Full-completion reserve: what is real and what is not

The deeper program has two independent proof lanes. They should not be publicly blended until both are clean.

### Lane 1: finite/frontier proof architecture

Current public surface:
- `1 <= n <= 14` finite certificate packet on Zenodo v3.1.0.
- Formal-conjectures #3958 records the statement but does not prove it in Lean.

Deeper internal material:
- `IN_SEARCH_OF_EHP114_LARGE_N_SMALL_N_STORY_2026-05-07.md` says an orthodox n=14 local-atlas proof lane had 48/64 local cell packets passing, with 16 cells still open at that point.
- The local-atlas lane is valuable, but it is not the same thing as the public DOI certificate packet.

Rule:
- Do not use the local-atlas narrative as public proof evidence unless the current cell accounting, source hashes, and global atlas integration are regenerated and frozen.

### Lane 2: Tao threshold effectivization

Current status from local packets:
- `EXP-MATH-EHP114-TAO-THRESHOLD-EXTRACTION-SKELETON-20260507-02`: dependency table ready, no numeric threshold.
- `EXP-MATH-EHP114-TAO-EFFECTIVE-THRESHOLD-CHECKER-20260507-02`: `TAO_THRESHOLD_REMAINS_OPAQUE`, 44 dependency rows and 44 opaque rows, candidate `N0 = null`.
- `EXP-MATH-EHP114-TAO-FINAL-SECTION-THRESHOLD-SUBCHECK-20260507-01`: `FINAL_SECTION_REMAINS_OPAQUE`, no global candidate `N0`, authorizes no higher-degree computation.

Rule:
- Do not claim a finite middle range.
- Do not run or cite `n >= 15` as authorized by Tao until a numeric threshold is extracted or a separate direct certificate is generated.

## What to prepare now

### Immediate

- Keep PR #3958 quiet unless a reviewer comments.
- Monitor erdosproblems moderation for post #6330.
- Patch the finite preprint DOI references to `10.5281/zenodo.20087919` before any deposit.

### If reception is positive

- Build `EHP114 finite note v3.1.0`:
  - one short preprint PDF,
  - Zenodo version DOI in abstract/data statement,
  - source release link,
  - n=13 audit paragraph,
  - all-degree caveat.
- Prepare a one-paragraph update for formal-conjectures linking the accepted forum/preprint record.
- Prepare short pointer posts only after the forum post is approved or preprint is live.

### If reception is skeptical

- Send a compact audit packet:
  - version DOI `10.5281/zenodo.20087919`,
  - GitHub release `v3.1.0`,
  - corrected n=13 SHA,
  - rerun log,
  - old zero-eval artifact archived as superseded,
  - exact claim ceiling.
- Do not introduce MDL, morphisms, physics, or platform framing in the reply.

### If reception asks about full completion

- Say the complete-proof program exists but is not yet a theorem.
- Offer the bridge map:
  - finite side: public `1 <= n <= 14` certificate,
  - large side: Tao sufficiently-large-`n`,
  - missing side: effective threshold or direct finite middle certificates.
- Do not publish the "complete proof" story until Tao threshold opacity is resolved.

## Safe social pointer after approval

Use only after post #6330 is approved or a preprint is live.

```text
I posted a finite-side update for Erdős #114: a corrected, SHA-pinned certificate packet for `1 <= n <= 14`, with the full all-degree conjecture still open. The corrected v3.1.0 packet is at version DOI 10.5281/zenodo.20087919; concept DOI 10.5281/zenodo.19184467 resolves to latest.
```

Avoid:

```text
Solved Erdős #114.
DeepMind accepted it.
Erdosproblems approved the proof.
Tao + our computation closes the conjecture.
```

## Next agent prompt

```text
You are working on MendozaLab EHP114 reception handling. Do not make public claims beyond the finite certificate-backed statement `1 <= n <= 14`. Check:
1. erdosproblems post #6330 moderation state,
2. formal-conjectures PR #3958 state/checks/comments,
3. Zenodo latest version DOI still `10.5281/zenodo.20087919`.

If post #6330 is approved, prepare a short status note and update the finite preprint DOI references. If a reviewer challenges the result, answer with artifacts only: Zenodo v3.1.0, GitHub release, n=13 SHA, rerun log, and open-conjecture caveat. Do not claim a full solution unless a separate audit shows Tao threshold extraction and finite middle coverage are complete.
```

