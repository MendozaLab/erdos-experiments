# PR Submission Guide — Erdős #114 to google-deepmind/formal-conjectures
**Date:** 2026-04-10
**Target:** https://github.com/google-deepmind/formal-conjectures
**File:** `FormalConjectures/ErdosProblems/114.lean`
**Status:** 114.lean does NOT exist in repo yet — net-new contribution

---

## Pre-Flight Checklist

### 1. Sign the Google CLA (one-time, blocking)
- Go to https://cla.developers.google.com/
- Sign the Individual CLA for Google projects
- If you've signed for any Google project before, you're already covered
- **The PR bot will block merge until this is done**

### 2. Verify the Lean file compiles locally

Before submitting, you MUST verify the file builds against the repo:

```bash
# Clone your fork
git clone https://github.com/KAmendoza/formal-conjectures.git
cd formal-conjectures

# Install dependencies (downloads Mathlib — takes 10-20 min first time)
lake update
lake exe cache get   # pulls precompiled Mathlib oleans

# Copy your file into place
cp /path/to/114.lean FormalConjectures/ErdosProblems/114.lean

# Build just your file
lake build FormalConjectures.ErdosProblems.«114»

# Or build everything (slower but confirms no conflicts)
lake build
```

**If it doesn't compile:** The most likely issues are:
- `Polynomial ℂ` may need explicit import — check if `ProblemImports` re-exports it
- `MeasureTheory.Measure.hausdorffMeasure` — verify this exists in their Mathlib version
- `‖p.eval z‖` norm notation — may need `open scoped NNNorm` or similar

### 3. Check your file against repo conventions

Your `114.lean` was reviewed against `123.lean` (reference file). Status:

| Convention | Your File | Status |
|------------|-----------|--------|
| Apache 2.0 header | ✅ Present | PASS |
| `import FormalConjectures.Util.ProblemImports` | ✅ Present | PASS |
| Module doc `/-! # Erdős Problem 114 ... -/` | ✅ Present | PASS |
| References section in module doc | ✅ [EHP58], [Po61], [EH99], [Ta25] | PASS |
| `namespace Erdos114 ... end Erdos114` | ✅ Matches pattern | PASS |
| `@[category research open, AMS 30]` on open conjecture | ✅ Present | PASS |
| `@[category research solved, AMS 30]` on small_n | ✅ Present | PASS |
| Definitions have docstrings | ✅ `levelCurveUnit`, `arcLength` | PASS |
| Theorems end with `by sorry` | ✅ Both theorems | PASS |
| No proofs >25 lines included | ✅ Both are `by sorry` | PASS |
| Zenodo DOI reference | ✅ doi:10.5281/zenodo.19480329 | VERIFY |

**One item to verify:** The Zenodo DOI `19480329` in the file — confirm this resolves correctly to the latest EHP preprint version. If you've published a newer version since writing this file, update the DOI.

---

## Submission Steps

### Step 1: Fork the repo
Go to https://github.com/google-deepmind/formal-conjectures and click "Fork"

### Step 2: Create a branch
```bash
git checkout -b erdos-problem-114
```

### Step 3: Add your file
```bash
cp 114.lean FormalConjectures/ErdosProblems/114.lean
git add FormalConjectures/ErdosProblems/114.lean
git commit -m "feat(ErdosProblems): add Problem 114 (Erdős–Herzog–Piranian conjecture)"
```

### Step 4: Push and open PR
```bash
git push origin erdos-problem-114
```
Then go to GitHub and open a PR against `google-deepmind/formal-conjectures:main`.

### Step 5: PR title and body

**Title:**
```
feat(ErdosProblems): add Problem 114 (Erdős–Herzog–Piranian conjecture)
```

**Body:** Use the content from `PR_BODY.md` in this directory. It's already been reviewed and formatted.

### Step 6: Wait for CI + review
- GitHub Actions will run `lake build` on your PR
- A maintainer (likely Moritz Firsching based on Zulip conversation) will review
- Respond to any review comments within 24h
- Expect 1-4 weeks for merge

---

## Zulip Announcement (after PR is open)

Post in `#Formal conjectures` on Lean Zulip (not #lean4):

```
Opened a PR for Erdős Problem 114 (EHP conjecture): [PR link]

Formalizes the statement plus a `solved` tag for n=3–14 via certified interval arithmetic. Separate from #3422.
```

Keep it short. One message. Moritz already knows this is coming.

---

## Cooley IP Filter: PASS

The 114.lean file contains ONLY:
- Standard Mathlib definitions (Polynomial, hausdorffMeasure)
- Public mathematical statements (EHP conjecture, known results)
- Zenodo DOI (public deposit)

It does NOT contain: ErdosAtlas internals, morphism index, MDL scoring, NatureLib, instruments, patent references.

---

## After Merge

1. Update `ERDOS_CREDIBILITY_SCORECARD.md` — add "formal-conjectures PR merged" under #114
2. Reference the PR in the Zenodo deposit description (next version update)
3. Cross-link from erdosproblems.com #114 page (once account access restored)
