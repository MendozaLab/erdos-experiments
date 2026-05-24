# EXP-MATH-EHP114-N14-BRANCH-COUNT-CENSUS-HARD-CELL-20260506-01 Report

## Verdict

- Status: `BRANCH_CENSUS_FAIL_AMBIGUOUS_ROOT_ISOLATION`
- Source: `EXP-MATH-EHP114-N14-SLAB-VALIDATED-LENGTH-HARD-CELL-20260506-03`
- Fine x-slabs: `7040`
- Active fine x-slabs: `2480`
- Active slabs with exactly one accepted branch and no unresolved tube: `0`
- Fine slabs with unresolved tubes: `2083`
- Fine slabs with multiple branch groups: `2441`
- Coarse buckets: `100`

## Distribution

Fine-slab total branch/tube count distribution:

```json
{
  "1": 39,
  "2": 1585,
  "3": 228,
  "4": 397,
  "5": 54,
  "6": 177
}
```

## Meaning

The z32 slab artifact already puts the certified accepted branch length below
the hard-cell cap, but this census shows the single-branch assumption is not
yet certified across the domain. The next proof-facing step is therefore not
another global length pass. It is an interval Newton or Bernstein isolation
resolver over the unresolved tubes named by this census.

## Claim Ceiling

This is a local branch-count diagnostic. It is not a proof of Erdős #114 and
not a global n=14 proof.
