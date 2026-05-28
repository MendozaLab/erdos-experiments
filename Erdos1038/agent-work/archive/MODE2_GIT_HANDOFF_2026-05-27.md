# Perplexity Computer Mode 2 Handoff

Status: `READY_FOR_MODE_2_INTENSE_SOLVE_PUBLIC_GIT`

Use this file as the public Git handoff for a Perplexity Computer Mode 2 round.

## Required Output Prefix

```text
Perplexity Computer mode: MODE_2_INTENSE_SOLVE
Git access check: VERIFIED_PUBLIC_GIT | UNABLE_TO_VERIFY_GIT_ACCESS
Freshness check: CURRENT | STALE_CONTEXT_STOP | UNABLE_TO_VERIFY
Status: WORK_PRODUCT_INTENT_ONLY
Assimilation: NOT_ASSIMILATED
```

## Public Git Access

Primary public repo:

```text
Remote: https://github.com/MendozaLab/erdos-experiments.git
Relevant subpath: Erdos1038/agent-work/
```

Current problem-at-hand bundle:

```text
Erdos1038/agent-work/problem-at-hand/
```

Related public context repos:

```text
https://github.com/MendozaLab/math-morphism-atlas
https://github.com/MendozaLab/mathlib-prs
```

Allowed Git actions:

```text
clone
fetch
checkout/read
status
log
show
grep
diff
```

Forbidden Git actions:

```text
push
commit
reset --hard
clean -fd
checkout -- paths
rebase
merge
remote branch creation
open pull request
close issue
publish
delete files
mutate repository state
```

## Freshness Guard

Before solving, inspect this directory and check whether `CURRENT_FRONTIER.md`,
`PUBLIC_CONTEXT.json`, or `WORK_QUEUE.jsonl` points to a newer frontier than the
one in this handoff. If there is a newer frontier, stop with
`STALE_CONTEXT_STOP` and report it. If public Git is unavailable, stop with
`UNABLE_TO_VERIFY_GIT_ACCESS` and do not invent file contents.

## Task

Return concrete work product for:

```text
EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-INTERVAL-BACKEND-IMPLEMENTATION-20260527-01
```

Design a public-safe backend plan for the missing directed-interval Stage-B
basis certificate.

Use `problem-at-hand/` as the primary source context. It includes the current
runner source, relevant sanitized packet receipts, Cargo surface, and reference
Rust/Inari backend patterns.

Also read `EVEREST_ROUTE_FRAME.md` before answering. Use its terms exactly:
Everest, current altitude, Route Confound, and next ridge. Do not upgrade
altitude or claim summit progress from Mode 2 work product.

Required sections:

1. Backend architecture.
2. Proposed Rust/Python file layout.
3. Output JSON schema.
4. B1 transform interval certificate.
5. B2 recovered-direction independence witnesses.
6. B3 transformed matrix certificate.
7. Verification commands.
8. Failure modes.
9. Orthogonal fallback if weighted QR is structurally unsound.
10. Claim ceiling.

## Public-Safe Boundaries

Use public papers and public repositories freely. You may use ErdosAtlas and the
Erdős Collider / PMF framing only as clue generators, not as proof evidence.
Distinguish literature-established facts from putative morphisms.

Do not request or expose credentials, private repo URLs, private Linear URLs,
local filesystem paths, protected scoring algorithms, private transfer
operators, dimensional decompositions, or unpublished atlas internals.

Do not claim local packets landed or tests passed. Return proposed work only.
