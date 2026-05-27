# Erdos #1038 Public Agent Context

Status: `PUBLIC_SAFE_AGENT_CONTEXT`

This file is the single public entry point for Perplexity Computer, Comet,
Claude, Codex, and human reviewers working from the public
`MendozaLab/erdos-experiments` repository.

It is a context packet, not proof evidence. Local packet artifacts, checksums,
Rust/Lean builds, and PR checks decide what is actually verified.

## Public Git Entry

```text
Remote: https://github.com/MendozaLab/erdos-experiments.git
Branch: agent/codex/erdos-experiments-preserve-2026-05-24
Current required head at time of writing: 124c87e
Public context root: Erdos1038/
Agent work root: Erdos1038/agent-work/
Current blocker root: Erdos1038/agent-work/problem-at-hand/
```

External agents must first check the current branch head. If the checked-out
head is older than the current required head, stop with:

```text
STALE_GIT_HEAD_STOP
```

and ask for a refreshed context.

## Problem

Erdos #1038 asks for the extremal measure of:

```text
{x in R : |f(x)| < 1}
```

for monic real-rooted polynomials with all roots in `[-1, 1]`.

The problem remains open. Nothing in this repository claims a solution.

## Route Frame

Use the Everest analogy in `Erdos1038/agent-work/EVEREST_ROUTE_FRAME.md`:

- Everest: a complete accepted proof of #1038.
- Current altitude: strongest local route evidence, not public proof status.
- Current local altitude marker: `8525 m`.
- Public SOTA: unchanged; the problem remains open.
- Route Confound: the missing bridge from fixed-projection/local certificates
  to a globally valid theorem.

Altitude does not move from route design, external opinion, Linear comments, or
synthetic scaffolds. It moves only from verified local packets with receipts.

## What Perplexity Can Read

Start here:

```text
Erdos1038/PERPLEXITY_FULL_CONTEXT.md
Erdos1038/agent-work/README.md
Erdos1038/agent-work/CURRENT_FRONTIER.md
Erdos1038/agent-work/EVEREST_ROUTE_FRAME.md
Erdos1038/agent-work/MODE2_GIT_HANDOFF.md
Erdos1038/agent-work/PUBLIC_CONTEXT.json
Erdos1038/agent-work/SALIENT_FILES_MANIFEST.json
Erdos1038/agent-work/WORK_QUEUE.jsonl
Erdos1038/agent-work/problem-at-hand/README.md
Erdos1038/agent-work/problem-at-hand/PROBLEM_AT_HAND_MANIFEST.json
```

Then inspect the current blocker implementation:

```text
Erdos1038/agent-work/problem-at-hand/backend-source/
Erdos1038/agent-work/problem-at-hand/runner-source/
Erdos1038/agent-work/problem-at-hand/private-route-artifacts/
Erdos1038/agent-work/problem-at-hand/reference-backends/
```

Use `SALIENT_FILES_MANIFEST.json` for the complete public list of files that
matter for #1038 agent work.

## Current Public Work Surface

Current pushed head includes:

```text
EXP-MATH-ERDOS1038-PHI-K-GAP-PERIOD-BASIS-INTERVAL-BACKEND-IMPLEMENTATION-20260527-01
EXP-MATH-ERDOS1038-PHI-K-DEPENDENT-VIETA-IMAGE-CONSUMER-20260527-01
```

Both are scaffold/harness-level public staging packets. They are useful because
they turn route assumptions into fail-closed code and typed packet artifacts.
They do not prove the real theorem.

## Active Next Receipts

The dependent Vieta consumer can only move beyond scaffold if these real
receipts are present and scope-matched:

```text
ROOT_BOX.json
ROOT_MULTIPLICITY_LEDGER.json
ORDERED_ROOT_INTERVALS.json
SCALED_VIETA_IMAGE_CONTRACT.json
FIXED_CLOUD_BOUND_CERTIFICATE.json
ATTAINED_WITNESS_TYPED_DUAL_MARGIN_RESULTS.json
```

If any receipt is missing, emit a missing-receipts ledger. Do not invent
receipt data.

## External Agent Modes

Mode 1: `MODE_1_REVIEW_RESEARCH_OPINE`

Use for literature review, route critique, orthogonal ideas, and falsifiable
next-gate proposals. It may pull public papers. It may not claim a route status
upgrade.

Mode 2: `MODE_2_INTENSE_SOLVE`

Use for concrete work product: patches, complete file contents, proof
scaffolds, theorem statements, algorithms, fixtures, schemas, and verification
commands. Mode 2 work becomes evidence only after local verification or a
GitHub PR check.

End-to-end mode: `MODE_2_INTENSE_SOLVE_END_TO_END`

Use only when explicitly requested. The agent may try to carry #1038 all the
way to a proof candidate by attacking multiple routes, opening PRs, writing
code, drafting Lean/theorem scaffolds, pulling public papers, and proposing
new morphisms or reductions. This is a solve attempt, not just a review.

Even in end-to-end mode, every claimed climb must be backed by receipts and
recorded in:

```text
Erdos1038/agent-work/EVEREST_ROUTE_PLAYBACK.jsonl
```

If a step is only an idea, mark it as claim level `0`. If it is a synthetic
scaffold, mark it as claim level `1`. Do not retroactively summarize route
moves without adding playback rows.

## Etiquette for Perplexity Computer

Allowed:

- clone or read the public repository;
- inspect the current branch head;
- run read-only searches;
- propose patches or open a PR if explicitly requested;
- return complete patches or file contents;
- pull public papers and cite them;
- propose orthogonal attacks, putative morphisms, and falsifiable gates.

Forbidden:

- claim #1038 is solved;
- claim public SOTA movement;
- claim altitude movement;
- claim KKT/global reduction closure;
- claim Lean proof status;
- invent receipts;
- use stale Git context;
- rely on Linear comments as evidence;
- expose credentials, tokens, private repo links, or protected atlas internals.

## PR System

Use GitHub PRs as the public review gate. PRs touching `Erdos1038/` should use:

```text
.github/PULL_REQUEST_TEMPLATE/erdos1038_packet.md
```

Every PR must state:

- packet id;
- claim level;
- changed paths;
- verification commands;
- receipts consumed;
- missing receipts;
- forbidden claims.

CI checks live in:

```text
.github/workflows/erdos1038-agent-work.yml
```

The expected workflow is:

```text
Linear/Perplexity idea -> public branch/PR -> CI and review -> local assimilation -> route map update.
```

Linear remains coordination. Git plus CI plus packet artifacts are the public
artifact truth.

For long solve attempts, one PR may contain multiple commits, but the PR body
must maintain a route ledger summary and point to appended rows in
`EVEREST_ROUTE_PLAYBACK.jsonl`. That makes the climb replayable later for a
human route review or an animation of route findings.
