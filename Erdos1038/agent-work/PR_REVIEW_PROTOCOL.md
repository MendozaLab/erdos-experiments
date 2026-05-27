# Erdos #1038 PR Review Protocol

Status: `PUBLIC_SAFE_PR_PROTOCOL`

This protocol gives #1038 a Formal-Conjectures-style public review queue while
preserving the route's evidence discipline.

## Identity

Every #1038 PR must name:

```text
Problem: Erdos #1038
Packet ID: EXP-MATH-...
Claim level: 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7
Scope: scaffold | local certificate | dependent-image theorem | ...
```

## Claim Levels

```text
0 = external work product / intent only
1 = local synthetic scaffold or harness pass
2 = real local certificate receipts consumed
3 = local theorem composition over fixed projection
4 = local theorem with legitimacy/attainment composition
5 = coefficient-facing theorem in explicitly named scope
6 = public SOTA comparison verified
7 = full #1038 theorem candidate
```

Do not skip levels. If receipts are missing, say exactly which receipts are
missing.

## Required PR Checks

At minimum, a PR touching `Erdos1038/` should pass:

```text
python runner compile
JSON parse checks
RESULTS.sha256 checks
cargo check for backend-source
claim-ceiling scan by reviewer
```

Lean claims require a current `lake build` receipt. CLEAN is not COMPILED.

## Required PR Body Sections

Use `.github/PULL_REQUEST_TEMPLATE/erdos1038_packet.md`.

The PR body must include:

- summary;
- packet id;
- claim level and claim ceiling;
- changed paths;
- receipts consumed;
- receipts missing;
- verification performed;
- external-review provenance, if any;
- forbidden claims.

## Review Rules

Approve only if:

- the packet id is immutable and versioned;
- generated artifacts are checksumed;
- numeric results carry a scope or an explicit scaffold-only ceiling;
- superseded packets have sidecars if applicable;
- no Linear comment is treated as evidence;
- no public theorem/SOTA/altitude claim outruns the receipts.

Request changes if:

- the PR claims a theorem from f64 or synthetic fixtures;
- the row count, dimension, atom count, or sign convention is implicit;
- missing receipts are not listed;
- a result leaks from fixed-projection scope into coefficient-box or global
  scope without a named lift packet;
- the PR relies on stale Git head context.

## Perplexity / Claude Use

External agents may submit PRs or patch proposals when explicitly asked. Their
output remains `WORK_PRODUCT_INTENT_ONLY` until CI passes and a local
assimilation packet accepts it.

## End-to-End Solve PRs

If an external agent is explicitly asked to take #1038 as far as possible, it
may open an end-to-end solve PR rather than one packet PR. In that case:

- use semantic mode `MODE_2_INTENSE_SOLVE_END_TO_END`;
- keep commits small enough that each route move is reviewable;
- append a row to `EVEREST_ROUTE_PLAYBACK.jsonl` for each route event;
- keep the PR body updated with the current route frontier;
- do not squash away intermediate route evidence until after review;
- do not claim the summit unless the PR includes a complete proof candidate
  and all required verification receipts.

End-to-end PRs still use the same claim ladder. A long PR can contain many
level-0/1/2 moves without being a theorem PR.
