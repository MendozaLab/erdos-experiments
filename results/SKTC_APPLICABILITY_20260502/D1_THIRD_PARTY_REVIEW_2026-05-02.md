# D1 Third-Party Consistency Review — SKTC Applicability Screen

**Reviewer role:** External technical reviewer, READ-ONLY posture.
**Date:** 2026-05-02
**Subject artifact:** `EXP-MATH-ERDOS-SKTC-APPLICABILITY-SCREEN-20260502-01`
**Mandated claim boundary (verbatim, enforced throughout):** *"We have an internal routing screen suggesting that the tensor/operator/cone-certificate strategy developed around EHP114 may be worth testing across three families of Erdős problems: polynomial-potential geometry, constrained-state additive geometry, and graph-spectral/discrepancy geometry. This is not a proof, not a solved-problem claim, and not a public morphism verdict."*

## Executive verdict

**`CONSISTENT_WITH_LIMITS`** with one **MATERIAL DRIFT** flag against any *public* downstream use of the Tier-1 list.

D1 was unreachable in this review environment, but local-DB consistency was substantially verifiable. The screen is internally consistent with the local `erdosatlas.db` source it actually reads (solvent index, physics-edge similarity, math-density / math-structure proxies, hand-anchored boost). However the screen **does not read `problems.math_status`**, and four of its top fifteen "TIER_1_PROTOTYPE_NOW" rows are already `solved` and one is `disproven` in that same local DB. For internal routing this is fine — the strategy itself can still be tested against a closed problem as a calibration / negative-control move. For any *public* claim, scorecard upgrade, morphism-verdict assertion, or downstream evidence-binding step the rows must be re-labeled before they leave the screen.

Hashes match: results JSON `10f3463de2480ee792866061c1ee60940486085bb792df086d6bde08664040ba`, report `843ea08e538f854c1079695299c48e619b917179d5a3e16c80db2f823bf05e65` — both reproduce exactly. No artifact tampering detected.

## D1 access summary

D1 (`research-hub-auth`, account `def73806d72730c63c3f19d95d1653f8`) was attempted and confirmed unreachable in this environment.

```
Attempt 1: wrangler d1 execute research-hub-auth --remote --command "SELECT name FROM sqlite_master WHERE type='table'"
  → Wrangler requires Node.js >= v20 (system node was v16.20.2)

Attempt 2: nvm use 20 + wrangler d1 execute …
  → wrangler 4.22.0 returned:
     "In a non-interactive environment, it's necessary to set a CLOUDFLARE_API_TOKEN
      environment variable for wrangler to work."
  → CLOUDFLARE_API_TOKEN not present in environment; CLOUDFLARE_ACCOUNT_ID alone is insufficient.
```

No alternate D1 path was available (no wrangler login session, no CF API token, no MCP D1 binding live in this thread). Per the task's own rule, the verdict therefore cannot be `CONSISTENT` outright — it must absorb that limit. Read-only confirmation: zero write attempts were made against any DB (D1, local SQLite, scoring DB) over the course of this review.

## Local source reality (what the screen actually consumes)

`erdos_problems.db` at `/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos_problems.db` is **0 bytes** (mtime 2026-04-11). It is dead and cannot be canonical for anything. Any session that thought it was reading from there is reading from an empty file.

The live local store is `erdosatlas-workbench/erdosatlas.db` (17.1 MB, mtime 2026-04-18). Counts at review time:

| table | rows |
|---|---:|
| `problems` | 5,496 |
| `morphisms` | 54,203 |
| `physics_objects` | 34 |
| `physics_edges` | 14,105 |
| `solvent_index_v2` | 1,190 |
| `spectral_operators` | **0** |

Two observations matter for the screen.

First, the screen's `rows_screened` is `1201`, which roughly matches `solvent_index_v2` (1,190) plus a small number of joined rows — i.e., the effective universe is not the 5,496 problems in `problems`, it is the ~1,200 problems that have solvent-index rows. Anything outside `solvent_index_v2` is silently absent from the ranking. That is a defensible scoping decision but worth recording.

Second — and this is the load-bearing local finding — the `spectral_operators` table is **empty**. The 7-feature applicability signature in the report header (extremal functional, finite state space, symmetry, boundary/singular layer, operator/tensor representation, calibration ladder, negative controls) reads like a check against operator-class data. It is not. Nothing in the scoring pipeline (`scripts/score_sktc_applicability.py`) consumes `spectral_operators` or `spectral_operator_edges`. The score is a weighted sum of `phys_sim`, `math_density`, `si_score` (solvent), `math_str`, `consensus`, a domain prior, and a hand-anchored boost on 11 curated problems. The feature signature in the report header is **anchor prose**, not column-checked feature credit.

## Candidate audit table

`db_status` is from `erdosatlas.db :: problems.math_status`. `db_evidence` is local solvent index rank + count of physics edges + max physics-edge similarity (where available). `morphism/proof` is "anchor prose" if the row was hand-anchored, otherwise "metadata-only." `alignment_verdict` is the row-by-row comparison against what the screen advertises.

| problem_id | local_tier | db_status | db_evidence (local) | morphism/proof | alignment_verdict | next_gate |
|---:|---|---|---|---|---|---|
| **114** | TIER_0_CALIBRATION | open | SI 0.9807 (rank uncomputed via SI; phys_sim 0.9463, 22 phys edges) | anchor: lemniscate length deficit, polynomial-potential geometry; **NO Lean COMPILED at this path** (local file is `Math-Problems/proofs/archive/EHP114_closed.lean` and DeepMind formal-conjectures `1145/1148/1141/1142.lean`) | **OK as calibration anchor** but the report's headline conflates EHP114 small-n compute with a "calibrated" certification; per Math/CLAUDE.md and memory.md the canonical status is "EHP114 = NO PROOF; small-n compute only" | statement-and-literature audit + clarify calibration scope |
| **755** | TIER_1_PROTOTYPE_NOW | **solved** (local DB) | SI 0.984 rank 11; phys_sim 0.7565, 6 phys edges | anchor: B_h[g], constrained-state geometry | **DRIFT — solved problem ranked #1.** Methods test still legitimate as negative control, but the row must not be cited as an open-problem opportunity | confirm published resolution via D1 + `notes` field; relabel as calibration / methods-only |
| **166** | TIER_1_PROTOTYPE_NOW | **solved** (local DB) | SI 0.991 rank 7; phys_sim 0.8707, 9 phys edges | anchor: sum-free, "useful negative-control role" | **DRIFT for public use; ANCHOR ALREADY ACKNOWLEDGES NEG-CONTROL** — internally fine, externally must be relabeled | same as #755 |
| **905** | TIER_1_PROTOTYPE_NOW | **solved** (local DB) | SI 0.992 rank 6; phys_sim 0.8754, 8 phys edges | anchor: additive bases / extremal | **DRIFT for public use** | same as #755 |
| **505** | TIER_1_PROTOTYPE_NOW | **disproven** (local DB) | SI 0.413 rank 92; phys_sim 0.782, 12 phys edges | anchor: covering deficit | **DRIFT** — a disproven covering result is a different kind of test; the screen treats it identically to an open one. Re-label or move to a separate "disproven-as-control" lane | confirm disproof citation in D1; either drop or re-anchor |
| **161** | TIER_1_PROTOTYPE_NOW | open | SI 0.994 rank 183 (NB: list rank vs si_score normalization differ); phys_sim 0.8615, 13 phys edges | anchor: discrepancy jumps, phase-boundary | **CONSISTENT** | proceed to statement+lit audit |
| **20** | TIER_1_PROTOTYPE_NOW | open | SI 0.993 rank 62; phys_sim 0.7426, 18 phys edges | anchor: sunflower, set-system tensors | **CONSISTENT** | proceed to statement+lit audit |
| **30** | TIER_1_PROTOTYPE_NOW | open | SI 0.876 rank 178 in solvent_v2; phys_sim 0.797, **only 4 phys edges** | anchor: Sidon, transfer-matrix surface | **CONSISTENT but soft signal** — the anchor boost (0.18) is doing most of the lifting; raw SI is rank 178 of 1,190 and physics edges are sparse. Score 0.83 is anchor-driven, not data-driven | proceed; note that EHP114-style transfer-operator work on Sidon (Erdos30 PMF) is the externally validated leg, not the screen's signal |
| **86** | TIER_2_AUDIT_THEN_PROTOTYPE | open | SI 0.865 rank 108; phys_sim 0.7062, 3 phys edges | anchor: C4-free hypercube, spectral graph geometry | **CONSISTENT** | proceed to statement+lit audit |
| **1120** | TIER_3_WATCHLIST | open | SI 0.486 rank 488; **phys_sim 0.9552 (highest in priority set), 23 phys edges (highest count)** | anchor: polynomial radial path, lemniscate-adjacent | **PARTIALLY INCONSISTENT signal** — physics-edge profile is the strongest of the entire priority set, but solvent-index rank (488) drags the composite down to TIER_3. The report's claim that #1120 is "negative-control / neighbor" of #114 is not actually borne out by the screen score; the physics-edge data alone would put it in Tier 1 | proceed as the user's stated negative-control / neighbor for Lane A; flag the SI/phys-sim disagreement |

Spot-check on the user's secondary list (#138, #765, #165, #975, #707, #234, #198, #231, #944, #188, #16, #217):

| problem_id | local_tier | db_status (local) | flag |
|---:|---|---|---|
| 138 | TIER_1_PROTOTYPE_NOW | open | OK |
| 765 | TIER_1_PROTOTYPE_NOW | **solved** | DRIFT |
| 165 | TIER_1_PROTOTYPE_NOW | open | OK |
| 975 | TIER_1_PROTOTYPE_NOW | open | OK |
| 707 | TIER_1_PROTOTYPE_NOW | **disproven** | DRIFT |
| 234 | TIER_1_PROTOTYPE_NOW | open | OK |
| 198 | TIER_1_PROTOTYPE_NOW | **disproven** | DRIFT |
| 231 | TIER_1_PROTOTYPE_NOW | **disproven** | DRIFT |
| 944 | TIER_1_PROTOTYPE_NOW | open | OK |
| 188 | (not in priority-20 in screen output) | open | OK |
| 16 | TIER_1_PROTOTYPE_NOW | open | OK |
| 217 | TIER_1_PROTOTYPE_NOW | open | OK |

Across the screen's top 20, **7 rows have a non-`open` `math_status` in the local DB** (#166, #198, #231, #505, #707, #755, #765, #905 — that's 8 across the top 20, of which 6 are in the top 12). The pool itself contains 445 solved + 117 disproven + 18 partially solved out of 5,496. Solved/disproven are not rare — the screen is just choosing not to filter.

## Three-lane verdict

**Lane A — EHP114 / 1120 polynomial-potential geometry.** The strongest of the three on local-DB evidence: #114 has 22 physics edges with top-sim 0.9463 and #1120 has 23 with top-sim 0.9552 — these are the densest physics-edge profiles in the entire priority set. Anchor narrative (lemniscate length deficit / polynomial radial path) is internally coherent. The downstream Lean state of #114 is **NO COMPILED PROOF** per the canonical Math/CLAUDE.md ledger; the screen's TIER_0_CALIBRATION label refers to small-n compute, not a verified proof. This is fine for a methods-applicability screen but the language should not migrate into any artifact that mixes "calibration" with "compiled".

**Lane B — Sidon / B_h[g] / sum-free constrained-state geometry.** Mixed quality. #30 (Sidon) is the externally legible leg and has Erdos30 PMF / transfer-matrix work behind it — that's the published-strong story. But #755 (B_h[g]) is solved per local DB and #166 (sum-free) is solved per local DB. As an internal routing question this lane is still useful (testing the certification machinery on a known-resolved constrained-state problem is a legitimate calibration), but the *external-facing* version of this lane has to be #30 only, with #755/#166 explicitly framed as control / replication targets, not opportunities.

**Lane C — hypercube / discrepancy / graph spectral geometry.** Cleanest of the three for opportunity framing: #161 (discrepancy jumps, open) and #86 (C4-free hypercube, open) are both `open` in local DB and have plausible operator/spectral framings. #905 (additive bases) is solved and shouldn't carry the lane on its own. Of the three lanes, this is the one where the methods-applicability claim has the least baggage — but it's also the lane where the screen has the *least* feature signal, since `spectral_operators` is empty. The "spectral graph geometry" framing is anchor prose, not data-checked.

Composite read across the three lanes: Lane A (polynomial-potential) has the strongest local data; Lane B is fine for internal calibration but needs careful external phrasing; Lane C is the cleanest "open-problem opportunity" lane but is the most data-thin in the screen as it currently runs. **All three lanes are still worth pursuing as internal prototype experiments** — the verdict is not that the strategy is wrong, just that the screen as built does not by itself license a public ranking.

## Red flags

1. **Generator does not read `math_status`.** `scripts/score_sktc_applicability.py` selects from `problems` joined with `solvent_index_v2` and `physics_edges` but never inspects `math_status`. Solved/disproven problems land in TIER_1 alongside open ones. **Six of the top twelve rows are not open problems** in the local DB. For internal routing this is recoverable; for any external surface it is a labeling defect.

2. **Empty `spectral_operators` table.** The 7-feature applicability signature implies operator-class verification. There is none — the table has 0 rows and the scorer never queries it. The signature is anchor prose. This is a "looks like more than it is" risk if the screen is shown to a third party as evidence the strategy maps to operator-level structure across 1,200 problems.

3. **`erdos_problems.db` at `/Math/erdos_problems.db` is 0 bytes.** Any pipeline still reading from that path is silently reading nothing. The `Math/CLAUDE.md` Attack Queue Pipeline section instructs running `build_morphism_index.py` against `erdos_problems.db` — that pipeline cannot be live with a 0-byte source. (Out of scope for this review, but it directly affects whether the SI scores in `erdosatlas.db` are themselves fresh.)

4. **D1 is the canonical store and was not reachable.** Per `Math/CLAUDE.md` (Lean Proof Status Pipeline) and the global rules, no Lean status / scoring number quotes are valid from local cache — they must come from D1. This review could not perform the canonical D1 reconciliation. A `CONSISTENT` verdict is not licensable until the D1 round-trip is run.

5. **Rank vs. score normalization mismatch in the screen output.** Some rows have a high `si_score` value (e.g., #161 at 0.994) but rank far down in `solvent_index_v2.si_rank` (rank 183). The screen uses the value, not the rank — that's a defensible choice, but it means readers should not interpret "high SI score" in the screen as "high SI rank." Document the scoring intent in the screen header.

6. **Anchor boost is doing heavy lifting on #30.** Raw SI is 0.876 (rank 178), physics edges are 4. The Tier-1 placement comes from the +0.18 anchor boost. This is fine if the anchor is honest prior knowledge; it is a problem if a downstream consumer reads the score as "data discovered Sidon is high-priority". For internal use, label the anchor-boosted scores as such.

7. **Token embargo (`Feynman` / `feynman`) is not at risk in this report**, and was not used. The `Math/CLAUDE.md` "Feynman-Plain" mention is in a project-local style header — flag for separate migration to "plain-language" per global rule, but that is out of scope here.

## Recommended next steps

**Narrow technical (one).** Modify `scripts/score_sktc_applicability.py` to (a) read `problems.math_status` and emit a `db_status` column on every output row, (b) add a `tier_modifier` that demotes `solved` and `disproven` rows out of TIER_1 unless they are explicitly tagged as calibration / negative-control, and (c) emit a "data vs. anchor" decomposition so readers can see how much of each score is anchor-driven versus DB-driven. Re-run, re-hash, archive the current results under a `_v1_anchor_only/` suffix and version the new run as `EXP-MATH-ERDOS-SKTC-APPLICABILITY-SCREEN-20260502-02`. Do not overwrite the existing `_RESULTS.json`.

**Registry cleanup (one).** Drain D1 (or schedule the drain when token is available) and reconcile the screen's `db_path` reality. Either delete the 0-byte `Math/erdos_problems.db` (with approval, per Rule 3 — never delete without confirming this is not historical state someone else is pointing to) or document loudly that it is a stub and the live local store is `erdosatlas-workbench/erdosatlas.db`. Add a session-start health check that fails if a DB referenced by any pipeline is 0 bytes.

**Public-claim boundary (one sentence, verbatim, mandatory in any external surface that derives from this screen):** *"We have an internal routing screen suggesting that the tensor/operator/cone-certificate strategy developed around EHP114 may be worth testing across three families of Erdős problems: polynomial-potential geometry, constrained-state additive geometry, and graph-spectral/discrepancy geometry. This is not a proof, not a solved-problem claim, and not a public morphism verdict."*

## Status vocabulary used in this review

- **Hypothesis** (this is not a hypothesis status update — the screen does not propose any hypothesis): n/a.
- **Scorecard claim:** the screen does not generate a scorecard claim; it generates a routing list. If a scorecard row is to be derived from it, the only honest current status is `UNTESTED`.
- **Evidence binding:** the screen does not produce evidence binding. If used to anchor a claim, the binding status would be `NEEDS_ASSEMBLY` (no prototype artifact yet exists for any of the three lanes besides the existing #30 PMF work and the small-n #114 compute).
- **Lean:** never quoted CLEAN as COMPILED in this review. #114 is **not COMPILED**; only small-n compute exists. #30 has multiple `pass` artifacts in the local scoring DB (`Erdos30_AdditiveEnergy.lean`, `Erdos30_BFR.lean`, `Erdos30_Singer.lean`, etc., last_verified 2026-04-12) but per the canonical pipeline a fresh D1 query is required before any of those can be quoted publicly today.

## Read-only confirmation

No writes were performed against any DB during this review — not D1 (unreachable), not `erdosatlas.db` (read-only `SELECT`), not `scoring/erdos_scoring.db` (read-only `SELECT`). No remote state was altered. No git commits, pushes, PRs, or social posts were attempted. The only file written by this review is the markdown report at this path.

## Artifact hashes (verified)

```
RESULTS.json   10f3463de2480ee792866061c1ee60940486085bb792df086d6bde08664040ba   ✅ matches manifest
REPORT.md      843ea08e538f854c1079695299c48e619b917179d5a3e16c80db2f823bf05e65   ✅ matches manifest
```
