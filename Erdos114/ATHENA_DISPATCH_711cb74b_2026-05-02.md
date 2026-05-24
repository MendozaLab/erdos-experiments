# Aristotle/Athena Dispatch — Erdős #114 Radial Spike

**Dispatched:** 2026-05-02T22:21:33Z
**Project ID:** `711cb74b-e169-4bfb-b2df-5766ddb3f5b6`
**Status at dispatch:** QUEUED (1 sec into queue)

## Submission environment

- **CLI:** `aristotle` v1.0.0 (`aristotlelib` Python package, harmonic.fun) — NOT the OHDSI `athena` CLI which is a name collision
- **Submission dir:** `/tmp/aristotle_spike_ehp114_2026-05-02` (24 KB)
- **Lean toolchain in submission:** `leanprover/lean4:v4.28.0` (Aristotle's preferred pin; scaffold dev tree is on v4.27)
- **Mathlib pin in submission:** `v4.28.0`
- **Auth:** `ARISTOTLE_API_KEY` env var (Ken's Gmail-bound key in `~/.zshrc`)

## Source-of-truth scaffold

`/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean`
- 290 lines, 4 axioms (citation-backed), 8 instances of `sorry` text (2 are real proof-body sorries; 6 appear in doc-comments per `grep -c "sorry"`)
- Pre-dispatch SHA-256: `d99575d9eca2a1c109a82be581b5c6b5cf5bd00ac259cd43da6f4a91495a8b3d`
- Local build: `lake build EhpRadialPuiseux` PASS in 39s, 7886 jobs, 2 expected sorry-warnings (lines 249 + 288)

## Targets

**Target 1: `theorem ehp_radial_L_form_n14`**
- Goal: $L(z^{14} - a) = 2\pi \cdot {}_2F_1(13/28, 13/28; 1; a^2)$ for $0 \leq a < 1$
- Athena's job: specialize `hypergeometric_lemniscate_radial 14 ...` and conclude via `norm_num` arithmetic on $(14-1)/(2 \cdot 14) = 13/28$

**Target 2: `theorem ehp_radial_puiseux_n14`**
- Goal: $\exists C > 0$ such that $L(z^{14} - 1) - L(z^{14} - (1 - \varepsilon)) \geq C \cdot \varepsilon^{1/14}$ for $\varepsilon \in (0, 1/4]$
- Athena's job: derive $C$ explicitly from `gauss_2F1_connection_at_one` + the Puiseux expansion at $a = 1$

## Constraints (in the submission prompt)

- No new axioms beyond the four declared
- No `sorry`/`admit` in output
- No tactic-smuggling (`nlinarith` on transcendentals, `norm_num` on analytic identities)
- Permitted tactics: `rfl`, `norm_num` (rational only), `ring`, `ring_nf`, `simp`, `field_simp`, `omega`, `linarith` (linear only), direct axiom invocation

## Empirical anchor

Radial-Puiseux fit confirms slope = $1/n$ at machine precision n=3..15 (`Math/erdos-experiments/Erdos114/RADIAL_PUISEUX_FIT_2026-05-02.md`). 95% CI at n=18 is three orders of magnitude inside the tolerance threshold. The analytic claim is iron; Athena's job is the constant extraction in Lean.

## Polling protocol

```bash
# Check status
aristotle list | head -3

# Get the result once status reaches COMPLETE or COMPLETE_WITH_ERRORS
aristotle result 711cb74b-e169-4bfb-b2df-5766ddb3f5b6 \
  --destination /tmp/aristotle_result_ehp114_2026-05-02

# Or block-and-wait (only if you want to babysit):
aristotle result 711cb74b-e169-4bfb-b2df-5766ddb3f5b6 \
  --wait \
  --destination /tmp/aristotle_result_ehp114_2026-05-02
```

## Success / failure metrics

**SUCCESS:** `grep -c "^[[:space:]]*sorry$\|by sorry" <output>/EhpRadialPuiseux.lean` returns 0 AND `grep -c "^axiom " <output>/EhpRadialPuiseux.lean` returns 4 (unchanged) AND `lake build EhpRadialPuiseux` exits 0 in our local v4.27 environment.

**PARTIAL_TARGET_1:** Step 1 closes (`ehp_radial_L_form_n14` body sorry-free) but Step 2 still has sorry. Worth a follow-up dispatch with explicit Puiseux-coefficient hint.

**PARTIAL_TARGET_2:** Step 2 closes but Step 1 doesn't. Less likely (Step 2 depends on Step 1). If this happens it suggests Athena introspected a different proof structure.

**FAILURE:** Both sorries remain, OR new axioms introduced, OR axiom count drops (Athena replaced an axiom with a "proof"), OR ≥5 days wall-clock with no progress.

## Telemetry to capture post-completion

| Field | Where | Meaning |
|---|---|---|
| Total wall-clock | from `aristotle list` CREATED → COMPLETE diff | Real spike duration |
| Status at completion | `aristotle list` STATUS column | COMPLETE / COMPLETE_WITH_ERRORS / FAILED / CANCELED |
| sorry-count delta | grep on output vs scaffold | Primary success metric |
| axiom-count delta | grep `^axiom ` on output vs scaffold | Integrity check |
| Lemmas added by Athena | diff scaffold vs output, count new theorem/lemma decls | Proof-architecture telemetry |
| Tactic palette used | grep tactic invocations in proof body | Methodology telemetry |
| Any dropped/replaced axiom | structural diff | If non-zero, halt — possible tactic-smuggling |

## Known risks

1. **Lean v4.27 vs v4.28 mismatch** — submission pinned to v4.28 (Aristotle's preferred), our dev tree is v4.27. Lemma renames or drift between v4.27 and v4.28 might surface. If output uses v4.28-specific names that don't exist in v4.27, the local-rebuild check will need a separate v4.28 build environment, or hand-translation.

2. **Mathlib hypergeometric infrastructure** — Athena may reach for ${}_2F_1$ lemmas that exist in v4.28 Mathlib but not in our v4.27 pin. Our scaffold uses a locally-defined `gauss_2F1` wrapper to avoid Mathlib drift; Athena may try to bypass it.

3. **Cooley filter** — this artifact stays internal until ErdosAtlas provisional is filed. The submission to Aristotle's hosted service IS an external disclosure of the proof scaffold — Aristotle's TOS / data retention should be checked before any v6 preprint cites the resulting Lean theorem.

## Post-completion analysis checklist

When the project finishes:

1. `aristotle list | grep 711cb74b` — confirm status
2. `aristotle result 711cb74b-... --destination <out>` — pull files
3. `diff <scaffold> <out>/EhpRadialPuiseux.lean` — see what changed
4. `cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30 && lake build EhpRadialPuiseux` — local v4.27 verify
5. Compute telemetry table → write `EXP-MATH-EHP114-ATHENA-SPIKE-V1-20260502_RESULTS.json` with all metrics
6. Update D1 `scoring_lean_artifacts` (new row `EhpRadialPuiseux.lean`, theorem_count, sorry_count, lake_build_status, last_verified_date)
7. Re-estimate the n=14..18 Lean program's true cost using this single datapoint as the anchor
