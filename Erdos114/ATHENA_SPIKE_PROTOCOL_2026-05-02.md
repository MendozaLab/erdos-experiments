# Athena Spike Protocol — EHP Radial Direction n=14

**Spike ID:** `EXP-MATH-EHP114-ATHENA-RADIAL-N14-V1-20260502`
**Date opened:** 2026-05-02
**Track:** Erdős #114 (Erdős–Herzog–Piranian) — radial-direction closure
**Status:** SCAFFOLD READY — awaiting Athena dispatch
**Cooley filter:** internal scaffolding; no public articulation
**Owner:** Ken Mendoza (sole inventor)

---

## 1. What this spike measures

Athena (Harmonic AI's Lean 4 prover, `athenalib` CLI) is asked to close two
arithmetic / substitution `sorry` stubs in a sorry-bearing scaffold for the
radial-direction closed form of the EHP lemniscate length at `n = 14`. The
deep classical analytic results required (Federer co-area, Fejér–Riesz,
Gauss connection at z=1, the radial closed form itself) are declared as
explicit `axiom`s with full citations and are explicitly **off-limits** to
Athena.

The spike answers a single operational question: **how much does Athena
cost (token spend, wall-clock) to close two `sorry` stubs that depend on
nothing more than `norm_num`/`ring`-class arithmetic plus a single
`rw` substitution from a hypothesis?** This datapoint anchors the
cost-curve for future Athena dispatches against the rest of the EHP
attack surface (n=15..18 closed forms, shape-mode bounds, Fejér–Riesz
feeder lemmas).

---

## 2. Scaffold file

**Path:** `Math/erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean`

**Build target:** `lean_lib EhpRadialPuiseux` (added to `Erdos30/lakefile.lean`
on 2026-05-02). Reasoning for placement: Erdos30 already has a fully resolved
Mathlib v4.27.0 build environment (~7.2 GB cached); creating a fresh package
for #114 would re-fetch Mathlib and burn a multi-GB on a 3–5 day spike with
~290 lines of Lean. Namespace `Erdos114.Radial` provides semantic separation;
no Erdos30 file imports the spike file.

**Build status (verified 2026-05-02):**
```
⚠ [7885/7886] Built EhpRadialPuiseux (29s)
warning: lean/EhpRadialPuiseux.lean:236:8: declaration uses 'sorry'
warning: lean/EhpRadialPuiseux.lean:275:8: declaration uses 'sorry'
Build completed successfully (7886 jobs).
```

**Inventory:**
| Count | Item |
|------:|------|
| 2     | `sorry` stubs (one per theorem) |
| 4     | `axiom` declarations (each with full bibliographic citation) |
| 2     | `theorem` declarations (`ehp_radial_L_form_n14`, `ehp_radial_puiseux_n14`) |
| 3     | `def` declarations (`lemniscate`, `ehp_lemniscate_length`, `gauss_2F1`) |

**Lean toolchain:** `leanprover/lean4:v4.27.0`
**Mathlib pin:** `v4.27.0` (rev `a3a10db0e9d66acbebf76c5e6a135066525ac900`)

---

## 3. Pre-spike issues to resolve

### 3.1 — Mathlib version mismatch (BLOCKING-IF-NOT-RESOLVED)

The aristotle-lean4 SKILL (`~/.claude/skills/aristotle-lean4/SKILL.md`)
states Athena is fixed to:
- Lean Toolchain: `leanprover/lean4:v4.24.0`
- Mathlib: `f897ebcf72cd16f89ab4577d0c826cd14afaafc7` (v4.24.0)

Our scaffold is on Mathlib v4.27.0. Two minor v4.24→v4.27 drifts the spike
file already adapts to:
1. `Complex.abs` deprecated (since 2025-08-26) → use generic `‖·‖` norm.
2. `Mathlib.Data.Complex.Norm` module deprecated as a whole.

**Resolution paths:**
- **(a) Run Athena with `--no-validate-lean-project --no-auto-add-imports`**
  flags. Athena treats the file as opaque text + lean syntax check, doesn't
  try to compile against its own pinned Mathlib. **Risk:** Athena's tactic
  attempts may use lemma names that exist in v4.24 but not v4.27, producing
  output that fails our local build. Mitigation: re-validate locally;
  if drift, patch the lemma names.
- **(b) Pre-flight via Athena on a reduced standalone file** containing
  only the two sorry-bearing theorems plus minimal context. Reduces the
  v4.24/v4.27 surface area Athena has to navigate.
- **(c) Backport the scaffold to v4.24** for the spike, then forward-port
  the closed proof. Costs an extra build dance but isolates the version
  question. **Recommended for the cost-curve datapoint** since it cleanly
  attributes Athena's effort to the mathematical content, not version-drift
  patching.

**Default for this spike:** Path (a). The two sorries are pure arithmetic
+ rewrite, so the lemma names Athena would use (`norm_num`, `ring`, `rfl`,
hypothesis rewrites) are stable across v4.24 → v4.27. Only escalate to
(b) or (c) if the first run produces v4.24-specific lemma names.

### 3.2 — Closed-form prefactor discrepancy (RESOLVED in scaffold)

The spike-task brief and `CROFTON_COAREA_DRAFT_2026-05-02.md` state
`L(z^n - a) = 2π · n · ₂F₁(...)`. The `RADIAL_PUISEUX_FIT_2026-05-02.md`
closed form (line 15) and the radial calibration JSON's `formula` field
state `L_n(a) = 2π · ₂F₁(...)` (no `n` factor).

**Resolution:** the `2π` form is correct. At `a = 0`, `p(z) = z^n` and
`Λ(p) = {|z| = 1}` is the unit circle (length `2π`); with `₂F₁(...; 1; 0) = 1`,
the `2π` form returns `2π` — correct. The `2π·n` form returns `2π·n` at
`a = 0`, which is wrong.

The Fryntov–Nazarov asymptotic `L(z^n - 1) = 2πn + o(n)` is recovered from
the `2π` form via the Gauss connection formula at `z = 1` — the `n`-prefactor
emerges from the Gamma quotient in the `a → 1⁻` limit, not from the
integrand.

**Action taken:** scaffold uses the `2π` form. The closed-form axiom
`hypergeometric_lemniscate_radial` is stated correctly. The target theorem
`ehp_radial_L_form_n14` matches the `2π` form. Diverges from the task
brief's stated theorem; documented inline in the scaffold's header comment
block. Per the H² Formalization Integrity Protocol §7, an axiom asserting
a mathematically false identity (the `2π·n` form) would constitute a
hidden integrity failure, so the resolution had to land on the
mathematically correct form.

### 3.3 — Perplexity Pre-Submission Gate (REQUIRED BEFORE PUBLIC PROMOTION)

Per `~/.claude/CLAUDE.md` Formalization Integrity Protocol §5, before any
Athena-completed version of this scaffold is promoted to a public artifact
(Zenodo, arXiv, formal-conjectures PR), the following Perplexity query
must be run and recorded in the file header:

> **Query:** Has the closed-form expression
> `L(z^n - a) = 2π · ₂F₁((n-1)/(2n), (n-1)/(2n); 1; a²)` for the lemniscate
> length of `p(z) = z^n - a` been published anywhere? Specifically check:
> Krishnapur–Lundberg–Ramachandran 2025, Fryntov–Nazarov 2025, the EHP
> survey literature, classical hypergeometric textbooks (Andrews–Askey–Roy,
> DLMF). Are there known formalizations of this identity in Lean 4 / Coq /
> Isabelle? What are the standard pitfalls?

If Perplexity reports the formula is published (e.g., in KLR2025), record
the citation in the axiom doc-comment and downgrade the axiom to a
`-- TODO: replace with Mathlib import` annotation. If novel, document as
"Mendoza Lab attribution; first published in v5 EHP preprint Theorem 1
(boundary case at `a = 1`); generic-`a` form derived from change of
variables `w = z^n` reducing to the elliptic period of `|w - a| = 1`."

The gate output goes into `EXP-MATH-EHP114-ATHENA-RADIAL-N14-V1-20260502_REPORT.md`
post-spike, alongside Athena's telemetry.

---

## 4. Dispatch — exact CLI invocation

### 4.1 — Pre-flight

```bash
# 1. Confirm the scaffold builds with current Mathlib (sorry warnings OK).
export PATH="$HOME/.elan/bin:$PATH"
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30
lake build EhpRadialPuiseux

# 2. Confirm Athena auth.
echo "$ARISTOTLE_API_KEY" | head -c 12   # should print "arstl_..."

# 3. Confirm athenalib is current.
pip install --upgrade athenalib
athena --version
```

### 4.2 — Primary dispatch (target: `ehp_radial_L_form_n14`)

The two sorries we want closed live in `ehp_radial_L_form_n14`. Athena's
default `prove-from-file` mode is the right tool.

```bash
export ARISTOTLE_API_KEY="arstl_..."
cd /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments

athena prove-from-file \
  Erdos30/lean/EhpRadialPuiseux.lean \
  --output-file Erdos114/ATHENA_OUTPUT_EhpRadialPuiseux_$(date +%Y%m%d-%H%M%S).lean \
  --no-validate-lean-project \
  --no-auto-add-imports \
  --polling-interval 30 \
  --max-polling-failures 240 \
  --context-folder Erdos114/
```

**Flag reasoning:**
- `--no-validate-lean-project` / `--no-auto-add-imports` — Athena's pinned
  Mathlib is v4.24, ours is v4.27. Skip the in-cloud project validation
  to avoid version-drift errors. Per skill SOP § Mode 1.
- `--polling-interval 30` — check every 30s.
- `--max-polling-failures 240` — Athena hard ceiling at this scaffold:
  240 polls × 30s = **2 hours wall-clock budget per call**. The 3–5 day
  budget is for *iteration* (multiple calls with refined hints), not a
  single mega-call.
- `--context-folder Erdos114/` — gives Athena the design docs (CROFTON,
  RADIAL_PUISEUX_FIT, TOEPLITZ_MOMENT_LIFT_DESIGN) as context so it
  understands what the axioms mean. Includes `*.md`, `*.json`, `*.tex`
  per skill SOP. (Excludes `.py` — Athena ignores those.)

### 4.3 — If Athena needs hints (Mode 1.b)

If the first run leaves `sorry` in `ehp_radial_L_form_n14`, edit the
scaffold to add `PROVIDED SOLUTION` doc-comment hints:

```lean
/--
Closed-form expression for the lemniscate length of z^14 - a (radial slice).

PROVIDED SOLUTION
Step 1: invoke `hypergeometric_lemniscate_radial 14 (by norm_num) a ha ha1`,
binding the result to `h_general`.
Step 2: rewrite the parameter `(14 - 1 : ℝ) / (2 * 14)` to `13 / 28` using
`norm_num` (this is pure arithmetic).
Step 3: After Step 2 and `rw [h_general]`, the goal is `rfl`-trivial.
Whole proof should be:
  `have h_general := hypergeometric_lemniscate_radial 14 (by norm_num) a ha ha1`
  `convert h_general using 3 <;> norm_num`
-/
theorem ehp_radial_L_form_n14 ... := by sorry
```

Then re-dispatch.

### 4.4 — Secondary dispatch (target: `ehp_radial_puiseux_n14`)

Only run after the primary target lands. The Puiseux corollary is harder
(needs `gauss_2F1_connection_at_one` + Gamma-function gymnastics for the
constant). Same CLI invocation; budget another 2-hour ceiling.

If Athena cannot close the Puiseux corollary within the spike's 5-day
budget, it stays `sorry` and the spike outcome is graded "partial" (see
Success/Partial/Failure table below).

---

## 5. Success / Partial / Failure metrics

### 5.1 — Success (full)

- `ehp_radial_L_form_n14` has zero `sorry`.
- `ehp_radial_puiseux_n14` has zero `sorry`.
- `lake build EhpRadialPuiseux` PASS with no warnings other than possibly
  unused-variable hints.
- Axiom count remains exactly **4** (Athena introduced no new axioms).
- All Athena-produced lemmas check at the local Mathlib v4.27.0.
- Wall-clock total ≤ 5 days.

### 5.2 — Partial credit

Two distinct partial paths:

| Partial outcome | What it tells us |
|---|---|
| **Step 2 closes, Step 3 doesn't** | Athena handles the arithmetic specialization `(14-1)/(2·14) = 13/28` but stalls on the substitution / `rfl`-conclusion. Likely indicates parser / definitional-unfolding issue with `gauss_2F1` notation. Worth one more dispatch with explicit hint. |
| **Step 3 closes, Step 2 doesn't** | Athena handles the substitution but balks at the rational arithmetic. Suggests `norm_num` extension or `ring` rather than direct numeric coercion. Worth one more dispatch with `norm_num` hint. |
| **Primary target closes, secondary doesn't** | Expected outcome at the easier end. Establishes the cost-curve baseline: Athena handles arithmetic + substitution but stalls on Gauss-connection-formula composition. |
| **Both sorries close but axiom count grew** | Failure of Formalization Integrity §7 — Athena introduced new axioms beyond the four declared. Roll back the run; do not promote. |

### 5.3 — Failure

- Any of:
  - Wall-clock > 5 days with no progress.
  - Axiom count exceeds 4 in the Athena output.
  - Athena uses a tactic to close any of the four declared axioms (this
    is hidden-sorry tactic-smuggling per H² Formalization Integrity
    Protocol §7).
  - Athena's output fails `lake build` against our local Mathlib v4.27.0
    and the failure is irreparable by trivial rename patches.
  - Athena returns a counterexample claim (would mean one of our axioms
    is mis-stated; recover the axiom statement and re-attempt).

In failure case: archive the output as `EXP-MATH-EHP114-ATHENA-RADIAL-N14-V1-20260502_FAILURE.lean`,
write a Feynman-skeptic post-mortem in the report, and re-evaluate whether
the spike should retry with the v4.24 backport path (option 3.1.c) or be
shelved pending a different prover.

---

## 6. Telemetry to capture

Each Athena call writes a server-side log accessible via the `athena status`
CLI. We capture:

| Field | Where | Why |
|---|---|---|
| `wall_clock_seconds` | from CLI start to output file written | Cost-curve datapoint |
| `polls_consumed` | CLI verbose output | Pacing of Athena's reasoning |
| `tokens_in / tokens_out` | post-call `athena status` | Direct token spend |
| `usd_cost_estimate` | computed from tokens × current rate | Compare to `aristotle-lean4` skill rate card if present; else use Anthropic API rates |
| `sorry_count_before / after` | `grep -c "sorry"` on input / output | Primary success metric |
| `axiom_count_before / after` | `grep -c "^axiom "` on input / output | Integrity check (must be unchanged) |
| `lemmas_added` | diff of `^theorem ` + `^lemma ` between input and output | Athena's added scaffolding |
| `lake_build_status` | `lake build EhpRadialPuiseux` exit code on output | Local re-validation |
| `tactic_inventory` | grep tactics used in Athena's output blocks | What Athena reached for |

Capture into `EXP-MATH-EHP114-ATHENA-RADIAL-N14-V1-20260502_RESULTS.json`
with SHA-256 sidecar at `_RESULTS.sha256`. JSON schema:

```json
{
  "experiment_id": "EXP-MATH-EHP114-ATHENA-RADIAL-N14-V1-20260502",
  "spike_target": "ehp_radial_L_form_n14",
  "scaffold_path": "Math/erdos-experiments/Erdos30/lean/EhpRadialPuiseux.lean",
  "scaffold_sha256": "<computed at dispatch>",
  "athena_calls": [
    {
      "call_index": 1,
      "wall_clock_seconds": 0,
      "polls_consumed": 0,
      "tokens_in": 0,
      "tokens_out": 0,
      "usd_cost_estimate": 0.0,
      "sorry_count_before": 2,
      "sorry_count_after": 0,
      "axiom_count_before": 4,
      "axiom_count_after": 4,
      "lemmas_added": [],
      "lake_build_status": "PASS",
      "tactic_inventory": []
    }
  ],
  "outcome": "SUCCESS|PARTIAL|FAILURE",
  "outcome_details": "<one paragraph>",
  "perplexity_gate_run": false,
  "perplexity_gate_result": null,
  "promoted_to_public": false
}
```

---

## 7. Post-spike analysis checklist

Before declaring the spike closed:

1. [ ] `lake build EhpRadialPuiseux` exits 0 with at most 2 sorry-warnings (success: 0).
2. [ ] Sorry count delta recorded.
3. [ ] Axiom count delta recorded (MUST be 0; non-zero is integrity failure).
4. [ ] No declared axiom closed by tactic in Athena's output (manual diff inspection).
5. [ ] Wall-clock and token cost recorded into RESULTS.json.
6. [ ] Cost-curve datapoint added to `Model-Mayhem/results/erdos-114/` (or
       wherever the cross-Athena cost curve lives — probably new file at
       `Model-Mayhem/results/athena-cost-curve.csv`).
7. [ ] If outcome SUCCESS: run Perplexity Pre-Submission Gate per §3.3 above
       before any public promotion.
8. [ ] If outcome SUCCESS: candidate for D1 update — `lake_build_status='PASS'`,
       `sorry_count=0`, `last_verified_date='<iso>'` against
       `scoring_lean_artifacts` row keyed on `lean_file_name='EhpRadialPuiseux.lean'`.
       Per `H2/Math/CLAUDE.md` Lean Proof Status Pipeline, status moves
       CLEAN → COMPILED.
9. [ ] If outcome FAILURE: post-mortem written; spike retry decision logged
       (continue, switch tactic, abandon).
10. [ ] Feynman-skeptic review of Athena's proof: did it teach us anything
        about the math? Or just close the symbol-pushing? Either is OK
        but we should know which.

---

## 8. One-line success/failure metric

**SUCCESS:** `grep -c "sorry" EhpRadialPuiseux.lean` returns `0` AND `lake build EhpRadialPuiseux` exits `0` AND `grep -c "^axiom " EhpRadialPuiseux.lean` still returns `4`, all within 5 days wall-clock.

---

## Provenance

- Scaffold built and verified `lake build` PASS on 2026-05-02.
- Toolchain `leanprover/lean4:v4.27.0`; Mathlib pin `v4.27.0` rev
  `a3a10db0e9d66acbebf76c5e6a135066525ac900`.
- Reference docs: `CROFTON_COAREA_DRAFT_2026-05-02.md`,
  `RADIAL_PUISEUX_FIT_2026-05-02.md`, `TOEPLITZ_MOMENT_LIFT_DESIGN_2026-05-02.md`.
- v5 preprint: `Math/erdosatlas-workbench/ehp_erdos114_preprint.tex`,
  Zenodo 10.5281/zenodo.19480329.
- Closed-form prefactor: `2π` (no `n`); discrepancy with task brief
  documented in §3.2.
- Cooley filter: this protocol document is internal; spike artifacts gain
  public articulation only after Perplexity Pre-Submission Gate (§3.3).
