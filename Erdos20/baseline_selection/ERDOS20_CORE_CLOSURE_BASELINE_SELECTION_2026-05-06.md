# Erdos #20 Core-Closure Baseline Selection

**Date:** 2026-05-06
**Worker:** #20-D
**Lane:** sunflower core-closure / Leg-4 precursor
**Scope:** local baseline-selection packet only. Writes are confined to `erdos-experiments/Erdos20/baseline_selection/`.
**Claim ceiling:** A-axis remains **A0**. This is a **shadow signature, not universal law**.

## Bottom Line

Use the floor-normalized fixed-core local closure cost for the next run:

```text
I_core_local(s,m) / I_floor(s,m)
```

The numerator should be the already measurable fixed-core local information channel:

```text
I_core_local(s,m)
```

Start with the saved best precursor:

```text
w = 3, fixed core size s = 2, m near m* and across m* +/- delta
```

The Abbott-Hansen-Sauer-normalized quantity should be reported as a secondary construction/literature control, not used as the primary denominator for the next run.

Meaning: the next run should stop asking only "is there local closure information?" and start asking "how far above the predeclared floor is this local closure channel?" That is the question Leg 4 actually needs. AHS remains essential because the April literature gate says the small-n measurements sit below the known asymptotic construction story; but AHS is a comparator, not the physical floor.

## Inputs Inspected

- `erdos-experiments/Erdos20/experimental_deepening/ERDOS20_CORE_CLOSURE_DEEPENING_PACKET_2026-05-05.md`
- `erdos-experiments/Erdos20/experimental_deepening/EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_RESULTS.json`
- `erdos-experiments/Erdos20/formal_core_closure/ERDOS20_FORMAL_CORE_CLOSURE_TARGETS_2026-05-05.md`
- `erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md`
- `erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md`
- `erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md`

## Diagnostic Artifact

I added and ran a saved-artifact diagnostic:

```text
python3 erdos-experiments/Erdos20/baseline_selection/run_baseline_selection_diagnostic.py
```

It produced:

- `EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.json`
- `EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_REPORT.md`
- `EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.sha256`

The diagnostic did not enumerate new families. It read the saved May 5 deepening result, the May 5 interpretation packets, and the April reports/literature gate, then ranked the three candidate quantities by immediate computability and Leg-4 suitability.

## Candidate Ranking

| Rank | Candidate | Quantity | Score | Interpretation |
|---:|---|---|---:|---|
| 1 | floor-normalized local closure | `I_core_local(s,m) / I_floor(s,m)` | 9 | Best Leg-4 target if `I_floor` is predeclared before the run. It preserves the measured per-core channel and adds the missing floor denominator. |
| 2 | observed local closure only | `I_core_local(s,m)` | 7 | Best immediate numerator, but it has no baseline. More of this alone keeps the lane in precursor mode. |
| 3 | AHS-normalized closure cost | `I_core_local(s,m) / I_AHS(s,m)` or closure excess over an AHS-style baseline | 5 | Necessary literature/construction control, but not ready as the primary quantity because no AHS construction artifact exists in this lane yet. |

## Why This Recommendation

The May deepening packet already found a real measurable channel: per-core local closure cost near jamming. The strongest saved precursor is the `w=3`, fixed-core-size `s=2` series. It is stable by the simple slope screen, with log-log slope `-0.076613` and late-window CV `0.182257`, but the packet correctly refuses to call this Leg-4 evidence because the signal is not floor-normalized.

Observed `I_core_local(s,m)` is therefore the right numerator but the wrong complete experiment. It tells us what the core has to remember, but not whether that cost is near a physical/information floor or merely another small-n geometric shadow.

The floor-normalized ratio is the clean next run because it keeps the numerator tied to the real per-core channel and forces the denominator to be declared before computation. That is the smallest move that can turn "closure signal present" into a real Leg-4 test candidate.

AHS normalization is not discarded. The April literature gate makes it unavoidable: Abbott-Hansen-Sauer dominates the small-n lower-bound story, so any serious #20 lane needs an AHS-style comparator. But using AHS as the primary denominator would change the question from "does the core channel approach a floor?" to "how far are we from a construction baseline?" That is useful, but it is not the Maxwell/Mendoza-floor Leg-4 question.

## Next Run Definition

Pre-register the next run as:

```text
primary_quantity = I_core_local(s,m) / I_floor(s,m)
numerator = fixed-core local closure information, I_core_local(s,m)
initial target = w=3, s=2, m in {m* - delta, m*, m* + delta}
baseline = predeclared Mendoza-floor denominator, I_floor(s,m)
secondary_control = AHS-style construction-normalized closure cost, reported separately
```

Minimum acceptance for that next run:

- `I_floor(s,m)` is defined before enumeration.
- The run records at least one defect/hysteresis window around `m*`, not only the jamming point.
- The report separates the floor ratio from the AHS construction control.
- Any verdict remains internal unless geometry-drift and construction-baseline controls are explicitly passed.

## Claim Ceiling

Safe statement:

> The #20 lane has selected `I_core_local(s,m) / I_floor(s,m)` as the next-run baseline. The current numerator is grounded in saved per-core closure artifacts, but no denominator has yet been executed. Status remains A0: shadow signature, not universal law.

Unsafe statements:

- The sunflower conjecture has been advanced.
- The lane has a Leg-4 pass.
- A Maxwell law governs sunflower closure.
- The current small-n values improve known lower bounds.
- The AHS construction has been reproduced or normalized here.

## Files Changed

- `erdos-experiments/Erdos20/baseline_selection/run_baseline_selection_diagnostic.py`
- `erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.json`
- `erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_REPORT.md`
- `erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.sha256`
- `erdos-experiments/Erdos20/baseline_selection/ERDOS20_CORE_CLOSURE_BASELINE_SELECTION_2026-05-06.md`

## Commands Run

```text
rg -n "Erdos20|sunflower|baseline|Abbott|Hansen|Sauer|core closure|#20" /Users/kenbengoetxea/.codex/memories/MEMORY.md
sed -n '1,220p' /Users/kenbengoetxea/.agents/skills/math-problems-manager/SKILL.md
sed -n '1,260p' erdos-experiments/Erdos20/experimental_deepening/ERDOS20_CORE_CLOSURE_DEEPENING_PACKET_2026-05-05.md
jq . erdos-experiments/Erdos20/experimental_deepening/EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_RESULTS.json
sed -n '1,260p' erdos-experiments/Erdos20/formal_core_closure/ERDOS20_FORMAL_CORE_CLOSURE_TARGETS_2026-05-05.md
sed -n '1,260p' erdos-experiments/Erdos20/Q1_LITERATURE_GATE_2026-04-17.md
ls -la erdos-experiments/Erdos20
find erdos-experiments/Erdos20 -maxdepth 2 -type f \( -name '*2026-04*' -o -name '*SUNFLOWER-001*' -o -name '*SUNFLOWER-002*' -o -name '*ABBOTT*' -o -name '*baseline*' \) | sort
sed -n '1,260p' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-001_REPORT_2026-04-17.md
sed -n '1,260p' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-SUNFLOWER-002_REPORT.md
ls -la erdos-experiments/Erdos20/baseline_selection 2>/dev/null || true
jq 'keys' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json
jq 'keys' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-20260505-01_RESULTS.json
jq 'keys' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-W3N8-20260505-02_RESULTS.json
jq '.experiments? // .observations? // .results? // .' erdos-experiments/Erdos20/EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-20260505-01_RESULTS.json | head -n 220
sed -n '1,220p' erdos-experiments/Erdos20/experimental_deepening/run_core_closure_deepening_diagnostic.py
mkdir -p erdos-experiments/Erdos20/baseline_selection
python3 erdos-experiments/Erdos20/baseline_selection/run_baseline_selection_diagnostic.py
ls -la erdos-experiments/Erdos20/baseline_selection
jq '.recommendation, .candidate_ranking, .recommended_next_run_quantity' erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.json
cat erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_REPORT.md
shasum -a 256 -c erdos-experiments/Erdos20/baseline_selection/EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.sha256
cd erdos-experiments/Erdos20/baseline_selection && shasum -a 256 -c EXP-MATH-ERDOS20-CORE-CLOSURE-BASELINE-SELECTION-20260506-01_RESULTS.sha256
sed -n '1,260p' erdos-experiments/Erdos20/baseline_selection/ERDOS20_CORE_CLOSURE_BASELINE_SELECTION_2026-05-06.md
find erdos-experiments/Erdos20/baseline_selection -maxdepth 1 -type f -print | sort
```

The first checksum command above was run from `Math/` and failed because the checksum file intentionally contains the result filename relative to `baseline_selection/`. The corrected command from the output directory passed.
