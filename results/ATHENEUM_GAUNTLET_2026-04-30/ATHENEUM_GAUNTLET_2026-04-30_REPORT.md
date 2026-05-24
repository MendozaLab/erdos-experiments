# ATHENEUM_GAUNTLET_2026-04-30 -- Atheneum Gauntlet

**Date:** 2026-04-30

## Meaning

This is a finite laboratory run, not a theorem claim. The Collider was asked to show whether its physics language produces measurable structure, transfers across related but distinct combinatorial systems, and demotes bad analogies.

The short answer: it passes the instrument test in the bounded PMF transfer setting, and it also records a real demotion. Sidon and B_2[2] expose field-sensitive extremal faces; sum-free stays rigid under the same observables; the hypercube Leg-4 candidate fails the stricter null test.

## Transfer Runs

| System | Rows | Parity / exactness | Ground states | Observable splits | Entropy ln range |
|---|---:|---|---:|---:|---|
| #30 sidon | 11 | h(n) + maximizer parity PASS | 9806 | 11 | 2.303 to 8.126 |
| #755 b2g2 | 11 | h-1 near-ground exact | 36234 | 11 | 1.792 to 9.702 |
| #166 sumfree | 11 | ceil(n/2) parity PASS | 28 | 0 | 0.693 to 1.099 |

## What Changed

- The same zero-temperature field probes picked different ground states in Sidon and B_2[2]. That means the object is often the whole extremal face, not one pretty witness.
- The same probes did not split the sum-free family over this range. That is a useful negative control: the instrument can say 'no response' under the same API.
- B_2[2] explicitly separated exact ground states from the h-1 near-ground layer with exact counts in this run, not capped lower bounds.
- The classifier backprop carried two ambiguous rows, but the stricter exp8_v2 run demoted the hypercube Leg-4 candidate to LEG4_FAIL.

## Pass / Fail / Ambiguous

| Case | Verdict | Why |
|---|---|---|
| Cross-problem transfer (#30 -> #755 -> #166) | PASS | All three engines ran fresh with bounded exact packets. Sidon and sum-free parity checks passed; B_2[2] produced exact h-1 near-ground counts. |
| Hypercube Leg-4 false positive | FAIL | exp8_v2 verdict is `LEG4_FAIL`; all three preregistered criteria failed. |
| Neutrino chiral row | AMBIGUOUS/DEMOTED | Classifier backprop kept it ambiguous, while exp9_v2 strict verdict is `CHIRAL_CLASS_REJECTED`; it passed only the stable-range criterion. |

## MDL Gate

MDL primitive tests returned code `0`. The raw pytest output is saved at `erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw/ATH-GAUNTLET-MDL-PREFLIGHT-2026-04-30.txt`.

## Raw Artifacts

- `erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw/ATH-GAUNTLET-SIDON-20-30-2026-04-30_RESULTS.json`
- `erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw/ATH-GAUNTLET-B2G2-D1-20-30-2026-04-30_RESULTS.json`
- `erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw/ATH-GAUNTLET-SUMFREE-D1-20-30-2026-04-30_RESULTS.json`
- `erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/ATHENEUM_GAUNTLET_2026-04-30_RESULTS.json`
- `erdos-experiments/results/CLASSIFIER_BACKPROP_2026-04-19` (bound historical demotion evidence, not rewritten)

## Rerun Commands

```bash
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30/rust-transfer-operator/target/release/erdos30-transfer-operator --n-min 20 --n-max 30 --frontier-k 3 --experiment-id ATH-GAUNTLET-SIDON-20-30-2026-04-30 --output-dir /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw
```
```bash
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30/rust-transfer-operator/target/release/b2g_transfer --n-min 20 --n-max 30 --g 2 --frontier-k 3 --near-ground-deficiency 1 --experiment-id ATH-GAUNTLET-B2G2-D1-20-30-2026-04-30 --output-dir /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw
```
```bash
/Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/Erdos30/rust-transfer-operator/target/release/sumfree_transfer --n-min 20 --n-max 30 --frontier-k 3 --prune-deficiency 1 --experiment-id ATH-GAUNTLET-SUMFREE-D1-20-30-2026-04-30 --output-dir /Users/kenbengoetxea/container-projects/apps/H2/Math/erdos-experiments/results/ATHENEUM_GAUNTLET_2026-04-30/raw
```
```bash
/Users/kenbengoetxea/miniconda3/bin/pytest /Users/kenbengoetxea/container-projects/apps/H2/Math/collider/tests/test_mdl_primitives.py -q
```
