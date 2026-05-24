# Erdos #20 Core-Closure Deepening Packet

**Date:** 2026-05-05
**Lane:** Sunflower core-closure / Maxwell Leg-4 precursor
**Scope:** Internal experimental deepening only. No public claim. No D1, scorecard, packet, git, or Downloads update.
**Claim ceiling:** A-axis remains **A0**. This is a **shadow signature, not universal law**.

## Bottom Line

The #20 lane now has a real measurable object: per-core local closure cost near the jamming point. That is stronger than the earlier aggregate lattice-gas story because it asks what a fixed shared core has to remember as petals accumulate.

It is still not a sunflower theorem, not lower-bound progress, and not a Leg-4 pass. The Maxwell language remains analogy: a displacement-current-style closure term is a useful way to look for missing bookkeeping, not a physical law acting on set systems.

## Existing Evidence

The existing evidence is experimental and diagnostic:

- `SUNFLOWER_CORE_CLOSURE_LEG4_PACKET_2026-05-05.md` defines the right observable and explicitly blocks public theorem or lower-bound framing.
- `EXP-MATH-ERDOS20-HESSIAN-CLOSURE-20260505-01_RESULTS.json` shows aggregate closure pressure and curvature near jamming, but classifies the signal as `HESSIAN_INCONCLUSIVE` because aggregate data is not per-core Leg-4 evidence.
- `EXP-MATH-ERDOS20-PER-CORE-CLOSURE-20260505-01_RESULTS.json` instruments exact cheap regimes and classifies `PER_CORE_SIGNAL_PRESENT`; its own ceiling says this is only a precursor, with no floor-normalized ratio and no theorem-status change.
- `EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-20260505-01_RESULTS.json` moves the per-core runner into Rust for `w=3,n=7` and `w=4,n=7`.
- `EXP-MATH-ERDOS20-RUST-PER-CORE-CLOSURE-W3N8-20260505-02_RESULTS.json` extends `w=3` to `n=8`, giving the first useful four-point exact per-core jamming series for `w=3`.

What remains analogy:

- Maxwell does not prove anything about sunflowers.
- Mendoza-floor language has not yet been bound to a predeclared numerator.
- Abbott-Hansen-Sauer-style construction baselines have not been executed here.
- The exact small-n sunflower-free counts do not improve the known lower-bound story.
- No Lean artifact or theorem beyond encoding was produced.

## Deepening Diagnostic Run

I ran a deterministic secondary diagnostic from saved artifacts only:

```text
python3 erdos-experiments/Erdos20/experimental_deepening/run_core_closure_deepening_diagnostic.py
```

Artifacts produced:

- `EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_RESULTS.json`
- `EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_REPORT.md`
- `EXP-MATH-ERDOS20-CORE-CLOSURE-DEEPENING-20260505-01_RESULTS.sha256`

The run did not enumerate new families. It read the saved Hessian, Python per-core, and Rust per-core artifacts, then computed a drift screen for jamming-time local closure cost against ambient state count `N = binomial(n,w)`.

## Diagnostic Result

Classification: `DEEPENING_GEOMETRY_DRIFT_PRESENT`.

Meaning: the saved data contain a real per-core closure channel, but ordinary geometry still explains too much of the pattern to call this Leg-4 evidence.

Key series:

| Series | points | log-log slope vs N | late CV | screen |
|---|---:|---:|---:|---|
| `w=2` strongest per-core | 5 | 0.625853 | 0.123378 | geometry drift |
| `w=3` strongest per-core | 4 | 0.061889 | 0.106989 | stable by slope only |
| `w=4` strongest per-core | 2 | n/a | n/a | insufficient points |
| `w=3`, fixed core size `s=2` | 4 | -0.076613 | 0.182257 | stable by slope only |

The useful signal is `w=3`, especially fixed core size `s=2`, because it stays relatively stable from `n=5..8`. But that is only a slope screen. It is not a floor-normalized result, and it is not enough to say the Maxwell/core-closure analogy survives Leg 4.

## Next Leg-4 Measurable

The next measurable should be:

```text
I_core_local(s,m) near m* and across m* +/- delta
```

In plain terms: fix a core size `s`, look near the jamming point `m*`, and measure how many bits the core-local channel needs to distinguish locally safe petals from blocked petals. Then repeat across a small defect/hysteresis window around `m*`, not just at one family size.

Why this quantity:

- It keeps the core fixed instead of letting the strongest core size switch with `n`.
- It separates local core closure from aggregate family growth.
- It gives a direct hook for defect/hysteresis: before jamming, at jamming, and after the channel begins to collapse.

What must be added before any Leg-4 verdict:

- A predeclared Mendoza-floor numerator.
- An Abbott-Hansen-Sauer or construction-normalized baseline.
- At least one three-point geometry-exhausting sweep beyond the current `w=3` exact series, or an honest reason why `w=3` is the only relevant core-closure class.
- A windowed Rust/symmetry-reduced run that records `m* - delta`, `m*`, and `m* + delta`.

## Claim Ceiling

Safe statement:

> The #20 sunflower lane has a measurable per-core closure-pressure channel. The current best precursor is the `w=3`, fixed-core-size `s=2` jamming series, which is stable by a simple slope screen but not yet Leg-4 evidence. The status remains A0: shadow signature, not universal law.

Unsafe statements:

- The sunflower conjecture has been advanced.
- A Maxwell law governs sunflowers.
- The lane has a power morphism.
- The small exact counts improve known lower bounds.
- The current data prove a thermodynamic floor.

## Next Blocker

The blocker is not computation alone. The blocker is the numerator/baseline definition. Before another heavy Rust run, pre-register exactly what is being normalized:

1. observed bits only: `I_core_local(s,m)`;
2. floor-normalized bits: `I_core_local / I_floor` after defining `I_floor`;
3. construction-normalized closure cost against an Abbott-Hansen-Sauer-style baseline.

Until that is fixed, more enumeration only deepens the precursor; it does not convert analogy into Leg-4 evidence.
