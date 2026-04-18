# EXP-MATH-ERDOS30-SIDON-001 — Erdős #30 (Sidon Set Size)

> **Status:** Experimental probe — NOT a formal proof.
> Perplexity site-verification gate passed on 2026-04-07.

## Thermodynamic anomaly record
- **anomaly_id**: ERDOS-30
- **system_temperature_bound**: 0.0
- **volume_expansion_limit**: N^{1/2} + O_eps(N^eps)
- **gravitational_collapse_metric**: 0.98183
- **weak_decay_cycle_target**: Sidon L4 spectral floor
- **status**: UNRESOLVED

## Thermal noise exponent α (the conjecture's control variable)

| method | α | r² | conjecture target | Erdős–Turán ceiling |
|---|---:|---:|---:|---:|
| greedy_em | +0.913308 | 0.9917 | 0.00 | 0.25 |
| erdos_turan | +0.483399 | 0.9999 | 0.00 | 0.25 |
| singer_pds | -0.025114 | 0.9682 | 0.00 | 0.25 |

α → 0 is consistent with Erdős #30; α → 1/4 would be consistent with the Erdős–Turán tail being tight.

## Growth law 1 — greedy EM binding crystals

| N | |A| | gap (|A| − √N) |
|---:|---:|---:|
| 64 | 8 | +0.0000 |
| 128 | 12 | +0.6863 |
| 256 | 16 | +0.0000 |
| 512 | 20 | -2.6274 |
| 1024 | 28 | -4.0000 |
| 2048 | 36 | -9.2548 |

## Growth law 2 — Erdős–Turán quadratic-residue scaffolds

| p | ambient N | |A| | gap |
|---:|---:|---:|---:|
| 3 | 20 | 3 | -1.4721 |
| 5 | 54 | 5 | -2.3485 |
| 7 | 104 | 7 | -3.1980 |
| 11 | 252 | 11 | -4.8745 |
| 13 | 350 | 13 | -5.7083 |
| 17 | 594 | 17 | -7.3721 |
| 19 | 740 | 19 | -8.2029 |
| 23 | 1080 | 23 | -9.8634 |
| 29 | 1710 | 29 | -12.3521 |
| 31 | 1952 | 31 | -13.1814 |
| 37 | 2774 | 37 | -15.6688 |
| 41 | 3402 | 41 | -17.3267 |
| 43 | 3740 | 43 | -18.1555 |
| 47 | 4464 | 47 | -19.8132 |
| 53 | 5670 | 53 | -22.2994 |
| 59 | 7020 | 59 | -24.7854 |
| 61 | 7502 | 61 | -25.6141 |

## Growth law 3 — Singer perfect difference sets

| q | n = q²+q+1 | |A| = q+1 | elements |
|---:|---:|---:|---|
| 2 | 7 | 3 | [0, 1, 3] |
| 3 | 13 | 4 | [0, 1, 3, 9] |
| 4 | 21 | 5 | [0, 1, 4, 14, 16] |
| 5 | 31 | 6 | [0, 1, 3, 8, 12, 18] |
| 7 | 57 | 8 | [0, 1, 3, 13, 32, 36, 43, 52] |
| 8 | 73 | 9 | [0, 1, 3, 7, 15, 31, 36, 54, 63] |
| 9 | 91 | 10 | [0, 1, 3, 9, 27, 49, 56, 61, 77, 81] |

## Branch A — Singer rigidity

Overlap fraction between each greedy EM crystal and its nearest algebraic template. High = observed crystal is algebraically rigid; low = chaotic.

| greedy |A| | template | template |A| | overlap | shift |
|---:|---|---:|---:|---:|
| 8 | singer_pds | 8 | 0.3750 | 1 |
| 12 | erdos_turan | 7 | 0.4286 | -13 |
| 16 | erdos_turan | 11 | 0.3636 | -62 |
| 20 | erdos_turan | 17 | 0.2353 | -310 |
| 28 | erdos_turan | 23 | 0.1739 | 1 |
| 36 | erdos_turan | 31 | 0.1290 | 95 |

## Branch B — Fourier L⁴ spectral defect

- Sidon L⁴ sum:   1165312.0000
- Random L⁴ sum:  28242944.0000
- Ratio:          0.0413

Ratio below 1 suggests the Sidon set is spectrally flatter than a same-density random set, which is the signature a Branch B Cauchy-Schwarz sharpening would need.

## Required next step before any external publication

Per H² Formalization Integrity Protocol §5, run the Perplexity Pre-Submission Gate with the complete theorem statement (h(N) = √N + O_ε(N^ε)) before claiming novelty in the Branch A / Branch B attacks. The engine records a reference date but does not bypass the gate.

Wall clock: 5.39s
