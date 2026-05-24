# ERDOS89_DISTANCE_PROJECTION_SHADOW_AUDIT_20260514-02

## Meaning

This is a finite projection-shadow scout for Erdos #89. It measures how
many point pairs collapse onto the same squared-distance value. The signal
is diagnostic: grids carry high multiplicity shadows, while deterministic
perturbations and random controls mostly restore pair-level distinguishability.

Claim ceiling: finite diagnostic only. This does not improve Guth-Katz,
close the sqrt-log gap, or change the open status of the conjecture.

## Results

| configuration | n | pairs | distinct d^2 | projection loss | energy | max mult | pinned min/median/max |
|---|---:|---:|---:|---:|---:|---:|---:|
| lattice_4x4 | 16 | 120 | 9 | 0.925 | 2072 | 24 | 5/8.0/9 |
| lattice_6x6 | 36 | 630 | 19 | 0.969841 | 28780 | 80 | 9/16.0/19 |
| lattice_8x8 | 64 | 2016 | 33 | 0.983631 | 182624 | 168 | 14/25.5/33 |
| perturbed_lattice_6x6 | 36 | 630 | 630 | 0.0 | 630 | 1 | 35/35.0/35 |
| random_integer_36_seed_8901 | 36 | 630 | 629 | 0.001587 | 632 | 2 | 35/35.0/35 |

## Same-n Comparison

For n=36, the lattice distinct-ratio is 0.030159,
the deterministic perturbation ratio is 1.0,
and the random integer control ratio is 0.998413.
That is the local projection-shadow signature: the lattice loses pair
information through repeated squared distances; generic controls mostly do not.

## Next Step

Use this artifact to guide a narrower #89/#94/#98/#217 triangulation pass.
Do not use it for solve, SOTA, or publication-readiness language.
