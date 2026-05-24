# Ahmes Story: The Sums Are Sorted Into Rooms

After the representation count entered the main BFR shelf, Eratosthenes began sorting the sums into rooms.

The new scout does not ask how sharp the final bound is. It asks a more basic Hilbertian question: if the sum range is the hallway, what is an interval, and how much representation mass can sit inside one interval?

The answer is now Lean-readable. A clipped interval lives inside the BFR support range. Its occupancy is the sum of the ordered representation function across that interval. One such room cannot contain more mass than the whole hallway, and the whole hallway has mass `|A|^2`.

That is not the BFR discrepancy argument yet. It is the floorplan needed before that argument can be written without handwaving.

## Receipts

- Run artifact: `EXP-MATH-ERDOS30-BFR-INTERVAL-OCCUPANCY-SCOUT-20260508-01`
- Target: Erdős #30
- Persona: Eratosthenes of Cyrene
- Scribe: Ahmes
- Mathematical lens: David Hilbert
- Scout file: `scratch/Erdos30_BFR_IntervalOccupancy_SCOUT.lean`
- Scout SHA-256: `75a23b2e37833197d165692a2174988caaa2ee9d4f1e401fdd29dfd1ca845f6b`
- Module build: `lake build Erdos30_BFR`
- Scout check: `lake env lean scratch/Erdos30_BFR_IntervalOccupancy_SCOUT.lean`
- Escape-hatch matches: `0`
- Generated scout targets: `6`
- Decision: `INTERVAL_OCCUPANCY_SCOUT_READY_FOR_REVIEW`
- Review warning: this is interval bookkeeping, not a coefficient change or public result.
