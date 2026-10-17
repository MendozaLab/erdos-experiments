# Verification record: Erdős 123, `powers_2_3_5_snug`

Checked 2026-10-17 (UTC) on Ken's Mac.

## What is proved
`SnugBridge.erdos_123_powers_2_3_5_snug` is the exact statement of
`Erdos123.erdos_123.variants.powers_2_3_5_snug` in google-deepmind/formal-conjectures
(main, FormalConjectures/ErdosProblems/123.lean). The local `IsSnug` definition and the
theorem statement were compared character for character with upstream main: identical.

The mathematics is Principia Math's. The bridge only specialises their band theorem
`erdos123_dcomplete_real` to (a, b, c) = (2, 3, 5) with rho = min(1 + eps, 3/2).

## Results
| Check | Result |
|---|---|
| Principia erdos123 build (standalone clone) | Build completed successfully (8590 jobs), exit 0 |
| Principia axioms (`erdos123_dcomplete_real`, `erdos123_dcomplete'`) | propext, Classical.choice, Quot.sound |
| This package build (Principia as a Lake git dependency) | Build completed successfully (8590 jobs), exit 0, 07:48:53 UTC |
| `#print axioms SnugBridge.erdos_123_powers_2_3_5_snug` | propext, Classical.choice, Quot.sound |
| `sorry`, `native_decide`, extra `axiom` in SnugBridge.lean | none |
| Warnings with `-DwarningAsError=true` | none |

## Pins
- Lean toolchain: leanprover/lean4:v4.31.0
- Principia-Math-Solutions: 30d7e5327cdf5cda55137d70abec8f4d3b7ff2b3 (subDir erdos123, Apache-2.0)
- Mathlib: db127794c79fdeb86f6b0cf6ff2c804026fbaff1

## File hashes (sha256)
- SnugBridge.lean: 8fd317597f61bda3d667b693c194f6d5de7596784fbe0171217ebdddb55c0aa5
- lakefile.toml: 4494da83c81a19bf9951c090970579ebcd1c4f8b7d18a5cb37a76aaae8945505
- lean-toolchain: efac0b94923b2d8b6840cd35be9177ad0fc5ab2332f4f4311c98712cee92fdee

## Caveats
- formal-conjectures builds on Lean v4.33.1 with a newer Mathlib. This proof is checked on
  Principia's pins, so it is linked from the upstream file, not merged into it.
- erdosproblems.com does not yet list the snug variant as solved. The claim rests on
  Principia's paper and this kernel check.

## Reproduce
```
lake exe cache get
lake build
```
