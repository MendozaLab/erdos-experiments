# v12 changelog (vs v11)

| File / area | v11 | v12 |
|---|---|---|
| Status string | `LOCAL_INTERVAL_CERTIFIED` | `IEEE-1788 interval-certified` |
| `n = 13` row | shown in public table with annotation `route exception flagged` | quarantined from public table; documented in plain language in README, REPORT, and Zenodo strategy |
| erdosproblems.com post draft | bold headers + 14-row table inline + multiple questions | prose under 300 words, table moved to README, single question |
| Formal-conjectures packet | prose-only PR body | prose PR body + Lean 4 stub file with per-degree axioms |
| Zenodo decision | not addressed (assumed v6 mint) | explicit memo: edit v5 description; v6 only if n=13 is re-run or new degrees added |
| Process gates | not enumerated | enumerated in README and REPORT (prepub-redteam → Publisher → Crackpot-Scrub before deposit) |
| `n = 2` attribution | "MacLane / Eremenko–Hayman literature row" (combined) | distinguished: MacLane is historical, Eremenko–Hayman is accessible source pin; both cited per row |
| File count | 9 | 12 (added: `LEAN_STUB.lean`, `ZENODO_STRATEGY.md`, `CHANGELOG.md`; renamed: `RESULTS.json` schema updated) |

## Issues v12 still does not resolve

- The `inari` IEEE-1788 interval-arithmetic implementation itself is not formally verified. A reviewer who does not trust the toolchain has no recourse from this packet alone.
- The Tao threshold extraction is not attempted. The bridge between the finite packet and the all-degree statement remains open.
- No external endorser path is established for an arXiv submission; per memory, Moree is dead, Baez is gated.
