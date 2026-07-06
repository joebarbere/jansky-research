# 62 — Merged GPPS+CRAFTS pulsar catalogue: the `lpt` compilation pattern at N≈1,000

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the GPPS bulk
table is live at zmtt.bao.ac.cn and that no merged FAST-discovery catalogue paper has appeared

## Context

FAST's two big pulsar surveys have discovered ~1,000 pulsars — 751 GPPS + >200 CRAFTS — but they
exist only as per-paper tables and the GPPS web table; no homogeneous, provenance-typed merge
exists with Galactic Z-heights, luminosity-function placement, or MSP fraction compared against
ATNF (fable-ideas F25). This is the `lpt` compilation pattern (per-row provenance flags,
discovery-paper spot-checks) at 60× the row count, reusing the merged `ppdot` (ATNF handling,
P–Ṗ plane) and `pulsarspec` (luminosity/spectral statistics) machinery. Data: the live GPPS bulk
table at `zmtt.bao.ac.cn/GPPS/GPPSnewPSR.html` (access pattern noted in earlier scans) + the GPPS
and CRAFTS paper tables + ATNF psrcat as the comparison population. Recover-a-known: the merged
catalogue must reproduce GPPS's own published faint-luminosity-tail claim before any new
population statement.

## Deliverables

- `src/jansky_research/gppsmerge.py`: `scrape_gpps_table` (bulk HTML table, `# pragma`),
  `parse_crafts_tables` (paper machine-readable tables, `# pragma`), `merge_catalogue`
  (provenance-typed rows: source survey, paper, flag column per the `lpt` discipline),
  `z_heights` (distance model + Galactic geometry), `luminosity_function` (placement vs ATNF,
  `pulsarspec` reuse), `msp_fraction` (P–Ṗ classing via `ppdot` reuse), `atnf_dedupe`
  (independent-rediscovery cross-match), `synthetic_merge` (fixture tables with planted
  conflicts → flag behaviour), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/gppsmerge/`; `survey/gppsmerge-findings.md`; wiring.

## Approach

0. GATE 0: verify the GPPS bulk table URL + column layout and locate the CRAFTS machine-readable
   tables; full-text pass on the GPPS series and CRAFTS discovery papers to confirm no merged
   catalogue (by the teams or others) has been published.
1. Tooling + synthetic recover-a-known: fixture tables with planted duplicates, unit conflicts,
   and missing fields; the merge must flag rather than silently resolve them.
2. Real leg: scrape + parse + merge; per-row spot-checks against discovery papers (the process
   that caught the `lpt` review's own data-file typo); Z-heights, luminosity function, MSP
   fraction vs ATNF; CPU, days.
3. GATE-2 science review: distance-model dependence (DM distances), survey selection-function
   differences between GPPS/CRAFTS/ATNF, transcription-error residual risk statement.
4. Paper: the merged catalogue (released as the artifact) + the population placement statistics.

## Verification

Merged catalogue reproduces GPPS's own faint-luminosity-tail claim before any new statistic is
reported; planted-conflict fixtures all flagged; spot-check sample clean; checks green; GATE-2
sign-off.

## Risks & mitigations

- **Transcription errors at N≈1,000** → keep the `lpt` flag discipline: every row carries
  provenance, every anomaly a flag; randomized spot-check sample against source papers is part
  of verification, and the residual error rate is stated in the paper.
- **Survey teams publish their own merge** → GATE-0 checks for in-press work; the compilation is
  fast (CPU, days) — ship bounded.
