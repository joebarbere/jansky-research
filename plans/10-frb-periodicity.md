# Plan 10 — FRB repeater activity-periodicity search 🚀

> Context: third backlog slice. Scope: small (reuses the FRB infrastructure). Status: validation.

## Context

Some repeating FRBs have periodic activity windows (FRB 20180916B, 16.35 d; CHIME/FRB 2020). This
slice searches the CHIME Catalog-1 repeaters for such periods with a phase-folding Rayleigh
periodogram — reusing the catalogue ingestion from the FRB burst-statistics slice. The known
16.35-day period is the validation target; the honest expectation is *recover-the-known + limits*,
not a new discovery (Catalog 1 is sparse).

## Deliverables

- `src/jansky_research/frbperiod.py` — `rayleigh_z2`, `period_search` (+ exposure-blind
  `false_alarm_prob`), `search_repeaters` (per-source), `run()` (writes metrics + per-repeater CSV +
  a periodogram figure), and `synthetic_periodic_arrivals` as the offline fixture. Tested to the
  85% floor.
- A real run over the 18 Catalog-1 repeaters → `survey/period_results.csv`.
- `survey/period-findings.md` — the honest assessment (recovery vs limits, exposure/aliasing caveats).

## Approach

Reuse `pipeline.build_catalog` for the catalogue (one row per event); group by `repeater_name`; run
the periodogram on sources with `>= min_bursts`. **Caveats to surface:** CHIME's transit exposure
(the FAP is an upper bound; the real confirmation is the match to the published period), daily
aliasing, and the sparse sampling (most repeaters unsearchable).

## Verification

- `make cov` ≥85% on synthetic fixtures (the offline run recovers an injected 16.35-day source).
- Real run recovers FRB 20180916B's $\sim$16.35-day period as the sanity check.
- **science-reviewer** gate on the method, the significance interpretation, and the citations.
