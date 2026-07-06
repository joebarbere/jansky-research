# 82 — Meteor influx → sporadic-E coupling: station counts vs the Wallops digisonde (~200 km)

Status: 📋 planned (hardware-gated) — needs the meteor-scatter station (plan 81) running
continuously for ≥6 months; the GIRO/Madrigal data harness and all statistics can be built and
tested offline now

## Context

Metallic ions from meteor ablation are a source population for sporadic-E, but the station sweep
found the observational literature is professional-radar and seasonal-scale only — no published
daily lagged correlation between a continuous amateur meteor count and a nearby digisonde's foEs.
The Wallops Island digisonde (~200 km from Philadelphia, public via GIRO) makes this station
uniquely placed; Madrigal TEC adds an ionospheric control. Rigor is the novelty axis: the naive
correlation is confounded by season, tides, and geomagnetic activity, so the deliverable is the
*controlled* daily-scale analysis. Depends on plan 81's station running; the TLE-based
satellite-glint excision falls out as a by-product. This is a ≥6-month campaign slice.

## Deliverables

- `src/jansky_research/sporadice.py`: `fetch_giro_foes` (Wallops foEs/hmEs from GIRO, pragma),
  `fetch_madrigal_tec` (pragma), `daily_meteor_index` (from plan 81's classified event logs),
  `lagged_correlation` (daily cross-correlation with block-bootstrap significance),
  `confound_model` (season/tidal harmonics + Kp/Ap regressors; correlation on residuals),
  `satellite_glint_mask` (TLE pass windows → event excision, the by-product),
  `synthetic_campaign` (injected lag + confounders → recovery), `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic/offline fixtures — no sky data needed for tests);
  `papers/sporadice/`; `survey/sporadice-findings.md`; wiring.

## Approach

0. GATE 0: plan 81's station on air and logging; GIRO Wallops foEs and Madrigal TEC access
   paths live-verified (standing fable-ideas GATE-0: no data path was verified in the scan).
1. Tooling + synthetic recover-a-known, ahead of the campaign: `synthetic_campaign` injects a
   known meteor→foEs lag on top of seasonal/tidal/geomagnetic confounders; the pipeline must
   recover the lag only after `confound_model` residualization — and report a null when the
   injection is absent (the false-positive control is load-bearing here).
2. Data harness: GIRO/Madrigal fetchers + the daily meteor index from plan 81's logs, with the
   TLE glint mask applied; all buildable before six months of data exist.
3. Campaign leg (`# pragma: no cover`, ≥6 months spanning shower and non-shower seasons):
   daily lagged correlation on confound-residualized series; Perseids/Geminids as natural
   high-influx events; report effect size, lag, and significance from block bootstrap.
4. GATE-2 science review → paper (an honest bounded null is publishable — the daily-scale
   amateur-count constraint does not exist in the literature either way).

## Verification

Synthetic campaign recovers the injected lag post-residualization and yields a clean null
without injection; real-leg significance comes from block-bootstrap against the confound model;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Confounders (season/tides/geomagnetics) dominate** → they are modelled, not ignored; the
  falsifiable claim lives in the residuals, and the no-injection null test gates the method.
- **Station uptime over ≥6 months** → plan 81 inherits the `station/operations.md` watchdog/SLO
  architecture; gaps enter the bootstrap honestly.
- **Single-station counts are a biased influx proxy** → state it; corroborate high-influx days
  against GMN/RMOB (plan 81's harness) before leaning on any positive.
