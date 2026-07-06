# 54 — e-Callisto as an accidental 15-year RFI observatory: the megaconstellation trend

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: read
Pérez+2020 (SoPh 295:11) full-text, confirm no post-2019 e-Callisto RFI-trend paper exists, and
verify Space-Track TLE access

## Context

The only RFI-trend precedent on this archive stops in 2019 — the year Starlink launches began
(Pérez+2020, SoPh 295:11; the fence noted in fable-ideas' F17 entry). Dedicated-instrument
satellite-RFI studies (LOFAR arXiv:2307.02316, SKA-Low arXiv:2506.02831) do not touch archival
spectrograph records, so the 2012–2026 e-Callisto record is an unexamined 15-year RFI
observatory spanning the entire megaconstellation era. Data: the e-Callisto FITS archive —
already streamed by the merged `ecallisto_pipeline` Airflow ingest, with coverage-correction
statistics from the merged `ecallisto_census` slice; infra reuse is total. Attribution layer:
Space-Track TLEs for satellite pass windows over each station. An honest null is still a
citable spectrum-management result.

## Deliverables

- `src/jansky_research/rfitrend.py`: `occupancy_metric` (robust quantile-based per-station/
  per-channel daily occupancy, immune to solar bursts), `select_stable_stations`
  (configuration-stability screen for 10–20 stations, 2012–2026), `trend_fit` (occupancy trend
  in the documented unintended-emission bands, with change-point flags), `tle_pass_windows`
  (Space-Track fetch + pass prediction, `# pragma`), `fixed_transmitter_controls` (FM/DAB
  bands), `synthetic_spectrogram_stack`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/rfitrend/`; `survey/rfitrend-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of Pérez+2020 (method + Spanish-site result to reproduce);
   confirm no post-2019 trend paper landed; verify Space-Track registration/TLE pull; standing
   full-text novelty pass; confirm the archive tree layout the Airflow ingest already uses.
1. Tooling + synthetic recover-a-known: build synthetic spectrogram stacks with injected solar
   bursts plus an injected secular occupancy trend; the quantile metric must recover the trend
   while remaining burst-immune (burst on/off makes no difference to the occupancy series).
2. Real leg: run the occupancy metric over 10–20 configuration-stable stations, 2012–2026;
   trend the documented unintended-emission bands; attribute candidate trends to constellation
   growth via TLE pass-window coincidence; report per-station and network-differential results.
3. GATE-2 science review naming the fable-ideas caveat: station hardware changes can masquerade
   as trends — the differential in-station band comparison must gate every trend claim.
4. Paper: 15-year occupancy trends + attribution + honest per-station systematics table.

## Verification

Reproduce Pérez+2020's Spanish-site occupancy increase from the same archive; FM/DAB fixed
transmitters serve as flat controls (no trend expected); injected synthetic trend recovered;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Station hardware changes masquerade as trends** → differential in-station bands (satellite
  band vs adjacent control band at the same station) mitigate; screen with
  `select_stable_stations` and drop stations with documented receiver swaps.
- **Sparse/uneven station uptime** → reuse `ecallisto_census` coverage corrections; report
  per-station exposure alongside every trend.
