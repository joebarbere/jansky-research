# 54 — e-Callisto as an accidental 15-year RFI observatory: the megaconstellation trend

Status: ✅ done 2026-07-10 — a SYSTEMATICS-LIMITED NULL. Burst-immune, gain-cancelling
narrowband-UEM-line metric over the continuous e-Callisto archive 2012–2026 (intrinsic Starlink
lines 125/135/150/175 MHz, Di Vruno+2023; GRAVES 143.05 excluded); synthetic recovers both
metrics. Real leg (286 usable station-months, HUMAIN/ALMATY/GLASGOW): the two line-sampling
stations trend in OPPOSITE signs (HUMAIN +0.45/yr rising r=0.48; ALMATY −0.13/yr falling r=0.04),
so the pipeline cross-station coherence test returns `coherent_rise=False` — no Starlink
attribution, HUMAIN flagged only for satellite-pass-gated follow-up. GATE-2 PASS (honest framing;
caught the 137→143.05/GRAVES line error, re-run). Module `rfitrend.py`, paper `papers/rfitrend/`,
findings `survey/rfitrend-findings.md`. — GATE 0 done 2026-07-10: novelty PASS. Anchor = Prieto/Pérez+2020
(SoPh 295:11): a two-epoch 2012-vs-2019 campaign at Spanish sites, ~**2× RFI rise** across
45–870 MHz, ending before Starlink scaled — a dedicated two-snapshot campaign, NOT the continuous
operational-archive time series this slice builds (methodologically distinct + extends past 2019).
No post-2019 e-Callisto RFI-trend/megaconstellation census exists (LOFAR/SKA-Low UEM studies are
dedicated-instrument, motivate not pre-empt). **Starlink UEM bands in 45–870 MHz:** broadband
110–188 MHz, **intrinsic** narrowband lines at **125/135/150/175 MHz** (Di Vruno+2023). The
143.05 MHz feature is reflected GRAVES radar (NOT intrinsic — excluded), corrected 2026-07-10 after
an initial 137.05 transcription error. **Control:** FM 87.5–108 MHz is the cleanest flat control; DAB Band III
(174–230) + UHF TV have analog-switchover/700 MHz-digital-dividend STEP changes (switchover-aware
masking; FM primary). **Attribution:** the public planet4589 Starlink-count-vs-date time series
(no Space-Track login needed). **Stations (45–870 MHz, span 2012–2026, sample UEM):** HUMAIN,
ALMATY, GLASGOW (+ BLEN7M anchor, partial UEM). Reuse: `solarbursts.fetch_ecallisto` +
`ecallisto_catalog` (archive) + `ecallisto_census.coverage_corrected_rate`. The differential
in-station UEM-vs-control comparison is the load-bearing systematics gate.

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
