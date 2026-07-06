# 81 — Geometry-corrected meteor-shower radio flux vs the Global Meteor Network (Perseids/Geminids)

Status: 📋 planned (hardware-gated) — needs the meteor-scatter station built (yagi + RTL-SDR +
Echoes on the Pi, currently planned); first shower data ~4 weeks after station start, with the
Perseids (Aug 12–13, 2026) as the commissioning target

## Context

Amateur forward-scatter meteor counting is a mature hobby (RMOB aggregates monthly counts), but
the station sweep found the literature stops at raw ping rates: nobody publishes single-station
radio rates geometry-corrected to a flux and validated against the Global Meteor Network's
calibrated optical km⁻²hr⁻¹ product. Geometry correction + public-survey ground truth is the
novelty axis, matching the repo's house style. The build is `station/meteor-scatter-station.md`
(FM forward scatter; Echoes per-event logs → CSV/SQLite). The honest ML add-on is a
ping/aircraft/RFI spectrogram classifier — head-echo work is physically out of reach for forward
scatter (see fable-ideas "Corrections") and is dropped. Perseids first, Geminids (Dec) second.

## Deliverables

- `src/jansky_research/meteorflux.py`: `read_echoes_log` (per-event ingest), `classify_event`
  (ping/aircraft/RFI spectrogram classifier — simple features first, small CNN only if needed),
  `observability_function` (single-station forward-scatter geometry vs radiant position/hour),
  `corrected_rate` (rate ÷ observability → flux-proportional series), `compare_gmn` (GMN flux +
  RMOB counts alignment), `synthetic_pings` (injected pings/aircraft/RFI + known activity
  profile → round-trip), `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic/offline fixtures — no sky data needed for tests);
  `papers/meteorflux/`; `survey/meteorflux-findings.md`; wiring.

## Approach

0. GATE 0: meteor station on air per `station/meteor-scatter-station.md` — frequency selected
   (locally silent, 50–100 kW station 500–1500 km away), Echoes logging, and the GMN flux
   product + RMOB access path re-verified (standing full-text/data GATE-0 from fable-ideas).
1. Tooling + synthetic recover-a-known, ahead of hardware: `synthetic_pings` injects a known
   shower activity profile through the observability function plus aircraft/RFI contaminants;
   classifier + correction must recover the injected profile and flux normalization.
2. Commissioning leg (`# pragma: no cover`): Perseids 2026 — per-event logs, classification,
   dawn-peak sporadic sanity check (built into the station design), monthly RMOB report.
3. Science leg: geometry-corrected shower flux vs GMN's calibrated optical flux for the same
   nights (Perseids, then Geminids); sporadic baseline subtraction; mass-index caveats stated.
4. GATE-2 science review → paper: method + two-shower radio/optical flux comparison.

## Verification

Synthetic round-trip recovers the injected activity profile through the full classify+correct
chain; the real anchor is GMN's calibrated km⁻²hr⁻¹ product as ground truth (plus RMOB
consistency); checks green; GATE-2 sign-off.

## Risks & mitigations

- **Station build slips past the Perseids** → software + synthetic legs ship regardless; the
  Geminids (Dec 2026) are an equally good first shower.
- **Single-station geometry correction is the hard part** → the observability function carries
  stated model uncertainties; compare *shapes* and relative fluxes before absolute numbers.
- **Urban RFI / aircraft contamination in Philadelphia** → the classifier is a deliverable, not
  a patch; report confusion rates from labelled events.
