# 72 — e-Callisto fine-structure census: network-scale spikes, J-bursts, and U-bursts

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

Solar radio fine structure — spikes, J-bursts, U-bursts — is studied on single high-cadence
instruments (the LOFAR spike-pair result, Nat. Comm. 2026, is the fresh anchor); no
network-scale census exists across the e-Callisto archive (fable-ideas F35). This repo already
streams that archive (`airflow/` e-Callisto ingest) and owns the census statistics and
coverage-correction machinery from the merged `ecallisto_census` slice — the new piece is a
sub-2-s fine-structure detector run network-wide. e-Callisto cadence (typically 0.25 s) sits at
the edge of the spike regime, so cadence-limited completeness is reported honestly as a design
constraint, not buried. The payoff is an independent-instrument, multi-station occurrence table
that the single-instrument literature cannot produce.

## Deliverables

- `src/jansky_research/finestructure.py`: `fetch_spectrogram` (archive tree, `# pragma`; ingest
  infra reused), `detect_fine` (sub-2-s detector: short-duration/narrow-band ridge + drift-sign
  classification into spike/J/U morphologies), `coverage_correct` (per-station duty-cycle
  correction, `ecallisto_census` reuse), `census_table` (occurrence vs time/frequency/station),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/finestructure/`; `survey/finestructure-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of the LOFAR spike-pair paper and a fine-structure-census
   literature sweep to confirm no network-scale e-Callisto census exists; verify the current
   archive tree layout (it has moved before — see the opportunity-scan corrections).
1. Tooling + synthetic recover-a-known: injected spikes/J/U shapes at known rates into real
   quiet spectrograms; the detector recovers rates and morphology labels, and the
   cadence-completeness curve is measured from the same injections.
2. Real leg: run over configuration-stable stations across a bounded interval spanning quiet and
   active periods; coverage-corrected occurrence census; independent-instrument comparison of
   spike statistics against the published LOFAR result on overlapping dates/bands.
3. GATE-2 (RFI masquerading as spikes — multi-station coincidence vetting; cadence-limited
   completeness stated as the headline caveat) → paper (census + methods note).

## Verification

Injection-measured completeness curve published alongside every rate; spike statistics compared
to the LOFAR spike-pair result as the independent-instrument anchor, with the cadence gap stated
plainly; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Cadence-limited completeness** (0.25 s vs ms-scale spikes) → measured by injection and
  reported honestly; the census claims only the ≥cadence-resolvable population.
- **RFI at spike timescales** → multi-station coincidence and per-station RFI controls
  (`ecallisto_census` discipline); single-station-only events flagged, not counted.
- **Heterogeneous station hardware** → restrict to configuration-stable stations; report
  per-station, not just pooled, statistics.
