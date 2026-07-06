# 84 — Calibrated 21-cm quiet-Sun monitor: daily solar flux vs F10.7/RSTN from a rooftop dish

Status: 📋 planned (hardware-gated) — needs the commissioned H-line receiver in daily drift mode
plus a working Cas A/Cyg A transit flux calibration; without the flux calibration the slice is
blog-grade and is killed (that is the explicit criterion)

## Context

The Sun transits the fixed dish's beam daily for free (`station/operations.md`), and the station
sweep found many amateur "solar SNR trend" plots — none flux-calibrated, so none comparable to
the professional record. The rigor axis: tie the daily transit amplitude to a real flux scale via
Cas A/Cyg A transits (known flux densities, including Cas A's secular decay), then track the
quiet-Sun 1420 MHz flux against F10.7 (DRAO) and RSTN through the solar cycle. Only the
calibrated version is worth doing — uncalibrated it is a dashboard, not a paper, and that is the
kill criterion. A flare-coincident burst caught in a transit is an opportunistic bonus, not a
deliverable. Depends on plan 78's pipeline; data accumulates alongside plans 79/80.

## Deliverables

- `src/jansky_research/quietsun.py`: `fit_transit` (drift-transit template fit → amplitude +
  error), `flux_scale` (Cas A/Cyg A transits + adopted flux densities with secular-decay
  correction → K-to-Jy factor with uncertainty), `solar_series` (daily calibrated flux with
  error bars and cal-epoch provenance), `compare_indices` (F10.7/RSTN correlation + scale
  check), `synthetic_transits` (known fluxes through a known beam/gain drift → round-trip),
  `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic/offline fixtures — no sky data needed for tests);
  `papers/quietsun/`; `survey/quietsun-findings.md`; wiring.

## Approach

0. GATE 0: plan 78 merged; receiver in daily drift mode; demonstrate that Cas A and/or Cyg A
   transits are detectable at usable SNR in stacked passes — if neither calibrator is
   recoverable, the flux scale fails and the slice is killed (blog-grade otherwise, by design).
1. Tooling + synthetic recover-a-known, ahead of hardware: `synthetic_transits` pushes known
   solar + calibrator fluxes through a known beam and injected gain drift; `flux_scale` +
   `solar_series` must recover the injected fluxes and track the drift via the calibrators.
2. Calibration leg (`# pragma: no cover`): stacked Cas A/Cyg A transits → K-to-Jy factor with
   stated uncertainty, re-derived on the weekly cal cadence from `station/operations.md`.
3. Monitoring leg: months of daily solar transits → calibrated flux series; compare level and
   variability against F10.7 and RSTN 1415 MHz; any flare-coincident burst is reported as an
   opportunistic aside, not claimed as a monitoring capability.
4. GATE-2 science review → paper: an amateur flux-calibrated quiet-Sun series with its scale
   uncertainty — the artifact the SNR-trend genre lacks.

## Verification

Synthetic round-trip recovers injected solar/calibrator fluxes under gain drift; the real anchor
is the Cas A/Cyg A flux calibration and agreement of the series' scale with F10.7/RSTN within
the stated uncertainty; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Calibrators too weak for this dish** → GATE-0 tests exactly this with stacked transits; a
  failed gate kills the slice cheaply and the finding is documented.
- **Without flux calibration it is blog-grade** → that IS the kill criterion; no uncalibrated
  fallback deliverable is planned.
- **Gain drift between calibrator and solar transits** → weekly cal cadence + per-epoch scale
  factors with provenance; drift enters the error bars, not the residuals.
- **Urban RFI / hardware delay** → transit fits reject non-beam-shaped events; software ships first.
