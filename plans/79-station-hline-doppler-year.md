# 79 — Earth's orbit from the roof: the annual ±142 kHz H-line Doppler sinusoid with an error budget

Status: 📋 planned (hardware-gated) — needs the commissioned H-line receiver plus plan 78's
calibrated pipeline, then a 12-month daily-spectra campaign (month-12+ deliverable; data complete
~late 2027 if first light lands late summer 2026)

## Context

Earth's ±30 km/s orbital velocity imprints a ±142 kHz annual sinusoid on the H-line frequency of
a fixed galactic pointing. The station sweep found the best amateur precedent is a *two-date*
proof of concept — nobody has published 12 months of daily spectra with a real uncertainty
budget. Rigor is the novelty axis, and this is the station's flagship measurement
(`station/operations.md`, year-1 plan). The limiting term is receiver gain/bandpass drift over a
year, not sensitivity: the design answer is the weekly calibration cadence (50 Ω load/cold sky)
and possibly a Dicke-switched second dish (`station/interferometry.md`, "alternative role for
dish B"). Depends on plan 78 (the calibrated pipeline is the substrate); sequence after it.

## Deliverables

- `src/jansky_research/hlinedoppler.py`: `barycentric_prediction` (expected line shift vs date
  for the fixed pointing, astropy), `centroid_series` (per-day line centroid + error from the
  plan-78 pipeline), `fit_annual_sinusoid` (amplitude/phase with covariance),
  `drift_budget` (gain/bandpass/TCXO/pointing terms → total), `residual_vs_hi4pi`
  (per-pointing falsifiable comparison), `synthetic_year` (injected sinusoid + drift),
  `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic year fixtures — no sky data needed for tests);
  `papers/hlinedoppler/`; `survey/hlinedoppler-findings.md`; wiring.

## Approach

0. GATE 0: plan 78 merged and the receiver producing calibrated daily spectra with the weekly
   cal cadence running per `station/operations.md`; data-continuity SLO instrumented.
1. Tooling + synthetic recover-a-known, ahead of the campaign: inject a ±142 kHz sinusoid plus
   realistic gain/bandpass drift and gaps into `synthetic_year`; the fit must recover amplitude
   and phase with honest, calibrated error bars, and `drift_budget` must attribute the noise.
2. Campaign leg (`# pragma: no cover`, 12 months): daily fixed-pointing spectra; weekly cal
   anchors the gain history; mid-campaign, evaluate the Dicke-switched second dish if drift
   dominates the budget as expected. Quarterly checkpoints (per the year-1 plan) de-risk it.
3. Analysis: measured centroid series vs `barycentric_prediction`; the falsifiable form is the
   residual vs HI4PI per pointing — the survey says what the line should do; deviations must be
   explained by the budget or reported as unexplained.
4. GATE-2 science review → paper: amplitude, phase, and the full uncertainty budget.

## Verification

Synthetic year round-trips recover injected amplitude/phase within stated errors under injected
drift; the real-data anchor is residual-vs-HI4PI per pointing; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Gain/bandpass drift is the limiting term** → weekly cal cadence from day one; Dicke-switched
  second dish as the escalation path; the drift budget is a deliverable, not an afterthought.
- **Hardware delay compresses the campaign** → the deliverable is month-12+; all tooling ships
  first, and quarterly partial fits are publishable progress markers, not the paper.
- **Data-continuity gaps** (weather, outages) → the SLO + watchdog architecture in
  `station/operations.md`; the fit tolerates gaps, the budget reports them.
- **Urban RFI near the line** → plan 78's flagging; sidereal/solar separation as the QC test.
