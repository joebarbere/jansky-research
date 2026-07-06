# 80 — Urban drift-scan HI strip vs HI4PI: beam-forward-modelled residuals, not "look, a map"

Status: 📋 planned (hardware-gated) — needs the commissioned H-line receiver in fixed-pointing
drift mode plus plan 78's calibrated pipeline; months of stacked sidereal passes on one
declination strip before the comparison is meaningful

## Context

Amateur drift-scan HI maps are common; the station sweep found they are uniformly "look, a map"
write-ups — none forward-models the beam before attributing residuals against a survey. The
rigorous version: pick one well-characterized declination strip, stack months of sidereal passes
(`station/operations.md` runs drift-scan continuously alongside the Doppler campaign), convolve
HI4PI with the measured dish beam, and only then interpret the residual — separating beam and
calibration systematics from anything real (and from Philadelphia's RFI). Rigor (survey ground
truth + forward modelling) is the novelty axis. Depends on plan 78's pipeline; runs in parallel
with plan 79 on the same hardware since drift-scan data accumulates for free.

## Deliverables

- `src/jansky_research/driftscan.py`: `stack_passes` (sidereal-day alignment + robust stacking
  of the strip), `beam_model` (parametric beam from drift transits of a bright source, with
  uncertainty), `forward_model_strip` (HI4PI cube → beam-convolved predicted strip),
  `residual_strip` (data − model with per-bin errors), `synthetic_strip` (known sky × known
  beam + drift/RFI → round-trip), `run`/`_figure`/macros.
- Tests to the 85% floor (synthetic/offline fixtures — no sky data needed for tests);
  `papers/driftscan/`; `survey/driftscan-findings.md`; wiring.

## Approach

0. GATE 0: plan 78 merged; receiver running fixed-elevation drift mode with the capture service
   logging averaged spectra; declination strip chosen (Galactic-plane crossing, bright anchor).
1. Tooling + synthetic recover-a-known, ahead of hardware: build the HI4PI strip extraction and
   the forward model offline; a synthetic sky convolved with a known beam plus injected gain
   drift and RFI must round-trip through stack → forward-model → flat residual.
2. Beam characterization leg (`# pragma: no cover`): measure the beam from Sun (or Cas A)
   transits; fit `beam_model` with honest width/sidelobe uncertainties.
3. Accumulation leg: months of sidereal stacking on the strip; sidereal-vs-solar separation
   flags human-clock RFI (the built-in QC from `station/operations.md`).
4. Comparison: beam-convolved HI4PI vs stacked strip; residual map with attributed error terms
   (beam, gain, RFI, survey resolution). GATE-2 → paper: a systematics note the genre lacks.

## Verification

Synthetic strip round-trips to a flat residual at stated noise; the real anchor is the
beam-forward-modelled HI4PI comparison on one strip with per-bin errors; checks green; GATE-2
sign-off.

## Risks & mitigations

- **Hardware delay** → all modelling/stacking software + synthetic legs ship first; sky data
  accumulates passively once plan 79's campaign starts, so marginal roof cost is near zero.
- **Urban RFI in Philadelphia** → sidereal/solar separation + per-pass robust statistics; report
  the excised fraction.
- **Beam knowledge dominates the residual** → measure, don't assume, the beam (step 2); if beam
  uncertainty swamps the residual, that bound is itself the honest headline.
