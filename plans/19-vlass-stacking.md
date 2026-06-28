# 19 — Sub-threshold radio stacking of a population in VLASS (with injection-recovery)

Status: ✅ done (tooling + real VLASS-SE stack + magnitude-binned radio-optical trend + paper)

## Context

Most members of an optically- or infrared-selected population are *fainter* than a radio survey's
single-source detection limit. Their *average* radio flux is still measurable by **image-plane
stacking**: at N known positions, thermal noise averages down as $N^{-1/2}$ while a coherent
sub-threshold signal adds, so the stacked image reveals the population mean (White et al. 2007; Karim
et al. 2011). The credibility of any stacked flux hinges on **calibrating the bias** — snapshot/CLEAN
and residual-vs-restored effects make raw VLASS Quick-Look stacks unreliable (Bonnassieux et al. 2023;
the VLASS QL Users' Guide) — so an **injection-recovery** step is mandatory, and the VLASS
**Single-Epoch** (self-calibrated, deeper-CLEAN) images are the right substrate, not Quick-Look.

This slice builds a small, tested image-plane stacking tool **with injection-recovery calibration**,
and applies it to a population, reporting the population-mean radio flux (or a calibrated upper limit).
It reuses the project's verified VLASS CADC-SODA cutout path (the `radio-cutout` skill) and the
`vlass` forced-photometry pattern.

## Reuse

- the `radio-cutout` skill / `vlass.fetch_vlass_cutouts` CADC-SODA path (validated: real VLASS FITS).
- `vlass.measure_image_flux` pattern (peak + annulus rms) for the per-stack measurement.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, AASTeX paper).

## Deliverables

- `src/jansky_research/stacking.py`:
  - `median_stack` — pixel-wise sigma-clipped median of N centred cutout stamps (robust to bright
    interlopers; White et al. preferred median over mean).
  - `measure_stacked_flux` — central peak + annulus rms + SNR of the stacked image.
  - `injection_recovery` — inject point sources of known flux into the (background) cutouts, stack,
    and measure the recovered/injected ratio: the **bias calibration** that makes the flux credible.
  - `synthetic_population` — N synthetic cutouts of a sub-threshold source + noise (individually
    undetected; the stack recovers the mean).
  - `fetch_se_cutouts` — VLASS Single-Epoch cutouts at a target list via CADC SODA.
  - `run(offline=...)`, `_figure` (the stacked stamp + the injection-recovery curve), `_write_macros`,
    `_main`.
- `tests/test_stacking.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** The stack + measure + injection-recovery core, validated on a synthetic
   population whose individually-undetected sources are recovered at high SNR by the stack, with the
   injection-recovery ratio near unity for the clean model.
2. **Real data (next).** A curated target population (e.g. Gaia DR3 / WISE AGN below the VLASS limit),
   VLASS-SE cutouts via CADC SODA, bright-neighbour exclusion, median stack, **injection-recovery on
   the real cutouts** to calibrate the bias, then the population-mean flux (or a 3σ upper limit).
3. **GATE-2** before write-up — the measurement survives the QL-vs-SE, CLEAN-bias, confusion, and
   bright-neighbour caveats; the injection-recovery calibration is shown.
4. **Write-up** as `papers/stacking/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the injected sub-threshold mean flux at high stacked SNR, and the
  injection-recovery ratio is near unity for the clean Gaussian model.
- `median_stack` / `measure_stacked_flux` match hand-computed values on textbook inputs.
- (Real-data, later) a calibrated population-mean flux or upper limit on VLASS-SE; GATE-2.
