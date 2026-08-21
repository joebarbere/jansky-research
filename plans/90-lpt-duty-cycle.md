# 90 — `lptduty`: how often are long-period transients actually on?

Status: 📋 planned 2026-08-20. Sourced from Rose et al. 2026 (Nature Astron. 10, 1166;
`10.1038/s41550-026-02882-x`), not from `fable-ideas.md`.

## Context

Rose et al. classify ASKAP J174508.9−505149 as an accreting white-dwarf binary and report,
qualitatively, that its bursts "turn off for several hours at a time". That statement is
made for one source, from one campaign, and is not quantified anywhere in the literature
for the class.

This repo already holds the data to quantify it. `results/lptv_vast_epochs.csv` (committed
evidence from the `lptv` VAST sweep) carries **966 snapshots across the 10 VAST-covered
LPTs**, each with `epoch_mjd`, `duration_s` (~12 min), forced `i_mjy`/`v_mjy` and their
per-epoch errors — plus the adjudicated detections (`results/lptv_vast_adjudication.json`).
Periods are in `data/lpt_sample.csv`. No new data acquisition is required; this is an
analysis slice over material already fetched, vetted and committed.

## What is actually measurable (and what is not)

A snapshot of length *T* on a source of period *P* with emitting-phase width *w* detects a
pulse with probability ≈ min(1, (w + T)/P) **if the source is active at all**. LPTs are also
known to switch off for hours to months. So a detection rate constrains the **product**

    f_active × (w + T)/P

and the two factors are **not separately identifiable** from detection counts alone. The
slice must say so plainly, quote the product, and separate the factors only where a
published ephemeris lets snapshots be assigned pulse phases (the `lptv` slice already did
this for J183950.5−075635, so the machinery exists).

**A null divides by sensitivity, not by epoch count** (the `frblens` lesson, and the reason
that slice's limit moved by 4×). Each epoch has its own noise, and a pulse of a given flux
is not detectable in all of them. The per-source denominator is therefore
Σ_epochs ε_i(S), the summed detection efficiency at an assumed pulse flux, **not** the
number of epochs. Epochs with ε ≈ 0 constrain nothing and must not inflate the denominator.
Report the efficiency-weighted exposure alongside the raw epoch count.

## Deliverables

- `src/jansky_research/lptduty.py` — tested pure logic: per-epoch efficiency from the
  committed noise, the binomial/Poisson constraint on the product above, phase assignment
  where an ephemeris exists, and an offline synthetic fixture. Network: none.
- `results/lptduty_metrics.json` — per source: N epochs, Σ duration, efficiency-weighted
  exposure, detections, the constrained product with a confidence interval (upper limit
  where zero detections), and an explicit `identifiable: false` flag on f_active where
  phase information is absent.
- `survey/lptduty-findings.md`, and a figure: efficiency-weighted exposure vs constraint
  per source.
- Write-up: most likely an **RNAAS note** (the genre fits a single quantitative statement
  over existing data), pending the size of the result.

## GATE 0 — before any code

1. **Novelty.** Full-text search for an existing LPT duty-cycle/activity-fraction census.
   Rose et al. and the discovery papers report per-source behaviour; a uniform multi-source
   constraint may already exist in a VAST or MWA transient paper. If it does, stop.
2. **Ephemeris audit.** For each of the 10 sources, record whether a published ephemeris of
   sufficient precision exists to phase snapshots years after the reference epoch. Phase
   drift over a multi-year baseline is the limiting systematic and must be propagated, not
   assumed away — `lptv` round 3 already had to fix exactly this for the 2026 epochs.
3. **Selection honesty.** VAST pointings are not randomly phased with respect to any LPT
   period; check for aliasing between the survey cadence and each period before quoting a
   probability that assumes independent sampling.

## Increments

1. Efficiency + constraint logic with tests (offline fixture).
2. Run over the committed CSV; write metrics + findings.
3. Phase-resolved leg for sources with usable ephemerides.
4. GATE-2 science review, then the note.

## Related

`lptv` (the source data and the phase machinery), `wdpulsar` (the sensitivity-weighted null
it inherits its method from), `frblens` (the lesson that sets the denominator).
