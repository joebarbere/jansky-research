# Findings — inner-heliosphere type III, tracking the beam beyond the corona (Wind/WAVES)

`jansky_research.windwaves` fits a type III burst's frequency drift in a Wind/WAVES dynamic spectrum
and inverts it, via the Leblanc heliospheric density model, to the electron beam's outward radial speed
and the heliocentric distance it reaches. It is the wider-distance companion to the coronal
`solarbursts` slice (#21): the same beam, tracked beyond the Newkirk corona with a heliospheric density
model and space-based data. Reuses `solarbursts`' dynamic-spectrum pipeline and `jansky.solar`.

## Data

Wind/WAVES Level-2 radio CDFs on NASA SPDF (public, no auth; needs `cdflib`; Bougeret et al. 1995).
RAD2 covers 1.075–13.825 MHz; RAD1 reaches down to 20 kHz (≈ the local plasma frequency at 1 AU). One
file per day per receiver. This run uses RAD2.

## Recover-a-known: the 2003-10-28 X17 flare type III

The strongest emission in the RAD2 file falls at **11:06 UT**, the time of the X17 flare. Across all
256 channels (1.075–13.825 MHz), the per-channel peak times trace a high-to-low drift; mapping each
frequency through the Leblanc model places the emission at heliocentric radii:

| quantity | value |
|---|---|
| frequency span | 1.075–13.825 MHz |
| heliocentric range (harmonic) | **2.4 → 10.2 R⊙** (≈ 0.011 → 0.048 AU) |
| outward radial speed (harmonic, 2f_p) | **0.083 c** (24 900 km s⁻¹) |
| outward radial speed (fundamental) | 0.045 c |
| height–time fit R² | 0.65 |

So the beam is tracked **beyond the upper corona into the inner heliosphere** — from a couple of solar
radii to ~10 R⊙, the latter near the **Alfvén surface** (the corona/solar-wind boundary, ~10–20 R⊙).
This is *not* yet the interplanetary (super-Alfvénic, ≳20 R⊙ to 1 AU) regime — RAD1 (to ~20 kHz) would
be needed to follow the beam to 1 AU. **The harmonic assumption is a convention** (we take 2f_p, as is
usual at these frequencies, but did not search for a fundamental–harmonic pair to confirm it).

**On the speed.** The recovered 0.045–0.083 c is *low*: studies at comparable distances find faster
beams (~0.3 c at 11–30 R⊙, Fainberg et al. 1972/74; 0.17–0.35 c, Reiner & MacDowall 2015; 0.02–0.35 c
at lower frequencies, Krupar et al. 2015). The most likely causes of our underestimate are the
**peak-time** ridge sampling (the peak arrives after the onset, biasing the speed low) and the **scatter
in the height–time fit (R² = 0.65)** — *not* projection (see caveats). The synthetic fixture only
confirms the inversion is algebraically self-consistent (injected 0.1/0.15/0.2 c recovered within 15%);
it is a round-trip code check, not evidence the real value is right.

## Honest assessment & caveats

- **Projection is negligible here, despite intuition.** At 2–10 R⊙ the Parker spiral angle is only a few
  degrees, so the radial-vs-field-aligned speed correction is <1% — it does **not** explain the low
  speed. (It would matter at 1 AU, where the spiral angle is ~45°.) The real biases are peak-time
  sampling and the moderate fit.
- **An average density model on an extreme event.** The Leblanc model is the *mean* solar wind; the
  2003-10-28 X17 flare drove a >2000 km s⁻¹ halo CME, which can locally enhance the density at 2–10 R⊙
  by a factor of several — shifting the inferred radii (and so the speed) by more than the harmonic
  factor of 2. This is the dominant systematic for this particular event.
- **Harmonic/fundamental factor ~2.** 0.045 c (fundamental) vs 0.083 c (harmonic); we report both. IP
  type III emission at these frequencies is commonly harmonic, so 0.083 c is the primary value.
- **The fit is moderate (R² = 0.65).** Near the Sun the band compresses (the beam crosses 2–3 R⊙ fast,
  so the high-frequency channels peak almost together) and the per-channel peak times carry scatter; the
  height–time track is real but not a tight line.
- **RAD2 only, single event.** A reproduction/method demo on the well-documented X17 benchmark, not a
  survey; combining RAD1 would extend the track toward 1 AU.
- **Reproducible:** `python -m jansky_research.windwaves --date 20031028 --receiver rad2` regenerates
  the metrics, the dynamic-spectrum + beam-track figure, and the macros from the public SPDF CDF.

## Full referee round (2026-08-25): MAJOR REVISION, 12 findings, one BLOCKER

The referee recovered the full 256-point ridge from the committed figure's marker coordinates
and re-fit it, reproducing EVERY committed metric to the printed precision -- nothing is
fabricated or stale. The problems are propagation and framing.

**BLOCKER: the self-declared "dominant systematic" is never propagated, and both abstract
framings flip inside it.** Scaling the Leblanc density by the "factor of several" the paper
itself invokes: x2 gives 0.119 c / 14.0 R_sun; x4 gives 0.173 c / 19.5 R_sun -- INSIDE the
reiner2015 range the paper says it falls below, and AT the >=20 R_sun "interplanetary regime"
the abstract says RAD2 cannot reach. And harmonic=1 with x4 density is numerically IDENTICAL
to harmonic=2 with x1 (f_p^2 ~ n), so the paper's two caveats are one axis. Fix: commit the
(harmonic x density-scale) grid and quote speed and reach as brackets.

**MAJORs:** no uncertainty anywhere while the honest error is ~2.5x the naive one -- the 256
points occupy 10 time columns (197 in the first three); column jackknife SE = 0.0138 c,
drop-last-column moves the speed -13% -- quote 0.083 +/- 0.014 c; the headline reach
(10.25 R_sun, "near the Alfven surface", uncited) rests on ONE ridge point in an isolated time
column and is really the band edge; the stated cause of R^2=0.651 ("band compression") is
contradicted by the paper's own forward model (0.995 +/- 0.001 at the real band and cadence)
-- the residuals have monotonic structure (apparent local speeds 0.034 -> 0.16 c across the
band), i.e. deceleration-like curvature, not compression; the "low speed" is
estimator-dependent (OLS r|t 0.083; t|r 0.128; TLS 0.103; column means 0.100) and the
published choice is the minimum of the family; "0.045 c" (abstract+Results) and "11:06 UT"
are in no committed evidence in a paper claiming pipeline provenance; no ridge CSV and no
analysis parameters committed (pad_s=1200, snr=5, clip sigma).

**MINOR:** the inherited unconverged _robust_linfit is DORMANT here (replayed: all-True mask
at every iteration; converge=True identical) -- but only the referee's replay proves it, so
pass converge=True and record clip provenance; the synthetic fixture runs at 3 s cadence over
a band to ~50 R_sun (5.4x finer than real) and its 15% tolerance is trivially met -- a
matched-sampling fixture shows cadence is NOT the R^2 story, which strengthens finding 4;
krupar2015 author list wrong (Kontar 2nd per Crossref -- same defect as swaves);
reiner2015 title drops "Solar"; "~0.3 c historically" uncited (Fainberg not in refs.bib); the
flare association is day-max-near-flare, guaranteed on an X17 day -- state the GOES offset and
the count of other >5-sigma ridges.

**NIT:** \wwTruth/\wwRatio placeholders unused; \wwFlo/Fhi presented as RAD2 coverage but
derived from the ridge; wind speed for the projection claim unstated.

**Status: fixes pending** (all fixes run on the existing ridge; no new download needed).

**Status: RESOLVED (2026-08-25).** The dominant systematic is propagated and the estimator is
rebuilt on the independent unit.

The (mode x density-scale) grid is committed: 0.054 c / 6.0 R_sun (fundamental, 1x) through
0.102 c / 10.25 (harmonic, 1x) to 0.215 c / 19.5 R_sun (harmonic, 4x) -- matching the
referee's own sweep, with the exact f_p^2~n degeneracy stated in code and paper, and the
abstract now quotes the bracket and claims neither distance regime exclusively. The headline is
a ONE-POINT-PER-TIME-SAMPLE (column median) fit -- the referee's estimator-family analysis made
the per-point OLS the family minimum, and en route we measured its converged-clip variant
collapse outright on the sibling STEREO ridge -- giving 0.102 +/- 0.011 c (leave-one-column-out
jackknife over the 10 samples), with the all-points fit (0.083) and inverse regression (0.122)
committed as the bracket. R^2 = 0.835 on the columns, and the paper now says the residuals'
band-monotonic structure is unexplained (deceleration-like curvature or density mismatch)
rather than blaming band compression, which the matched-cadence fixture disproves. The ridge is
committed (results/windwaves_ridge.csv) with pad/snr provenance in the metrics; the burst peak
epoch (2003-10-28T11:06:29) is committed and the flare association stated as location, not
validation; the fundamental values come from the grid macros (no hand-typing); the reach is
described as the band edge through the model, carried by one channel; the Parker correction is
the path-weighted 1-3%; macros are namespaced wwSyn*/wwReal*; the figure draws the exact quoted
fit with per-channel points, per-sample medians, and the fitted line distinguished; krupar2015's
author list and reiner2015's title fixed against Crossref; the uncited "0.3 c historically"
dropped; the fixture gained matched-cadence variants whose quantisation bias is asserted in CI.
