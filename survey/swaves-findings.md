# Findings — type III to true interplanetary distances (STEREO/WAVES HFR)

`jansky_research.swaves` fits a type III burst's frequency drift in a STEREO/WAVES dynamic spectrum and
inverts it, via the Leblanc heliospheric density model, to the beam's outward speed and the heliocentric
distance it reaches. It completes what the Wind/WAVES slice (#24) could not: that slice's RAD2 stopped at
1 MHz (~10 R⊙, the Alfvén surface), whereas STEREO/WAVES' HFR reaches **0.125 MHz**, tracking the beam
to **genuinely interplanetary** (super-Alfvénic, ≳20 R⊙) distances — ~0.4 AU. Reuses the Wind/WAVES
Leblanc tooling and the `solarbursts` dynamic-spectrum pipeline; the only new code is the STEREO ASCII
parser.

## Data

STEREO/WAVES (SWAVES; Bougeret et al. 2008) one-minute-averaged HFR dynamic spectra (0.125–16 MHz) on
NASA SPDF, as daily **ASCII** (public, no auth, **no CDF** — unlike Wind/WAVES). One file per day per
spacecraft (STEREO-A/B). This run uses STEREO-A.

## Recover-a-known: the 2013-05-15 type III

From the May 2013 active period (AR 11748, X-class flares), STEREO-A/WAVES HFR. Across all 319 channels
(0.125–16.025 MHz, 310 kept after sigma-clipping) the drift ridge maps through the Leblanc model to:

| quantity | value |
|---|---|
| frequency span | 0.125–16.025 MHz |
| heliocentric range (harmonic) | **2.3 → 82.6 R⊙** (0.011 → **0.38 AU**) |
| outward radial speed (harmonic, 2f_p) | **0.150 c** (45 100 km s⁻¹), *mean over 2–82 R⊙* |
| outward radial speed (fundamental) | 0.075 c (to 0.19 AU) |
| height–time fit R² | 0.906 (over 310 ridge points but only **11 independent 1-min time bins**) |

The beam is tracked from the corona out to **~0.38 AU (≈0.4 of the way to Earth), genuinely
interplanetary** (well past the Alfvén-surface region, ~10–20 R⊙) — exactly the regime the Wind/WAVES
slice flagged as out of reach. The recovered **0.150 c is a *mean* radial speed averaged over 2–82 R⊙**,
not a speed at any one radius: IP type III beams **decelerate** (Krupar et al. 2015 measure
−12 km s⁻² over 0.1–1 MHz), so the true trajectory is concave — faster near the corona (~0.3–0.5 c) and
slower far out (~0.05–0.1 c) — and 0.150 c sits at the high end of the purely-interplanetary range
(median 0.09–0.16 c harmonic, Krupar 2015) precisely because the fit also spans the faster coronal part.
What the wide 0.125–16 MHz baseline robustly constrains is the **slope (the mean speed)**, *not* the
R²: the latter is high mostly because the points span two decades in distance, and it rests on only 11
unique time samples (the 1-min cadence bins ~280 high-frequency channels into ~2 timestamps), so it is
an indicative goodness-of-fit, not evidence of 310 independent measurements. The synthetic fixture
confirms the inversion round-trips for a *constant-speed* beam (within 15%); it does not test robustness
to deceleration.

## Honest assessment & caveats

- **A reproduction/method demo, not a survey.** One well-placed STEREO-A type III, on the public
  one-minute HFR ASCII product; the contribution is a tested, reproducible pipeline that reaches true
  interplanetary distances.
- **Harmonic/fundamental factor ~2.** 0.075 c (fundamental, to 0.19 AU) vs 0.150 c (harmonic, to
  0.38 AU); IP type III emission is commonly harmonic, so we take 0.150 c / 0.38 AU as primary, but did
  not search for a fundamental–harmonic pair to confirm the mode.
- **Average density model.** The Leblanc model is the *mean* solar wind; a specific event's density can
  differ, shifting the inferred radii. This event is less extreme than the X17+halo-CME used in the
  Wind/WAVES slice, so the average model is a *better* (not perfect) approximation — X-class events can
  still drive factor-~2 density excursions.
- **A burst-storm period.** May 2013 (AR 11748) produced several X-class flares on May 13–15, so the
  spectrum could in principle hold overlapping bursts; the clean single-ridge fit (310 of 319 channels
  kept after sigma-clipping, no double ridge) indicates one dominant drift track, but this is the main
  qualitative uncertainty.
- **Peak-time, radial speed.** The ridge uses intensity peaks (biases the speed low vs onset), and the
  result is the *radial* speed. At 0.4 AU the Parker-spiral angle is **~16–23°** (for 600–400 km s⁻¹
  wind), so the field-aligned speed exceeds the radial value by ~6–9% — a modest, no-longer-negligible
  correction (it *was* negligible at the inner-heliosphere radii of the Wind/WAVES slice).
- **Reproducible:** `python -m jansky_research.swaves --date 20130515 --spacecraft a` regenerates the
  metrics, the dynamic-spectrum + beam-track (in AU) figure, and the macros from the public SPDF ASCII.

## Full referee round (2026-08-25): MAJOR REVISION, 12 findings, one BLOCKER

The referee reproduced the Leblanc inversion to the digit and forward-modelled the
fundamental-mode value independently; nothing is numerically fabricated. The revision is about
uncertainty, provenance, and claim strength.

**BLOCKER: the headline 0.1503 c has no uncertainty of any kind and was produced by the
non-converged legacy fit** -- and here the exposure is LIVE (n_used 310 < n_ridge 319: the mask
changed, so slope, mask and R^2 come from two different iterations). The paper also asserts
"the much wider baseline substantially reduces the slope uncertainty" without computing one.
Fix: converge=True; jackknife over the 11 INDEPENDENT time bins (not the 310 points -- the
rmstructure wrong-unit trap); quote 0.150 +/- sigma.

**MAJORs:** "0.075 c / 0.19 AU (fundamental)" is hand-typed halving with no committed
harmonic=1 run, in a paper whose prose claims every number is pipeline-written -- the referee
verified it is numerically right to 2% here, but the sibling's committed grid shows the halving
assumption CAN fail badly (solarbursts ratio 0.689), so commit a speed_grid and macro it; the
title's "to 0.4 AU" is a deterministic function of the band edge + harmonic + density model --
any burst reaching the bottom HFR channel returns 82.6 R_sun -- so qualify the title and state
that only the SPEED is event-specific; 18 of 319 channels carry 92.5% of the fit variance
(3 channels carry 2/3) and a 3-sigma clip cannot reject a high-leverage point -- report the
leverage and jackknife the low-frequency channels; the evidence block is 13 keys vs the
refereed sibling's full set (no ridge CSV, no snr/clip/pad/converged, pad_s=2400 undocumented);
the findings file's self-declared main uncertainty (storm-period burst confusion) never reached
the paper and the 9/319-clip evidence offered for it cannot do that work; the offline fixture
runs at 7.5 s cadence vs the real 60 s -- quantisation is a measured -3.9% BIAS at the real
cadence, outside the fixture's reach -- and the recovery ratio reaches no committed evidence
(unnamespaced macros structurally unfillable).

**MINOR/NIT:** the Parker-spiral correction is wrong twice (4-9% not 6-9%; and the
path-weighted track mean is 1.3-2.8%, ~3x smaller than the endpoint value applied);
krupar2015's author list is wrong (Kontar is 2nd, per Crossref); "spans all 319 channels" is
unverifiable (n_channels_total unrecorded) and 100% retention at 5 sigma over an 80-minute
argmax window is the stokesv_discovery searching-measurement shape -- commit per-channel peak
SNR; "None of these undermine..." contradicts its own list item (ii); \swHarmonic unused;
astropy cited but not imported; hardcoded 16 in the band prose; stale arxiv-submission.

**Status: fixes pending.**
