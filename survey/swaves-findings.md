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
| outward radial speed (harmonic, 2f_p) | **0.150 c** (45 100 km s⁻¹) |
| outward radial speed (fundamental) | 0.075 c (to 0.19 AU) |
| height–time fit R² | **0.906** |

The beam is tracked from the corona out to **~0.38 AU — about halfway to Earth, genuinely
interplanetary** (well beyond the Alfvén surface). This is exactly the regime the Wind/WAVES slice
flagged as out of reach. Two things improve markedly over that RAD2-only run: the speed is **0.150 c
(harmonic), squarely in the canonical type III range** (0.1–0.5 c), rather than the low 0.083 c; and the
much wider distance baseline gives a **tight fit, R² = 0.906** (vs 0.65 for RAD2). The synthetic fixture
confirms the inversion round-trips (injected speed recovered within 15%).

## Honest assessment & caveats

- **A reproduction/method demo, not a survey.** One well-placed STEREO-A type III, on the public
  one-minute HFR ASCII product; the contribution is a tested, reproducible pipeline that reaches true
  interplanetary distances.
- **Harmonic/fundamental factor ~2.** 0.075 c (fundamental, to 0.19 AU) vs 0.150 c (harmonic, to
  0.38 AU); IP type III emission is commonly harmonic, so we take 0.150 c / 0.38 AU as primary, but did
  not search for a fundamental–harmonic pair to confirm the mode.
- **Average density model.** The Leblanc model is the *mean* solar wind; a specific event's density can
  differ, shifting the inferred radii. (This event is less extreme than the X17+CME used in the
  Wind/WAVES slice, so the average model is a better approximation here.)
- **Peak-time, radial speed.** The ridge uses intensity peaks (biases the speed low vs onset), and the
  result is the *radial* speed (a lower bound on the field-aligned speed — though at these larger radii
  the Parker-spiral projection, ~10–20° by 0.4 AU, is no longer negligible, unlike the inner heliosphere).
- **Reproducible:** `python -m jansky_research.swaves --date 20130515 --spacecraft a` regenerates the
  metrics, the dynamic-spectrum + beam-track (in AU) figure, and the macros from the public SPDF ASCII.
