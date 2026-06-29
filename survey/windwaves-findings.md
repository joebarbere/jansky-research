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
