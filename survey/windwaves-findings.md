# Findings — interplanetary type III, tracking the beam into the heliosphere (Wind/WAVES)

`jansky_research.windwaves` fits an interplanetary type III burst's frequency drift in a Wind/WAVES
dynamic spectrum and inverts it, via the Leblanc heliospheric density model, to the electron beam's
outward speed and the heliocentric distance it reaches. It is the interplanetary companion to the
coronal `solarbursts` slice (#21): same beam, but tracked *out of* the corona into the solar wind with
a heliospheric (not Newkirk) density model and space-based data. Reuses `solarbursts`'
dynamic-spectrum pipeline and `jansky.solar`.

## Data

Wind/WAVES Level-2 radio CDFs on NASA SPDF (public, no auth; needs `cdflib`). RAD2 covers
1.075–13.825 MHz; RAD1 reaches down to 20 kHz (≈ the local plasma frequency at 1 AU). One file per day
per receiver.

## Recover-a-known: the 2003-10-28 X17 flare type III

The strongest emission in the RAD2 file falls at **11:06 UT**, the time of the X17 flare. Across all
256 channels (1.075–13.825 MHz), the per-channel peak times trace a high-to-low drift; mapping each
frequency through the Leblanc model places the emission at heliocentric radii:

| quantity | value |
|---|---|
| frequency span | 1.075–13.825 MHz |
| heliocentric range (harmonic) | **2.4 → 10.2 R⊙** (≈ 0.011 → 0.048 AU) |
| outward speed (harmonic, 2f_p) | **0.083 c** (24 900 km s⁻¹) |
| outward speed (fundamental) | 0.045 c |
| height–time fit R² | 0.65 |

So the beam is tracked **out of the corona into interplanetary space** — from a couple of solar radii to
~10 R⊙ within RAD2's band (RAD1 would follow it on toward 1 AU). The inferred radial speed
(0.05–0.08 c) sits at the low end of the canonical type III range (0.1–0.5 c), which is expected for an
*interplanetary* radial-speed inference (see caveats). The synthetic fixture independently confirms the
inversion is algebraically self-consistent (injected 0.1/0.15/0.2 c recovered within 15%).

## Honest assessment & caveats

- **Radial speed is a lower bound on the beam speed.** We measure the *radial* outward speed; the beam
  travels along the (Parker-spiral-winding) magnetic field, so the field-aligned speed is higher, and
  projection makes the inferred value an underestimate — part of why 0.05–0.08 c is below the coronal
  type III norm.
- **Density model + harmonic span a factor ~2.** Fundamental vs harmonic emission moves the speed from
  0.045 to 0.083 c; the Leblanc model is one of several heliospheric profiles. We quote both, not one.
- **The fit is moderate (R² = 0.65).** Near the Sun the band compresses (the beam crosses 2–3 R⊙ fast,
  so the high-frequency channels peak almost together), and the per-channel peak times carry intrinsic
  scatter; the height–time track is real but not a tight line. Peak-time (not onset) sampling further
  biases the speed low, as in the coronal slice.
- **RAD2 only — the inner-heliosphere segment.** This run uses RAD2 (1–14 MHz ≈ 2–10 R⊙); combining
  RAD1 (down to 20 kHz) would extend the track to ~1 AU. A reproduction/method demo, not a survey.
- **Reproducible:** `python -m jansky_research.windwaves --date 20031028 --receiver rad2` regenerates
  the metrics, the dynamic-spectrum + beam-track figure, and the macros from the public SPDF CDF.
