# Findings — 3D triangulation of a type III source with STEREO-A + STEREO-B

`jansky_research.triangulate` locates a solar type III radio source in 3D **geometrically**, from the
two STEREO spacecraft's direction-finding, and cross-checks the result against the **independent**
plasma-frequency distance used by the `swaves`/`windwaves` slices. It is the geometric complement to
those drift-to-distance slices: where they map frequency → distance *through a density model*,
triangulation needs no density model at all — two lines of sight fix the source.

## Method

STEREO/WAVES (S/WAVES; Bougeret et al. 2008) Level-3 HFR **direction-finding** (goniopolarimetry;
Cecconi et al. 2008; Krupar et al. 2012) gives, per time and frequency, the **direction of arrival** of
the emission — the direction toward the source — as an azimuth and colatitude in the heliocentric HEEQ
frame, plus the spacecraft's HEEQ position. For each frequency we:

1. **Intensity-weight vector-average** the per-sample direction over a drift-tracking burst window on
   each spacecraft (scalar angle averaging is wrong near the azimuth wrap; vector averaging is the
   correct circular mean).
2. **Triangulate**: the source is the least-squares closest point of the two rays
   (spacecraft position + direction). We keep a channel only if both rays point *forward* (`t>0`), the
   **miss distance** (shortest segment between the rays) is below a threshold, and each spacecraft
   contributed enough good samples.
3. Compare the geometric heliocentric radius to the Leblanc plasma-frequency radius (same harmonic),
   and read off the source's heliographic **longitude and latitude** — which the drift method cannot give.

## Recover-a-known: the 2013-05-15 type III (STEREO-A + STEREO-B)

| quantity | value |
|---|---|
| spacecraft baseline (Sun-centred angle) | **82°** (A ahead, B behind) |
| channels triangulated | **38** (0.125–1.975 MHz) |
| geometric heliocentric range | **15.3 → 106.1 R⊙** (0.07 → 0.49 AU) |
| median miss distance (ray consistency) | **17.1 R⊙** |
| source longitude / latitude (HEEQ, median) | **169° / +4.7°** (near-ecliptic) |
| corr(r_geom, r_plasma) | **0.989** |
| ratio r_geom / r_plasma (median) | **2.18** |

**The headline is the correlation, not the absolute scale.** The geometric distance — built purely
from two pointing directions — and the plasma-frequency distance — built purely from a density model —
track each other in shape at r = **0.99** across two decades. Part of any such correlation is trivially
expected: both estimators decrease monotonically with frequency by construction, so even a *linear*
ramp in frequency already correlates with the Leblanc curve at r ≈ 0.75. What the 0.989 value adds is
that the geometric distances follow the **correct log–log curvature** of the density model across two
decades — a real cross-check of two fully independent distance estimators, not merely a shared trend.
The geometry additionally pins the source to HEEQ longitude ≈169°, latitude ≈+5° (near the ecliptic),
which the 1-D drift method cannot do.

The absolute geometric radii run **~2× the average-Leblanc plasma radii** (ratio 2.18). Two effects,
both honest, push the same way and we do **not** claim to separate them:
- **Outward triangulation bias from direction-finding noise.** A single type III has a large apparent
  source size (~60° FWHM; Krupar et al. 2014), so per-sample directions scatter by tens of degrees, and
  noisy near-grazing rays bias the closest-point *outward*. Quantified at the **real 82° baseline**, our
  synthetic fixture (zero model error) gives ratio ≈1.08 for 9° per-sample scatter, ≈1.30 for 18°, and
  ≈1.58 for 25° — so at the physically motivated scatter (~25°) this bias plausibly accounts for a large
  part of the 2.18, but not all of it. (The default fixture uses a *wider* 135° baseline and so is not a
  direct control for the real geometry; the 82° numbers above are the relevant ones.) The geometric
  distances are therefore upper-biased and the median miss (17 R⊙) sets the per-channel uncertainty scale.
- **Density enhancement.** 2013-05-15 was an active period (AR 11748, X-class flares); a denser-than-
  average corona/wind puts a given plasma level *farther out* than the mean Leblanc model, which drives
  the residual r_geom > r_plasma beyond the noise bias. The harmonic assumption is the favourable one
  here — fundamental emission would make r_plasma *smaller* and the ratio *worse*.

## Honest assessment & caveats

- **A reproduction/method demo, not a survey.** One well-placed event on the public L3 DF product; the
  contribution is a tested, reproducible *geometric* localisation + a clean independent cross-check of
  the density-model distance, with the noise budget surfaced.
- **The robust output is the correlation and the (longitude, latitude); the absolute radial scale is
  upper-biased** by DF noise (and possibly density enhancement) by a factor ~2 — we report it but do not
  interpret the factor as a measurement.
- **Only the low band (0.125–~2 MHz) triangulates** for this event: higher HFR channels lack DF
  solutions in the window (below the goniopolarimetry SNR / flagged), so this reaches the *interplanetary*
  source, not the coronal one.
- **Two spacecraft, one solution per frequency** — no third view, so each channel's position is a single
  closest-point estimate with the miss distance as its consistency check, not a fitted error ellipsoid.
- **Needs STEREO-B**, lost in 2014, so this exact two-view geometry is only reproducible on 2007–2014
  data. The tool runs on any A+B day.
- **Reproducible:** `python -m jansky_research.triangulate --date 20130515` regenerates the metrics,
  the (geometric-vs-plasma distance + HEEQ geometry) figure, and the macros from the public SPDF CDFs.
