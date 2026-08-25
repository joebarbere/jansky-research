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
| source longitude / latitude (HEEQ, median) | **169° / +4.7°** — superseded: wrap-safe median is **179.5° ± 7.6** |
| corr(r_geom, r_plasma) | **0.989** |
| ratio r_geom / r_plasma (median) | **2.18** |

**[SUPERSEDED 2026-08-25 by the round-8 referee — see below.]** This section originally
presented r = 0.99 as the headline and claimed the geometric distances "follow the correct
log-log curvature" across "two decades"; the round-8 review showed the correlation cannot fail
(any monotone power law scores >0.95; a bare 1/f matches the geometric points as well as
Leblanc), the log-log slopes actually DIFFER (−0.60 vs −0.93 — the curvature claim is falsified
by the same data), the band is 1.2 decades, and the "≈169°" longitude was a scalar median taken
across the ±180° wrap (wrap-safe: 179.5°). The honest headline is the ADDITIVE comparison: the
triangulated distances sit a constant 13.4 ± 3.3 R⊙ outside the Leblanc level (slope 1.12),
i.e. one ~4° pointing bias, after which the model reproduces the geometry to 3.4 R⊙ rms.

[Superseded 2026-08-25: the block below reads the discrepancy as a ratio and invokes a density
enhancement; the committed channels show a constant ADDITIVE offset with OLS slope ~1, which
leaves no room for an enhancement — see the round-8 sections below.]
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
- **The robust output is the source direction and the additive comparison** (superseded wording:
  this bullet previously named the correlation, which the round-8 review showed cannot fail; the
  absolute radial scale carries one constant ~4° pointing-bias offset, measured by the committed
  calibration grid, and is not read as calibrated).
- **Only the low band (0.125–~2 MHz) triangulates** for this event: higher HFR channels lack DF
  solutions in the window (below the goniopolarimetry SNR / flagged), so this reaches the *interplanetary*
  source, not the coronal one.
- **Two spacecraft, one solution per frequency** — no third view, so each channel's position is a single
  closest-point estimate with the miss distance as its consistency check, not a fitted error ellipsoid.
- **Needs STEREO-B**, lost in 2014, so this exact two-view geometry is only reproducible on 2007–2014
  data. The tool runs on any A+B day.
- **Reproducible:** `python -m jansky_research.triangulate --date 20130515` regenerates the metrics,
  the (geometric-vs-plasma distance + HEEQ geometry) figure, and the macros from the public SPDF CDFs.

## A citation can pass a coordinate check and still be the wrong paper (2026-08-23)

`krupar2014` was **internally consistent**: DOI `10.1007/s11207-014-0522-x` matched its recorded
volume 289, issue 8, pages 3121--3135 exactly. The Crossref adjudication procedure that catches
constructed identifiers would have passed it. But the claim it was cited for is not in that paper.

Krupar et al. published two same-year Solar Physics companions sharing a title stem:

| | 0522-x (was cited) | 0601-z (correct) |
|---|---|---|
| subtitle | Radio Flux Density Variations with Frequency | Goniopolarimetric Properties and Radio Source Locations |
| vol/iss/pages | 289 / 8 / 3121--3135 | 289 / 12 / 4633--4652 |
| band | 125 kHz -- 16 MHz | **125 kHz -- 2 MHz** |
| source size | not measured | "apparent source size gamma is very extended (~60 deg) for the lowest analyzed frequencies" |

The paper cites it for "apparent source size is ~60 deg FWHM". That is 0601-z's result, and its
band is *exactly* this paper's `\triFlo`--`\triFhi` (0.125--1.975 MHz), so the correct citation is
also the better-matched one. **When two same-year companions share a title stem, the discriminator
is the subtitle and the claim, not the identifier** --- coordinates alone cannot separate them.

Second defect in the same sentence: 0601-z reports ~60 deg **at the lowest analyzed frequencies**,
with gamma expanding linearly with radial distance below 1 MHz. The paper stated it as a constant
across the band. Now qualified.

## The miss-distance cut was never stated (2026-08-23)

The Method said a channel is kept when the miss distance "is below a threshold" and never gave it.
It is `max_miss_rsun = 60.0` R_sun, a bare default in `triangulate.py`. That matters here because
the distances being measured are 15.3--106.1 R_sun: at the high-frequency end the cut admits
channels whose two rays miss by four times the inferred distance, and the committed `\triMiss`
= 17.1 R_sun is the *median* miss, larger than the closest source distance. The value is now stated,
with the note that it is permissive rather than selective.

**Outstanding:** sweeping `max_miss_rsun` (15/30/60/100) to see how `\triCorr` and `\triRatio` move
requires re-running the real leg. It cannot be done from committed evidence, because
`results/triangulate_metrics.json` keeps only summary scalars, not the per-channel `miss`, `lon`,
`lat` and `r_geom` arrays `triangulate_track` returns --- the `innerrc` lesson (a results file
omitting the numbers its own headline is computed from). No STEREO/WAVES cache exists locally, so
the sweep means re-fetching the 2013-05-15 L3 DF data. **Commit the per-channel arrays on that run**
so the cut becomes auditable and any future sweep is offline.

### Resolved 2026-08-24 (real re-run): the sweep is measured, and the cut is not load-bearing

`miss_sweep` re-applies the threshold to a track built with the cut open (pure filtering --- the
per-channel misses were already computed), and the real 2013-05-15 leg was re-fetched from SPDF to
run it. Purely additive: every previously committed value is unchanged.

| cut (R_sun) | n | corr | ratio |
|---|---|---|---|
| 15 | 12 | 0.977 | 2.20 |
| 30 | 36 | 0.975 | 2.21 |
| 60 (analysis) | 38 | 0.989 | 2.18 |
| 100 | 38 | 0.989 | 2.18 |

**Neither headline depends on the cut.** Tightening 60 -> 15 keeps only 12 of 38 channels and moves
the correlation by 0.012 and the ratio by 0.02; no kept channel misses by more than 60, so the
100 row is identical. The cut is permissive, but it is not what makes the numbers --- now stated in
the paper with the sweep macro-backed (`\triSweep*`).

The per-channel arrays (freq, r_geom, r_plasma, miss, lon, lat, cut open) are committed as
`results/triangulate_channels.csv` (38 channels), so the 60 R_sun choice is auditable and any
future sweep is offline.

## Full referee round (2026-08-25): MAJOR REVISION, 15 findings, two BLOCKERs

The evidence discipline is better than average (all 38 channels committed, a real miss sweep,
correct merge wiring on JSON+macros) — but the two numbers the interpretation is built on are
respectively wrong and the wrong summary of the data. Both are recoverable from the committed
CSV alone, and the corrected version is a STRONGER paper.

**BLOCKER 1: `\triLon` = 168.9 is a scalar median taken across the ±180° wrap.** 16 of 38
committed longitudes are negative (−179.8…−169.2), 22 positive (159.5…+179.8); the scalar
median lands at the 7.9th percentile of its own unwrapped sample. Wrap-safe values: unwrapped
median 179.5°, circular mean 178.7° (R = 0.993) — the published longitude is ~10° wrong. The
module's own docstring states the principle ("scalar angle averaging is wrong near the azimuth
wrap") and then the output summary violates it. The fixture longitude (35°) is 145° from the
cut, so no test can see it.

**BLOCKER 2: the discrepancy is additive and nearly constant, not the multiplicative "factor
2.18".** From the committed CSV: r_geom − r_plasma = 14.2 ± 3.3 R⊙ (constant while r_plasma
varies 13×); OLS r_geom = 1.118·r_plasma + 12.3 (slope ≈ 1 → NO room for a density
enhancement, which multiplies); additive model rms 3.4 R⊙ vs 17.9 for "2.18×" — 5.3× better
with one parameter each. The ratio runs 1.29→3.65 monotonically with frequency (Spearman
0.90), so 2.18 (a median, never labelled one) describes no channel. The paper's own DF-bias
mechanism predicts the OPPOSITE trend (largest fractional bias at low frequency; measured:
smallest). The implied cause of 13–14 R⊙ at a 205 R⊙ lever arm is a ~3.7° constant DF error; a
(−2.3°, +2.8°) constant azimuthal bias pair reproduces the whole r_geom track to 0.055 dex
with ZERO density enhancement. The "denser-than-average corona" inference describes the
residual of a model the data reject (and even at face value needs density_scale ≈ 20,
unstated).

**MAJORs:** r = 0.989 cannot fail (any monotone power law scores 0.95–0.9999; bare 1/f
correlates with r_geom at 0.9894, i.e. BETTER than Leblanc; the log-log slopes are −0.60 vs
−0.93 — the data fail the "correct log-log curvature" property the statistic is advertised to
test, and the computed space is never stated: linear 0.989, log-log 0.975); no uncertainty on
r or the ratio (channel jackknife 0.010; block bootstrap on the ratio ±0.15; harmonic 1 →
3.89; a 2–3° DF bias → ~1.0 — the systematics dwarf the statistics); the harmonic × density
degeneracy footnoted, never gridded (h=1,scale=4 ≡ h=2,scale=1 verified numerically), though
the tool supports it; the A/B burst windows are never aligned in absolute time
(`times = ep − ep[0]` per file; a few-minute offset produces exactly the constant few-degree
bias of Blocker 2, the fixture gives both spacecraft identical time arrays, and no committed
epoch exists to check — a human must run the provided CDF-epoch check); the triangulation is
unauditable (pos_a/pos_b and per-channel ua/ub not committed — no reader can reproduce one
r_geom); the channels CSV and figure have NO clobber guard and the bare CLI default
(`python -m jansky_research.triangulate`) overwrites the real CSV with synthetic channels
while the JSON keeps saying "STEREO-A+B" — a marker that lies; the committed arXiv package
still says "two decades" (retracted in-paper: 1.20 decades) and leaks the raw key
"leblanc1998" into the abstract.

**MINOR/NIT:** the miss-based uncertainty proxy is never compared to the effect it bounds
(median miss/r_geom = 0.60; 31/38 channels individually consistent with ratio 1 — though all
38 offsets are one-signed, so the offset is real); the offline validation cannot fail
(fixture built AT the Leblanc radii → r→1 by construction; sep_deg=135 locked in by a test
the findings file itself says is not a control for the real 82°); four hand-typed calibration
numbers (r≈0.75, ≈1.1, 9°, ≈1.6, ~25°) under a "none typed by hand" header; figure caption
"track in shape" is what the panel disproves; krupar2012 authors 2/3 wrong (same wrong triple
as type3synthesis); the sweep quotes its endpoints and skips its largest excursion (30 R⊙
row), and the 60/100 rows are vacuous (max miss 31.65).

**Checked and clean:** all committed scalars reproduce from the CSV; krupar2014 companion fix
correct and its band matches \triFlo–\triFhi exactly; leverage check PASSES (drop the largest
point: 0.979; drop ten: 0.896 — not the swaves shape); five other citations exact vs Crossref.

**Status: fixes pending.** The single change: replace the multiplicative framing with the
additive one — r_geom − r_plasma = 14.2 ± 3.3 R⊙, slope 1.12, one ~3.7° instrumental offset —
which retires the density enhancement, demotes r = 0.989, and turns the result into "Leblanc
confirmed by pure two-spacecraft geometry to ~12% in radius from 0.07 to 0.5 AU."

**Status: RESOLVED (2026-08-25).** One real re-run; every referee-side number reproduced in the
committed pipeline, and the two blockers closed the way the referee predicted.

**Blocker 1 (wrap-broken longitude):** `circular_median_deg` (centre on the circular mean, take
the median of wrapped deviations) replaces the scalar median everywhere: the published 168.9°
becomes **179.5° ± 7.6** (the referee's unwrapped median was 179.5 exactly). A fixture ON the
branch cut at the real 82° separation now exists and fails under the old code; the latitude
gains its own scatter (4.7 ± 4.6).

**Blocker 2 (multiplicative → additive):** the paper's headline comparison is now
`additive_vs_multiplicative`, all committed: diff 13.4 ± 3.3 R⊙ (block-jackknife ±1.1 on the
median), OLS slope 1.119 / intercept 12.3, additive rms 3.37 vs multiplicative 17.83 (5.3×),
ratio range 1.24–3.66 with the median 2.18 demoted to "the contrast, since a single scale
factor describes no channel". The regression slope ≈ 1 retires the density-enhancement
narrative outright; the implied constant pointing bias is 4.0° at the 194 R⊙ median lever. The
result is reframed as a recover-a-known for the Leblanc shape: one constant instrumental
offset, then 3.4 R⊙ rms (~12%) over 0.07–0.5 AU.

**MAJORs:**
- r = 0.989 demoted with its nulls stated (any monotone power law > 0.95; 1/f ≈ Leblanc);
  the statistic that can fail — the log-log slope of r_geom on r_plasma, 0.653 vs 1.0 — is
  committed and quoted as what forced the additive framing. Block-jackknife errors on r
  (±0.034) and the ratio (±0.19) committed.
- The (harmonic × density) grid is committed (`harmonic_density_grid`): h=1 gives diff 17.8 /
  ratio 3.89; h=1,scale=4 ≡ h=2,scale=1 exactly (verified in a test); the harmonic is stated
  as a convention inside the grid.
- The hand-typed calibration numbers are a committed measurement (`noise_bias_calibration` at
  the real 82° baseline, source at the unmodified Leblanc radii): 9°/18°/25° scatter →
  additive offsets 3.1/10.5/16.9 R⊙ — the observed 13.4 R⊙ needs no new physics, and the
  calibration produces ADDITIVE offsets (ratio 1.11–1.68), confirming the framing.
- Time is now carried as seconds since the file date's UTC midnight on both spacecraft; the
  committed epochs show the two per-file origins differed by 25 s (harmless inside ±900 s) —
  the hazard is closed and measured, and a shifted-origin fixture documents the sensitivity.
- Auditability: pos_a/pos_b in the JSON; per-channel mean direction vectors and sample counts
  in the CSV — every r_geom is now reproducible offline.
- Clobber guard: a synthetic run refuses to overwrite the channels CSV / figure when the
  results JSON on disk is real (tested) — the marker can no longer lie.
- The stale arXiv package regenerated (1917/1920 chars); the \citeyearpar leak fixed by
  switching the abstract to \citet, which the assembler resolves ("Leblanc et al. (1998)").

**MINOR/NIT:** the miss-vs-effect comparison is stated (median miss ≈ the offset; channels
individually consistent, ensemble one-signed on all 38); the sweep quotes its largest
excursion (30 R⊙ row) and states that the 60/100 rows are vacuous (max miss 32 R⊙ — the
analysis cut never binds); krupar2012's author list fixed (Santolik/Cecconi restored); the
findings file's superseded claims marked in place. Headline values that were supposed to be
stable stayed stable (n 38, sep 82.1, miss 17.1, corr 0.989).
