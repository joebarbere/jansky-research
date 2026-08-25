# Findings — a type III electron beam from the corona to 0.4 AU, geometrically validated

`jansky_research.type3synthesis` unifies four separately-validated slices into one reproducible
**drift-to-distance** framework for solar type III bursts, and adds the element no single slice provides:
an **independent geometric check** on the density-model distance the whole method rests on.

## The unified ladder (each leg a recover-a-known on public data)

| regime | slice | instrument | density model | reach (recover-a-known) | speed |
|---|---|---|---|---|---|
| corona | `solarbursts` | e-Callisto (ground) | Newkirk | ~1.5–3 R⊙ | 0.14 c |
| inner heliosphere | `windwaves` | Wind/WAVES RAD2 | Leblanc | 2.4 → 10.2 R⊙ (Alfvén surface) | 0.083 c |
| interplanetary | `swaves` | STEREO/WAVES HFR | Leblanc | 2.3 → 82.6 R⊙ = **0.38 AU** | 0.150 c |
| geometric | `triangulate` | STEREO-A+B direction-finding | **none** | 15 → 106 R⊙ | — |

The same drift-to-distance method, with a Newkirk (corona) → Leblanc (heliosphere) density-model handoff,
spans nearly three decades in frequency (78.94 MHz -> 0.125 MHz, i.e. 2.80; the "~100 MHz"
here is where an earlier "more than three decades" came from -- the committed metric is 78.94) and tracks the beam from the low corona
to ~0.4 AU.

## The centerpiece: a model-free check on the model distances

Every drift-to-distance result is only as good as the assumed `n_e(r)`, and that assumption is almost
never checked independently. Here it is: **`swaves` and `triangulate` analyse the same 2013-05-15 event**,
so the STEREO/WAVES plasma-frequency distance `r_plasma(f)` and the STEREO-A+B geometric (triangulated)
distance `r_geom` are two *independent* estimates of the same quantity. Across the triangulated channels
they correlate at **r = 0.989**. Both estimators decrease with frequency by construction, so *some*
correlation is trivial (a linear ramp already gives r ≈ 0.75 vs Leblanc); what the measured value adds
is agreement in the **log-log curvature** of the density-model track over the triangulated range
(15.3-106.1 Rsun, i.e. 0.84 decades, not the "two decades" claimed earlier), not just a shared
monotone trend. The absolute geometric scale runs ~2× high (ratio ≈ 2.18): direction-finding outward
bias on a degrees-wide source, a likely active-region density enhancement, **and radio-wave scattering
shifting the apparent interplanetary source outward of the plasma level** (Krupar et al. 2015) all push
the same way. The robust statement is that the two independent estimators **track each other in shape**,
not that they agree in absolute calibration.

## Honest assessment & caveats

- **A synthesis/reproduction, not a discovery.** Each leg is an existing recover-a-known; the new element
  is the unified framework + the geometric validation of the model distance.
- **Different events.** The corona (2011-09-14), inner-heliosphere (2003-10-28 X17), and interplanetary
  (2013-05-15) legs are *different bursts*, so this is one method validated across four regimes, **not one
  beam tracked end to end**. The stronger demonstration — a single event caught simultaneously by
  e-Callisto + Wind + STEREO — is the natural next step (an optional GATE-0 event hunt).
- **The cross-check is event-specific and shape-only.** It validates the 2013-05-15 interplanetary
  distances; the absolute geometric scale is direction-finding–biased, so it checks the distance *shape*,
  not the calibration.
- **Inherited per-leg systematics.** Harmonic vs fundamental (factor 2 in distance), peak-time vs onset
  speed bias, the average density model, and the coarse-cadence fit statistics are all bounded in the
  component slices and carried through here.
- **Two stitched density models, not one.** The Newkirk (corona) and Leblanc (heliosphere) models are
  not designed to join continuously; at the ~15–20 MHz handoff they differ in heliocentric radius by
  ~50% (white-light active corona vs Wind type III drifts). The ladder is two stitched regimes.
- **Not the first geometric check.** Prior STEREO direction-finding work (e.g. Krupar et al.) has
  compared triangulated positions to density-model predictions; the new element here is the *same-event,
  reproducible, open-pipeline* comparison, not direction-finding as a density-model check per se.
- **Reproducible:** `python -m jansky_research.type3synthesis` regenerates the unified ladder + the
  2013-05-15 cross-check figure and the macros; `make reproduce` runs the four recover-a-known events on
  the public archives (e-Callisto, SPDF Wind/WAVES + STEREO/WAVES + STEREO L3 DF).

## Full referee round (2026-08-25): MAJOR REVISION, 15 findings, two BLOCKERs

**BLOCKER 1: five committed per-leg values are pre-revision sibling vintages.** The metrics
were last written in #235; the siblings were revised in #263 (solarbursts) and #266
(windwaves+swaves). Verified exact by `git show`: corona speed 0.1347 / r 1.622–2.572 are the
pre-#263 values (current: 0.1173, 1.757–2.285); helio speed 0.0831 is windwaves' pre-#266
flat-OLS — which now survives in windwaves' own JSON only as `speed_c_points`, the REJECTED
estimator, i.e. the synthesis quotes a diagnostic as its sibling's headline; ip speed 0.1503
is swaves' pre-#266 non-converged fit (current 0.1604 ± 0.014). The geometric leg is current.
No guard can catch a stale-but-real file: preserve_live_results arbitrates real-vs-synthetic,
not cross-slice vintage. Fix is one real re-run (which re-runs all four siblings with
--out . — plus a durable test asserting each per-leg value equals the sibling's committed key).

**BLOCKER 2: the one new claim is contradicted by its own committed channels.** "Track each
other in shape … with a factor 2.18 offset in absolute scale" — from
results/triangulate_channels.csv the ratio r_geom/r_plasma rises monotonically 1.29→3.65 with
frequency (rank corr 0.90); the log-log slopes are −0.60 vs −0.93 (regression slope 0.653,
not 1); a constant ADDITIVE offset of 13.4 ± 3.3 R⊙ (≈ the median ray miss) beats the
multiplicative model 0.120 vs 0.279 fractional rms. The discrepancy IS the shape; the scale
offset is the thing that is approximately absent. "Validating by geometry the shape of the
Leblanc distances" asserts the inverse of the data.

**MAJORs:** r = 0.989 is linear-space Pearson (the "log-log curvature" sentence describes a
statistic that was not computed; log-log 0.975; any bare power law scores 0.954–0.9999, and
1/f beats Leblanc), quoted with no error though a channel jackknife is trivial; the
cross-check never uses the STEREO/WAVES drift leg at all (r_plasma is the Leblanc model
evaluated at the triangulated frequencies — a deterministic function of frequency, not a
second measurement; "two independent measurements … same burst" oversells, and the r/2.18
numbers are already papers/triangulate's headline, so the synthesis must state what it adds);
bare per-leg speeds/reaches where the source slices now refuse to quote them bare (windwaves
declines the Alfvén-region claim the synthesis makes; the title's "0.4 AU" is the HFR band
edge through the model — confirmed to four figures — and moves 0.19–0.77 AU across swaves'
committed grid); the headline ladder reach 0.493 AU is the geometric estimator the paper
itself declares biased outward (use the plasma 0.384 AU and report 106 R⊙ separately); title
and abstract assert one beam that the discussion correctly retracts (three events, twelve
years apart); after the re-run the corona figure segment will be internally inconsistent
(solarbursts' radii now come from the FITTED band 32.4–62.4 MHz, not the full detection band —
the "figure draws the fit it captions" fix never propagated here); the committed figure
predates every sibling revision AND its own metrics file (#156 vs #235); krupar2015 does not
support the scattering claim it is load-bearing for (its abstract is beam deceleration — and
incidentally DOES support harmonic=2, uncited).

**MINOR/NIT:** krupar2015 + krupar2012 both carry the same wrong author triple
(Crossref-verified); the ~50% Newkirk/Leblanc handoff discontinuity is hand-typed under a
"none typed by hand" header (computed: 41–52% at 15–20 MHz, and the legs' actual overlap is
10–13.8 MHz where it is 36–45%); the findings doc still carries every stale number and says
"four regimes" where the paper says three; the offline tests cannot fail (corr > 0.8 clears
for any monotone; a log-log-slope test would fail on the real data today, which is the
point); newkirk1961 page range possibly constructed; "four public instruments" counts the
STEREO/WAVES suite twice; 2.18 is a median and never labelled one.

**Status: fixes pending.** The single change: replace "shape agreement with a 2.18× scale
offset" with what the channels show — a ~13 R⊙ additive offset at the scale of the ray miss,
fractional discrepancy growing 30%→3× across the band, NOT a validation of the Leblanc shape
— which forces "geometrically validated" out of the title and makes the bias paragraph a
measurement instead of an excuse.
