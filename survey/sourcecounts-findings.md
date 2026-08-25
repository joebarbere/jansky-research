# Findings — 1.4 GHz Euclidean-normalised source counts from NVSS

`jansky_research.sourcecounts` builds the differential radio source count dN/dS from a public NVSS
region and compares it, Euclidean-normalised, to the canonical 1.4 GHz counts (the Hopkins et al. 2003
polynomial fit). It is the first slice to exercise the `jansky.sourcecounts` helpers, and a clean,
deterministic recover-a-known: NVSS *should* reproduce the well-established counts, and it does.

## Method

The differential count dN/dS is the number of sources per unit flux per unit solid angle; a static
Euclidean universe gives dN/dS ∝ S^−5/2, so the **Euclidean-normalised** count S^5/2 dN/dS divides out
that slope and any real structure stands out as a departure from a flat line. We fetch every NVSS
source in a cone, cut at a completeness limit, bin in log-flux, form dN/dS with Poisson errors
(`jansky.sourcecounts.differential_counts`), divide by the cone solid angle 2π(1−cos θ),
Euclidean-normalise (`euclidean_normalised_counts`), and compare bin-by-bin to the Hopkins 2003
reference. The `jansky` helpers do the counting; the new code is the NVSS fetch, the solid-angle
normalisation, and the published reference curve.

## Recover-a-known: an 8° NVSS cone at (180°, +30°), b ≈ +80°

| quantity | value |
|---|---|
| NVSS sources (> 3.5 mJy, in cone) | **7428** |
| solid angle | 0.0611 sr |
| sample flux range | 3.5 mJy – 2.98 Jy |
| comparison range (Hopkins valid) | 3.5 mJy – 1 Jy, 10 bins (N ≥ 5) |
| differential slope (3.5 mJy–1 Jy, log–log) | **−1.91** (Euclidean = −2.5) |
| median ratio to Hopkins 2003 | **1.021** |
| scatter about Hopkins | **0.061 dex** (~14%) |

The NVSS counts reproduce the canonical 1.4 GHz Euclidean-normalised counts across ~2.5 decades in flux
(3.5 mJy–1 Jy) at the **0.061 dex** level — the pipeline and NVSS agree with the published Hopkins fit
to ~14%, bin to bin. The differential slope **−1.91** is *flatter* than the Euclidean −2.5: a slope
shallower than −2.5 means the Euclidean-normalised count S^5/2 dN/dS **rises** with S, so the counts
climb from faint fluxes toward the ~1 Jy bright-end peak. This is the well-known **sub-Euclidean**
behaviour of the 1.4 GHz counts below ~1 Jy (first measured by Condon 1984), a real
cosmological-evolution signature — not a flat Euclidean line. The synthetic fixture (fluxes drawn from
the Hopkins differential count) round-trips to ratio ≈ 1.00 with 0.03 dex scatter, confirming the
binning/normalisation is unbiased.

## Honest assessment & caveats

- **A reproduction/method demo, not a new measurement.** The Hopkins 2003 polynomial is itself a fit
  to a deeper multi-survey compilation; agreement validates the *pipeline* and the *NVSS data*, it does
  not measure new counts. The contribution is a tested, reproducible Euclidean-count pipeline plus a
  recover-a-known.
- **Single field → cosmic variance.** One 0.06 sr cone; large-scale clustering makes the counts vary
  field-to-field, so the 0.073 dex scatter folds in cosmic variance on top of Poisson — a wider-area or
  multi-field run would tighten it.
- **Faint-end systematics near the NVSS limit.** NVSS is ~50% complete near ~2.5 mJy and its 45″ beam
  resolves some extended sources into multiple catalogue components, both of which bias the lowest
  flux bins; we cut at 3.5 mJy to stay above the worst of it, but the faintest bin still carries
  Eddington bias and the resolution caveat.
- **Bright end is Poisson-starved and outside the Hopkins fit.** The Hopkins 2003 polynomial is
  formally valid only to 1 Jy; above that the cone holds only a handful of sources per bin, so we
  restrict both the Hopkins ratio and the slope fit to bins with ≥5 sources **and below 1 Jy** — the
  brightest sample sources (to 2.98 Jy) are plotted but excluded from the comparison.
- **Off-plane field by design.** (180°, +30°) is at b ≈ +80°, away from the Galactic plane, so the
  sample is extragalactic; a low-latitude field would mix in Galactic sources and a different count.
- **Reproducible:** `python -m jansky_research.sourcecounts --ra 180 --dec 30 --radius 8` regenerates
  the metrics, the Euclidean-normalised count figure with the Hopkins reference, and the macros from
  the public VizieR NVSS catalogue.

## Full referee round (2026-08-25): MAJOR REVISION, 14 findings

No fabricated number and no macro hole: the referee re-fetched the cone from VizieR TAP and
reproduced every committed value to the last digit, verified the Hopkins coefficients against
the e-print itself (astro-ph/0211068 line 684), and confirmed the estimator is unbiased
(Hopkins-through-the-estimator: 0.994, 0.0022 dex). What needs revision is claim strength and
precision — and most fixes make the paper stronger.

**MAJORs:**
1. None of the per-bin numbers behind the headline are committed (nine scalars; the figure is
   the only record of the bin-by-bin comparison). compute_counts already returns them; run()
   drops them. The innerrc lesson.
2. The headline 1.021 is quoted to a precision it lacks: the binning choice alone spans
   0.983–1.022 (n_bins 8–20) and Poisson is ±0.024 — honest form: 1.02 ± 0.03, distance from
   unity unresolved.
3. The entire binning hangs on the single brightest source (bins = geomspace to s_max·1.001):
   dropping the two brightest — both excluded from the comparison anyway — swings the quoted
   scatter 0.037–0.078 dex. Fix: fixed 0.2-dex grid.
4. The scatter is attributed to the wrong terms: cosmic variance is coherent across flux bins
   (it moves the MEDIAN — where the measured 2.1% excess is ~1σ against the referee's computed
   1.7–2.1% clustering + 1.2% Poisson) while the bin-to-bin 0.061 dex is Poisson (0.038) ⊕ the
   Hopkins fit's own quoted 0.04-dex residual = 0.055. The per-bin errors are computed,
   plotted, and never used: χ²/dof = 2.80 (p ≈ 0.002) against ratio 1 errors-only, falling to
   ~1.2 with the reference residual folded in. The "expected scatter" is asserted, never
   computed.
5. "Complete to ~2.5 mJy" is measurably wrong in this footprint (2.1–2.5 mJy bin sits at
   0.59 of Hopkins; the findings doc itself says ~50%), and the 3.5 mJy cut does NOT clear
   the roll-off: a monotone 1.28→0.93 threshold artifact persists at 9–16% inside the first
   used bin, which carries 32% of the sample — hidden by the 0.244-dex bin averaging its two
   halves to 0.994. No cut sweep exists; the referee's sweep (3.5→20 mJy: 1.021→0.976) shows
   stability the paper is entitled to claim only after running it. Bonus clean result: surface
   density uniform to ±1.3% across four annuli.
6. Component vs source is a headline systematic, not a faint-end caveat: pair counts show a
   1.5–1.6× excess at 60–90″ (~2.4% of components are extra components of counted sources);
   friends-of-friends merging moves the median 1.021→0.986 (60″) →0.949 (100″) — larger than
   the departure from unity. And Hopkins is FIRST-anchored above 2.5 mJy (their §4), so the
   whole comparison range is NVSS(45″) components vs FIRST(5″) counts — a convention match
   the paper never states.
7. The offline recover-a-known cannot fail: generator and comparator share the same reference
   and the same solid-angle expression, so a 2×-wrong normalisation and even a FLAT Euclidean
   reference pass every test; the 2π(1−cosθ) normalisation — the one new piece of code — has
   no coverage (verified correct by hand, but untested). The "published anchors" test asserts
   the polynomial against its own outputs. Fix: closed-form power-law test (analytic k S^-2.5
   over a known area must return k).

**MINOR/NIT:** slope −1.91 has no error (±0.03) and no comparand (Hopkins through the same
bins gives −1.87; and the slope is interval-dependent: −1.91→−2.07 over cuts 3.5→20 mJy);
the figure draws the fit 0.35 dex past its 1 Jy validity limit and does not distinguish used
from excluded bins (the −2.5σ point at 0.417 Jy IS used); "first measured by Condon (1984)"
overstates (an evolution-modelling synthesis, not the first sub-Euclidean measurement); the
findings doc's 0.073-dex scatter is stale (0.061 everywhere else); the packaged arXiv
abstract lost its \citealt (the curve attributed to nobody) and carries "TODO pages"; the
~14% prose figure is hand-derived from \scScatter; the N≥5 cut never fires in this run.

**Status: fixes pending.** The single change: replace the asserted agreement sentence with
the computed budget — bin scatter 0.038 Poisson ⊕ 0.040 Hopkins residual = 0.055 vs 0.061
measured; cosmic variance moved to the median (2.1% ≈ 1σ); χ²/dof 2.80 → ~1.2 with the
reference error — which forces the per-bin table, the uncertainties, and the fixed grid in on
the way through.
