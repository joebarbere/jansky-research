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
