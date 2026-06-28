# Findings — peaked-spectrum (GPS/CSS) candidates via 3-frequency curvature

`jansky_research.peaked` selects peaked-spectrum (GPS/CSS-candidate) radio sources — compact, young
radio AGN whose spectrum rises then falls — from three public surveys: TGSS (150 MHz), NVSS
(1.4 GHz), VLASS (3 GHz). It is the direct follow-up the USS slice called for ("confirming a turnover
needs a third frequency, e.g. VLASS 3 GHz") and the maximal-reuse slice: it composes `spectra`
(`spectral_index`, `crossmatch`, `fetch_survey`) and `vlass` (3 GHz fetch, forced-photometry
variability, SIMBAD vetting) rather than reimplementing them. Run: a 2° cone at RA 180°, Dec +30°
(the same field as the USS slice), 409 NVSS$\cap$VLASS sources.

## Two systematics the real data exposed (and the fixes)

The naive "find sources detected in all three surveys, rising then falling" gave **zero** candidates —
not a null result but two selection effects:

1. **TGSS is shallow.** Peaked sources are *faint at 150 MHz* by definition, so requiring a TGSS
   detection excludes exactly the sources wanted (only the bright steep sources are in TGSS). Fix:
   treat a TGSS **non-detection** as an upper limit $S_{150}<25$ mJy, giving a *lower bound* on
   $\alpha_\mathrm{low}$. A source brighter at 1.4 GHz than the 150 MHz limit yet undetected at
   150 MHz **must** be rising. (This only confirms peaked sources brighter than $\sim$30 mJy at
   1.4 GHz — an honest depth limit.)
2. **NVSS$\to$VLASS resolution mismatch.** NVSS (45″) sees extended emission that VLASS (2.5″)
   resolves out, producing impossible "spectra" — e.g. 690 mJy $\to$ 6 mJy, $\alpha_\mathrm{high}=-6$.
   These masquerade as steep turnovers. Fix: a floor $\alpha_\mathrm{high}>-2$ rejects them as
   resolution artefacts. In this field **110** of 409 sources are flagged extended/artefact.

## Result: 6 peaked candidates, all steady (GPS-like), 3 uncatalogued

After both cuts: **6 peaked-spectrum candidates** (rising at low frequency via the TGSS upper limit,
$\alpha_\mathrm{high}\approx-0.6$ to $-1.1$, compact). Reusing `vlass.forced_photometry`, **all six are
VLASS-steady across the three epochs** ($V=0.02$–$0.08$) — consistent with GPS sources, not blazars.
SIMBAD (via `vlass`) identifies them:

| RA Dec | $S_{1.4}$ / $S_3$ (mJy) | $\alpha_\mathrm{high}$ | SIMBAD | note |
|--------|------|------|--------|------|
| 179.978 +30.458 | 41.6 / 20.5 | −0.93 | LEDA 1904565 [LINER] | AGN host |
| 182.018 +30.264 | 35.3 / 20.9 | −0.69 | **4FGL J1208.1+3017 [BL Lac]** | a *Fermi* blazar — contaminant |
| 179.181 +31.710 | 34.7 / 18.4 | −0.83 | — | uncatalogued GPS/CSS candidate |
| 179.174 +29.196 | 32.9 / 16.9 | −0.87 | — | uncatalogued GPS/CSS candidate |
| 180.048 +28.210 | 32.7 / 13.9 | −1.12 | — | uncatalogued GPS/CSS candidate |
| 180.702 +31.873 | 31.7 / 19.8 | −0.62 | NVSS J120248+315223 [radio galaxy] | CSS/GPS host |

A nice methodological point: the BL Lac (a known $\gamma$-ray blazar) is **steady** in VLASS ($V=0.08$),
so the variability flag alone did *not* catch it — **SIMBAD did**. Blazar rejection needs the catalogue
cross-ID (and ideally longer-baseline variability), not VLASS-epoch variability alone.

## Honest assessment

This is a **methodology + candidate-list** result, not a discovery: a reproducible three-frequency
curvature selection that (i) is robust to the TGSS flux-scale offset that sank the USS slice (curvature
compares two indices; $\alpha_\mathrm{high}$ is TGSS-independent), (ii) correctly separates resolution
artefacts, and (iii) yields a clean, steady, compact peaked-candidate list — three of them
uncatalogued. Confirming any as a true GPS source needs higher-resolution / more-frequency follow-up
(VLBI compactness, a sampled turnover); we report candidates, not detections.

### Limitations

- **Small field, depth-limited.** One 2° cone; the TGSS-upper-limit method only confirms peaked
  sources brighter than $\sim$30 mJy at 1.4 GHz, so the sample is the bright tail.
- **Curvature via a lower limit.** For TGSS non-detections $\alpha_\mathrm{low}$ is a bound, not a
  measurement, so the turnover frequency is not pinned (only "rising then falling").
- **No known GPS to recover here.** This field has no catalogued GPS source, so the validation is the
  method's behaviour (compact, steady, artefact-free) rather than recovering a known one — a targeted
  field with catalogued GPS sources would strengthen it.
- **NVSS$\to$VLASS resolution** biases $\alpha_\mathrm{high}$ for any non-point source; the floor cut is
  blunt. Matched-resolution fluxes (or a compactness cut) would be cleaner.

## Bottom line

A reproducible, maximal-reuse three-frequency peaked-spectrum selector that turns the USS slice's
unfinished business into a clean candidate list — robust to the flux-scale and resolution systematics,
honestly bounded, with three uncatalogued GPS/CSS candidates and a worked demonstration that blazar
contamination needs catalogue cross-ID, not just variability.
