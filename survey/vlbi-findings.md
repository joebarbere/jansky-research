# Findings — multi-decade VLBI flux variability of compact AGN (Astrogeo)

`jansky_research.vlbi` builds per-source, multi-decade **VLBI light curves** from the Astrogeo database
(Petrov) and characterises each source's variability with the same transient-survey statistics we built
for VLASS (`vlass.variability_metrics`: η = weighted reduced χ² vs. constant; V = coefficient of
variation) plus a dual-band **S/X spectral index** (`spectra.spectral_index`). This is the real-data +
recover-a-known leg.

## Data access (verified, no auth)

Astrogeo exposes **per-epoch correlated-flux files** at
`http://astrogeo.org/images/{J2000}/{J2000}_{BAND}_{YYYY}_{MM}_{DD}_{analyst}_cfd.tab` — one tiny
two-line file per source per session per band. We read the integrated total flux density `Fl_int` (Jy)
and the image noise `Fl_noi`; the per-source directory listing enumerates the epochs. No astroquery/VO
interface exists — plain `requests` against the public HTTP tree. **Key assumption:** geodetic/astrometric
VLBI has no per-observation flux error, so we adopt a calibration-fraction floor `err =
sqrt((0.05·Fl_int)² + Fl_noi²)` (`VLBI_CAL_FRAC = 0.05`; Petrov & Kovalev 2025). Everything the
variability significance rests on hangs on this number — see the central caveat below.

## The validation set (a recover-a-known, not a blind survey)

18 well-known, well-observed compact AGN (`VALIDATION_SOURCES`): 14 Doppler-boosted blazars expected to
vary strongly, plus **4 compact symmetric objects (CSOs)** — OQ 208, 2021+614, 0108+388, NGC 3894 —
which lack a boosted core and are the **steady negative controls** (OQ 208 is the textbook stable VLBI
source). X-band light curves span 1995–2022 with 4–123 epochs each.

## Result 1 — the known variables are recovered

Every famous blazar comes back highly significantly variable, led by the recover-a-known anchor:

| source | name | X epochs | η | V | mean S$_X$ (Jy) | α$_{SX}$ |
|---|---|---|---|---|---|---|
| J0854+2006 | **OJ 287** | 123 | **208.6** | 0.48 | 3.79 | +0.42 |
| J0238+1636 | AO 0235+164 | 52 | 104.0 | 0.41 | 1.44 | +0.19 |
| J2253+1608 | 3C 454.3 | 9 | 93.7 | 0.35 | 13.69 | +0.09 |
| J2202+4216 | **BL Lac** | 58 | 47.5 | 0.54 | 3.23 | +0.19 |
| J1512-0905 | PKS 1510-089 | 30 | 58.7 | 0.39 | 2.65 | −0.16 |

OJ 287 and BL Lac — two of the most famous variable blazars in the sky — are recovered as the most
variable, with flat/inverted S/X spectra (α$_{SX}\simeq0$ to +0.4) as expected for self-absorbed
compact cores. The population median α$_{SX}$ = +0.03 (flat), the compact-core signature.

## Result 2 — the steady controls expose the central caveat

The 4 CSO controls are **also** formally "variable" by the absolute χ² test — OQ 208 has η = 79.6,
p ≈ 0. That is **not** intrinsic variability: a CSO is steady, so its apparent scatter is amplitude
calibration + (u,v)-coverage / resolved-structure differences between sessions. This is the smoking gun
that **the 5% error floor is too optimistic** and the absolute η/χ² massively over-rejects — *every*
source, steady or not, is "significantly variable" against a 5% error.

**What survives is the amplitude V, benchmarked against the controls.** The CSOs set an empirical
variability floor **V = 0.193** (their median); above it sit **13 of the 14** non-control AGN. The
boosted blazars have median V = 0.32 vs the controls' 0.19 — **~1.7× more variable in amplitude** — the
genuine, calibration-robust signal. The one non-control below the floor is **3C 273** (V = 0.13): bright
but resolved, with only 4 epochs, so its core variability is diluted in the total flux — an honest edge
case, not a failure.

## Honest assessment

- **The method works; the absolute significance does not.** The recover-a-known is clean (known
  blazars rank top; CSOs are the least variable), but η/χ² is uninformative here because the per-session
  flux error is not 5% — it is ~19% (the CSO floor). **V relative to a steady-control floor is the
  trustworthy discriminant**, not absolute η. The tool reports both and the figure draws the floor.
- **VLBI total flux ≠ intrinsic flux.** `Fl_int` folds (u,v) coverage and resolved structure; the CSO
  controls quantify how large that structural floor is (V ≈ 0.19). Any "variability" below it is noise.
- **A curated validation set, not a survey.** 18 hand-picked sources demonstrate the pipeline and the
  recover-a-known; a blind catalogue would need a large, magnitude-limited input list and the
  control-floor calibration applied per declination/epoch-count bin.
- **Reproducible:** `python -m jansky_research.vlbi --online` regenerates `results/vlbi_candidates.csv`,
  the metrics JSON, and the η–V figure from the public Astrogeo files.
