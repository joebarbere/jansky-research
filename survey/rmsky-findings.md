# Findings — the Galactic Faraday rotation sky (Taylor+2009 RM catalogue)

`jansky_research.rmsky` maps the Galactic **Faraday rotation sky** — the line-of-sight integral
$\mathrm{RM}=0.81\int n_e B_\parallel\,\mathrm{d}l$ — through the rotation measures of extragalactic
sources, reusing `jansky.polarization` for the underlying $\chi(\lambda^2)$ measurement. (Isolating
$B_\parallel$ would need an electron-density model, which we do not apply.) This is the tooling +
real-data + recover-a-known leg (the real fetch is a single reliable VizieR query, so all three are done
together). The synthetic offline fixture is a **round-trip code check** (it injects a $\csc|b|$ disk +
sign field that the analysis then recovers); the actual physics check is the real-data run below.

## Data

The **Taylor, Stil & Sunstrum (2009)** NVSS rotation-measure catalogue — RM for **37,543** polarised
extragalactic sources at 1.4 GHz — on VizieR (`J/ApJ/702/1230`, fully public, no auth). We fetch
`(RAJ2000, DEJ2000, RM)`, convert to Galactic $(l, b)$, and analyse the RM sky.

## Recover-a-known: the two textbook Galactic-RM signatures

**1. Plane enhancement.** Sightlines near the Galactic plane traverse more magneto-ionic disk
($\propto\csc|b|$), so $|\mathrm{RM}|$ rises sharply toward $b=0$:

| $\lvert b\rvert$ (deg) | median $\lvert\mathrm{RM}\rvert$ (rad m⁻²) | N |
|---|---|---|
| 0–10 | **62.1** | 3 821 |
| 10–30 | 34.1 | 11 868 |
| 30–60 | 16.5 | 15 587 |
| 60–90 | **11.5** | 6 267 |

The plane/pole **enhancement ratio is 5.4 ± 0.15** (bootstrap) — a clean, monotonic, highly significant
disk signature. It is *softer* than the bare $\csc|b|$ limit (~11 for these bin centres) for two real
reasons: a flat $\sim$11 rad m⁻² extragalactic+intrinsic floor that does not scale with $|b|$, and
Faraday **depolarisation** in the thick disk, which suppresses the detectable source count near the
plane (only 3 821 sources at $|b|<10°$) so the *detected* near-plane RMs under-represent the true
population. The ratio is therefore a lower bound on the path-length effect.

**2. Sign organisation.** The mean (not median — we want the net ordered-field sign) RM is sign-organised
across the plane:

| region | mean RM (rad m⁻²) | N |
|---|---|---|
| inner, north ($b>0$) | **+9.3 ± 0.6** | 9 585 |
| inner, south ($b<0$) | **−23.9 ± 0.9** | 5 248 |
| outer, north | +5.1 | 12 602 |
| outer, south | −6.4 | 10 108 |

The net sense is positive above the plane, negative below, at high significance — consistent with the
large-scale Galactic field reported by Taylor+2009 and Sun et al. (2008). **Caveats on the precise
numbers:** (i) the true structure is a *quadrupole* — the $l<90°$ and $l>270°$ halves of our "inner"
mask carry **opposite** sign at a given $b$, so this coarse mask conflates them and the recovered means
are partial cancellations, a net-sign indicator rather than a field measurement; (ii) the inner
north/south count imbalance (9 585 vs 5 248) is largely an NVSS **coverage artefact** ($\delta>-40°$
cuts the inner-south sky), so the larger $|{-}23.9|$ vs $|{+}9.3|$ should *not* be read as a physical
north/south amplitude asymmetry. Both signatures reproduce the established Galactic RM sky.

## Honest assessment & caveats

- **A reproduction, not a discovery.** The tool recovers two well-known Taylor+2009 results from the
  public catalogue; the contribution is a tested, reproducible pipeline, not new astrophysics.
- **Two-band $n\pi$ ambiguity.** Taylor+2009 RMs come from only the two NVSS IF bands, so individual
  large $|\mathrm{RM}|$ values can be aliased by $\pm652\ \mathrm{rad\,m^{-2}}$ (Brentjens & de Bruyn
  2005); the ~50 aliased sources are 0.13% of the catalogue, so **medians** and **means over thousands**
  are insensitive to them.
- **The polar floor is mostly intrinsic extragalactic RM.** The $\sim$11 rad m⁻² median at the poles is
  dominated by the sources' own (redshift-diluted) RM (~7 rad m⁻² rms; Schnitzeler 2010; Mao et al.
  2010), with a smaller high-latitude Galactic foreground and only a minor ($\sim$1–2 rad m⁻²)
  ionospheric residual. This floor adds scatter but does not bias the large-scale pattern the
  medians/means isolate.
- **Northern-sky catalogue.** NVSS covers $\delta>-40°$, so the southern Galactic sky is incomplete;
  the quadrant split is a coarse probe of the field geometry, not a full harmonic decomposition.
- **Reproducible:** `python -m jansky_research.rmsky` regenerates the metrics, the Aitoff RM-sky +
  $|b|$-profile figure, and the macros from the public VizieR catalogue.
