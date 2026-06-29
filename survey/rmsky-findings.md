# Findings — the Galactic Faraday rotation sky (Taylor+2009 RM catalogue)

`jansky_research.rmsky` maps the large-scale Galactic magnetic field through the rotation measures of
extragalactic sources, reusing `jansky.polarization` for the underlying $\chi(\lambda^2)$ measurement.
This is the tooling + real-data + recover-a-known leg (the real fetch is a single reliable VizieR query,
so all three are done together).

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

The plane/pole **enhancement ratio is 5.4** — a clean, monotonic factor-of-five disk signature.

**2. Sign antisymmetry.** The hallmark of the organised Galactic field is a sign flip of the *mean* RM
across the plane in the inner Galaxy:

| quadrant | mean RM (rad m⁻²) | N |
|---|---|---|
| inner, north ($b>0$) | **+9.3** | 9 585 |
| inner, south ($b<0$) | **−23.9** | 5 248 |
| outer, north | +5.1 | 12 602 |
| outer, south | −6.4 | 10 108 |

The inner Galaxy is **positive above the plane, negative below** (and the outer Galaxy shows the same,
weaker, sense) — the antisymmetric disk/halo field of the Milky Way (Taylor+2009; Sun et al. 2008;
Mao et al. 2010). Both signatures reproduce the established Galactic RM sky.

## Honest assessment & caveats

- **A reproduction, not a discovery.** The tool recovers two well-known Taylor+2009 results from the
  public catalogue; the contribution is a tested, reproducible pipeline, not new astrophysics.
- **Two-band $n\pi$ ambiguity.** Taylor+2009 RMs come from only the two NVSS IF bands, so individual
  large $|\mathrm{RM}|$ values can be aliased by $\pm652\ \mathrm{rad\,m^{-2}}$; we use **medians** and
  **quadrant means**, which are robust to the rare aliased outliers, rather than individual extremes.
- **Intrinsic + extragalactic RM.** Each measured RM includes the source's own (redshift-diluted)
  contribution and the ionosphere; these add scatter (the $\sim$10 rad m⁻² floor seen at the poles) but
  do not bias the large-scale Galactic pattern, which is what the medians/means isolate.
- **Northern-sky catalogue.** NVSS covers $\delta>-40°$, so the southern Galactic sky is incomplete;
  the quadrant split is a coarse probe of the field geometry, not a full harmonic decomposition.
- **Reproducible:** `python -m jansky_research.rmsky` regenerates the metrics, the Aitoff RM-sky +
  $|b|$-profile figure, and the macros from the public VizieR catalogue.
