# Findings — radio–optical position offsets of AGN (ICRF3 × Gaia DR3)

`jansky_research.offsets` reproduces the well-established result that the VLBI radio and Gaia optical
positions of AGN are systematically displaced: the *normalised* offset
$X=\sqrt{(\Delta\alpha^*/\sigma_\alpha)^2+(\Delta\delta/\sigma_\delta)^2}$ has a heavy tail far beyond
the Rayleigh expectation for pure Gaussian astrometric noise (Mignard et al. 2016; Petrov & Kovalev
2017; Kovalev et al. 2017; Lindegren et al. 2018; Plavin et al. 2019). It is catalogue-only and
maximal-reuse (cross-match + the project conventions), with no blocked archives.

## Data path (and the traps)

- **ICRF3** is VizieR **`J/A+A/644/A159`** (Charlot et al. 2020), not `I/367`; the S/X catalogue
  (4536 sources) is `table10`. Two unit traps: `RAICRS`/`DEICRS` are **sexagesimal** (RA in *hours*) —
  parsed with `SkyCoord(unit=(hourangle, deg))` — and `e_RAICRS` is in **time-seconds**
  ($\sigma_{\alpha^*}=e_\mathrm{RA}\times15000\cos\delta$ mas) while `e_DEICRS` is in **arcsec**
  ($\times1000$ mas).
- **Gaia DR3** match via the CDS **X-Match** service (`vizier:I/355/gaiadr3`), nearest within 0.5″;
  Gaia `e_RAdeg`/`e_DEdeg` are already mas on $\alpha\cos\delta$ and $\delta$.

## Result: the offset excess is reproduced (24×)

Over **3502** ICRF3∩Gaia AGN:

| quantity | value |
|---|---|
| median raw radio–optical offset | **0.58 mas** |
| fraction with $X>3$ | **26.7%** |
| Rayleigh expectation (pure Gaussian noise) | 1.11% |
| **excess** | **24×** |

The $X$ distribution tracks the Rayleigh curve at small $X$ and departs sharply in the tail — exactly
the structural-offset signature of Mignard et al. (2016) and Lindegren et al. (2018): a large AGN
population whose optical photocentre is displaced from the VLBI core by milliarcsecond-scale structure.

## Honest assessment & caveats

- **Reproduction, not discovery.** This recovers a known, much-studied result with a small tested
  tool — a reproducibility/tooling contribution (a validated offset + normalised-offset statistic and
  a reproducible matched catalogue), on-brand for this repo.
- **The exact $X>3$ fraction depends on the error model.** Gaia and VLBI *formal* errors are mildly
  underestimated (Gaia DR3 by $\sim$5–30% depending on magnitude/colour), which inflates $X$, so the
  26.7% is an **upper bound** on the structurally-offset fraction (the literature's structure-only
  estimates are nearer $\sim$9%; Petrov & Kovalev 2017). The **excess over Rayleigh is robust**
  regardless — no error inflation turns a 24× excess into none.
- **Catalogue-only scope.** The established physical interpretation — offsets aligning with the
  parsec-scale VLBI jet (Kovalev et al. 2017) and the disk–jet decomposition (Plavin et al. 2019) —
  needs VLBI jet position angles beyond this match; the module records the offset PA for that future
  step. A magnitude-resolved error-inflation model would tighten the significant fraction.
