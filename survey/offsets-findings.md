# Findings — radio–optical position offsets of AGN (ICRF3 × Gaia DR3)

`jansky_research.offsets` reproduces the well-established result that the VLBI radio and Gaia optical
positions of AGN are systematically displaced: the *normalised* offset
$X=\sqrt{(\Delta\alpha^*/\sigma_\alpha)^2+(\Delta\delta/\sigma_\delta)^2}$ has a heavy tail far beyond
the Rayleigh expectation for pure Gaussian astrometric noise (Mignard et al. 2016; Petrov & Kovalev
2017; Kovalev et al. 2017; Lindegren et al. 2018; Plavin et al. 2019) — and reproduces the deeper
result that the offset *direction* aligns with the parsec-scale jet (adding MOJAVE jet PAs). It is
catalogue-only and maximal-reuse (cross-match + the project conventions), with no blocked archives.

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

## Result: the offset direction aligns with the parsec-scale jet (Kovalev+2017 / Plavin+2019)

The magnitude excess shows the offsets are real; their **direction** shows what they trace. Adding
**MOJAVE XVIII** (Lister et al. 2021; VizieR `J/ApJ/923/30/mojave18`, "mean innermost jet position
angle"), positionally matched within 1″ to the ICRF3×Gaia AGN → **414 sources** with both an offset and
a jet PA. Jet-axis angle = `min(Δ, 180−Δ)` between offset PA and jet PA (0° = along the jet axis;
random median 45°, random frac<30° = 1/3):

| sample | n | median jet-axis angle | frac within 30° of jet | frac downstream (<45°) | KS p vs uniform |
|---|---|---|---|---|---|
| **all matched (primary)** | 414 | **24.3°** | **0.57** (rand 0.33) | **0.50** (rand 0.25) | **3.2×10⁻²²** |
| X>2 (consistency check) | 252 | 18.9° | 0.66 | 0.57 | 4×10⁻²⁵ |
| delPA<45° (sensitivity) | 318 | 22.3° | 0.60 | — | 1.1×10⁻²² |

The offsets point **along the jet axis, predominantly downstream** — exactly Kovalev, Petrov & Plavin
(2017) and Plavin, Kovalev & Petrov (2019): the VLBI core sits upstream (opacity/core-shift), the
optical centroid downstream along the extended optical jet. The alignment **tightens with offset
significance** (weak offsets are noise with random PA) and **survives the jet-wobble cut** (delPA<45°),
so it is not a tunable-cut artefact. Three independent datasets (ICRF3 / Gaia / MOJAVE) — no circularity.

## Honest assessment & caveats

- **Reproduction, not discovery.** This recovers a known, much-studied result with a small tested
  tool — a reproducibility/tooling contribution (a validated offset + normalised-offset statistic and
  a reproducible matched catalogue), on-brand for this repo.
- **The exact $X>3$ fraction depends on the error model.** Gaia and VLBI *formal* errors are mildly
  underestimated (Gaia DR3 by $\sim$5–30% depending on magnitude/colour), which inflates $X$, so the
  26.7% is an **upper bound** on the structurally-offset fraction (the literature's structure-only
  estimates are nearer $\sim$9%; Petrov & Kovalev 2017). The **excess over Rayleigh is robust**
  regardless — no error inflation turns a 24× excess into none.
- **MOJAVE selection bounds the alignment claim.** MOJAVE is flux-limited (>1.5 Jy at 15 GHz) — strongly
  beamed blazars with well-defined jets near the line of sight, the population where the alignment is
  strongest. The result applies to that population, not AGN in general (the source papers are likewise
  scoped). The jet PA is 15 GHz; the radio position is the ICRF3 S/X core — standard to combine, and the
  inner jet is essentially straight on parsec scales.
- **Two reproductions, framed as such.** Both the excess (Mignard/Lindegren) and the jet alignment
  (Kovalev/Plavin) are recovered, not discovered; the value added is a fully open, end-to-end
  ICRF3×Gaia×MOJAVE pipeline. The full 414-source sample is the primary test; the X-cut rows are a
  qualitative "strengthens with significance" check (nested subsets), not independent confirmations.
- A magnitude-resolved error-inflation model would tighten the significant-offset fraction (a future
  refinement on the magnitude side).
