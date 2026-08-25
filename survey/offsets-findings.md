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

## Referee round on the style conversion (2026-08-23)

Verdict *minor revision*; all 15 macros re-verified against `results/offsets_metrics.json`. Every
substantive caveat (MOJAVE flux limit, 15 GHz vs S/X, "nested subsets, not independent
confirmations") came through byte-identical. Three fixes:

1. **The retitled section asserted its result nowhere.** `\section{The offset direction aligns with
   the jet}` became a noun phrase, which is the right register --- but nothing in the body carried
   the assertion the old title made, so "The alignment is **also** directional" back-referenced a
   header and "This reproduces \citet{kovalev2017}" had no antecedent. Added one declarative
   sentence after the KS result, at the strength the numbers support.
2. **A collapsed appositive attached the 24.0x excess to the wrong noun.** With an em-dash the
   appositive was set off from the whole comparison; with a comma it attaches to "the Rayleigh
   expectation of 1.11%", which is not what 24.0 is (it is the ratio 26.7/1.11). Reworded to name
   the referent.
3. **A concessive link was cut, so a scoped robustness claim became a general one.** "... is an
   upper bound (structure-only estimates are nearer ~9%) **--- but** the excess over Rayleigh is
   robust" became two sentences, letting the robustness carry onto the 26.7%. It is the *existence*
   of the excess that is robust, not its size; now stated that way.

Same appositive collapse in the abstract ("a well-established result" attaching to the noise model
rather than to the heavy tail) and a lost comma before "since" were also fixed.

## Full referee round (2026-08-25): MINOR REVISION, 16 findings

Every headline number reproduced EXACTLY from an independent re-fetch of both catalogues
(N=3502, 0.584 mas, 26.70%, 24.0×, 414, 24.26°, KS p to the last digit); the ICRF3 unit traps
verified against the VOTable header; both merge guards wired; the offline suite is a genuine
recover/reject-null pair; the reproduction framing and MOJAVE scoping are honest. The
revision is for four statements not backed by the evidence — one changes the abstract.

**MAJORs:**
1. **The Rayleigh null ignores the RA/Dec error correlation ICRF3 ships in the same table**
   (the un-requested `Corr` column: mean ρ = −0.069, 29% with |ρ| > 0.2). The correlated-
   normal tail is always heavier than χ²₂: marginalising over the empirical ρ distribution
   gives P(X>3) = 1.273%, not 1.11%, so the abstract's 24.0× is ≈ **21×** (Gaia's
   `ra_dec_corr` still un-folded). The bias is one-signed. Fix: fetch both correlations, form
   X as a Mahalanobis distance, re-run the real leg.
2. **The "directional" test uses a null the axis result already excluded**: 50% within 45° of
   downstream vs "random 25%" — but given axis alignment the sign-symmetric null is
   frac(axis<45°)/2 = 35.1%. The honest test still passes (among 291 axis-aligned: 208
   downstream vs 83 upstream, binomial p = 1.5e-13), but the quoted contrast conflates two
   results and frac_axis_lt45 is not committed.
3. **The error-model paragraph is misattributed and contradicted by the data**: the "~9%"
   comparand is not in petrovkovalev2017 (the right number is 6% at the 99% level from the
   companion letter MNRAS 467, L71 — not in refs.bib); the 26.7%-vs-6% factor is DR1→DR3
   precision, not error underestimation; and the noise-core Rayleigh scale fit gives
   f = 0.97–1.00 (KS p ≈ 0.95) — the formal errors are NOT detectably underestimated, so the
   paper hedges in the wrong direction and its data support a stronger claim than it makes.
4. **A robustness claim whose numbers exist only in this findings file**: the delPA<45°
   variant (n=318, 22.3°, p=1.1e-22 — verified exact) is in no committed artifact and run()
   throws `delpa` away. Same for the X>2 row's directional numbers.

**MINOR:** "robust to any plausible error inflation" is never computed (it is true —
inflation sweep 1.2×→19.3×, 2×→11.5×, 3×→7.7×: four lines retire it permanently); 26.7%
varies 0.23–0.38 across formal-error quartiles (a population-convolved-with-errors quantity
presented as a population one); no uncertainty on any headline (binomial 0.75pp; a
declination jackknife gives 1.4× that with a monotone N–S gradient 22.3×→27.3× — known
zonal/scanning-law systematics; both hemispheres individually significant, worth reporting);
no chance-coincidence estimate (≈9 false matches, 0.27%, all landing in the tail) and 25
matches with offsets >100 mas (max X = 652 — false matches or host centroids; dropping >50
mas: 26.70→25.87%); no per-source table committed while the abstract promises "a reproducible
matched catalogue" (claimed but not shipped); "disk–jet signature" cites Plavin+2019 while
the 83 upstream (29%) — the disk half — go unreported; the figure's "tracks Rayleigh at small
X" is not what the density-normalised panel shows (a factor ~2 low below X≈2 from the
truncated normalisation, though the claim is true conditionally).

**NIT:** lister2021 third author wrong (Kellermann omitted — the lawrance2024 pattern);
synthetic_field's sigma_mas parameter documented and unused; x_cut dropped from the metrics;
the figure has no clobber guard (an offline --out . run leaves a synthetic figure under a
real JSON); \offJetKsP prints 3e-22 where the JSON carries 3.2e-22.

**Status: fixes pending.** The single change: replace the error-model paragraph with the two
measurements it needs (core Rayleigh scale 1.00 ± 0.02; the inflation sweep), cite the right
Petrov & Kovalev letter with the right number, attribute the difference to DR3 precision —
the paragraph the entire 26.7% rests on becomes the paper's strongest.
