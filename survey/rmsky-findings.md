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

## Full referee round (2026-08-25): MAJOR REVISION, 15 findings

The traceability is intact (every macro resolves, the figure matches the committed profile
point for point, the quadrant-conflation caveat is at the right strength) — but the paper's
stated contribution is "honest uncertainties", and every error bar in it is an i.i.d.-source
estimator applied to a manifestly correlated field. `rmsky` is where `_ratio_bootstrap_se` is
*defined*; `rmstructure` imports that exact function and has already measured the same
statistic's i.i.d. bootstrap understating the block jackknife by 11× (±0.10 vs ±1.1). The
definer never got the fix the importer got.

**MAJORs:**
1. The headline **5.4 ± 0.15** is the i.i.d. bootstrap. If the SPICE-RACS inflation transfers
   even partially, the honest figure is nearer ±1 and "highly significant" changes strength.
   `rmstructure.spatial_block_jackknife` is public and takes exactly this interface.
2. "Each many standard errors from zero" (region means at 15.5σ/26.6σ) assumes 10⁴ independent
   sources where the exchangeable unit is a ~10² sky patch. A 5× inflation leaves 3σ; 11×
   leaves 1.4σ — the paper cannot currently tell those worlds apart, and the verb must follow
   the jackknife.
3. The nπ dismissal is argued for medians and applied to a mean: k aliased sources shift a
   region mean by k×652.9/n, so the paper's own k≈50 is 6.9 quoted SEs on inner-south — and
   aliasing is not sign-random and concentrates at low |b| where the sign signal lives. The
   sibling paper (`rmstructure`) treats the same effect as first-order. The "~50/0.13%" figure
   is uncommitted and uncited; ±652.9 is a Taylor+2009 number, not Brentjens & de Bruyn. Fix:
   |RM| < 300 variant + an alias-immune sign-fraction statistic.
4. A recover-a-known with no known: "match the literature" has no comparand anywhere — no
   published plane/pole ratio or quadrant mean is quoted. Quote Taylor+2009/Schnitzeler or
   downgrade to "consistent in sign and order of magnitude".
5. No per-source evidence committed: the whole result is 15 scalars + 4 bins; the catalogue is
   <1 MB gzipped. Commit `data/rmsky_taylor2009.csv.gz` (l, b, RM, e_RM) + fetch metadata.
6. The synthetic fixture is deterministic-signal + i.i.d. noise, so the i.i.d. bootstrap is
   *correct on the fixture by construction* — no offline test can expose finding 1. Import
   `rmstructure.synthetic_rm_screen` (correlated field) and assert jackknife > bootstrap.
7. The Aitoff longitude tick labels are −l (the map negates l and never relabels): the region
   printed 120° is l = 240°, and the paper's second result is a longitude-quadrant claim.

**MINOR/NIT:** outer-region SEs computed then discarded at the write step (a significance
claim with no error attached — the `innerrc` shape) and the outer counts exist only in this
file; the latitude profile is plotted with no errors; zero robustness variants and `e_RM` is
never even fetched; deduplication asserted by omission (self-match, one line); the λ²-fit
framing overstates `rm_from_angles`, which `run()` never calls; `\rm*` macros mode-dependent
and unnamespaced (CLAUDE.md's own recorded clobber, `\rmRatio` 5.4→8.4); `\rmTruth` ships as
`--` in the real macro file; 5.4 ± 0.15 mixes display precisions.

**Checked and clean:** all six DOIs vs Crossref; macro traceability; the internal arithmetic
(62.1/11.5 = 5.400, counts sum, csc|b| ≈ 11.08); the quadrant-conflation caveat.

**Status: fixes pending.** The single change that most improves the paper: block-jackknife
every quoted uncertainty (ratio + four region means) and keep the bootstrap only as the
contrast — whatever the jackknife returns, the paper stops being unable to say which world
it is in.

**Status: RESOLVED (2026-08-25).** One real re-run; the jackknife settled which world the paper
was in, and both signatures survive with honest errors.

1. **Every quoted uncertainty is now a 10° sky-block jackknife** (571 blocks;
   `rmstructure.spatial_block_jackknife`, the fix that had reached the importer but never the
   definer). Measured: the headline ratio is **5.4 ± 0.6** against the bootstrap's ±0.15 — a
   factor-3.9 understatement (not the SPICE-RACS 11×, and the paper quotes the gap as itself a
   measurement of the sky's correlation). The i.i.d. bootstrap is retained only as the contrast.
2. **"Many standard errors" is dead; the verbs now follow the jackknife.** Inner-north
   9.3 ± 3.3 (2.8σ), inner-south −23.9 ± 4.8 (5.0σ) — "significantly non-zero and of opposite
   sign ... at a few standard errors, not the tens the per-source errors would claim." The
   outer regions ship their SEs and counts (previously computed and discarded): outer-north
   5.1 ± 2.3 (2.2σ), outer-south −6.4 ± 4.4 (1.5σ, stated as not individually significant).
3. **The nπ exposure is bounded, not argued away**: the ±652.9 is attributed to Taylor+2009
   (not Brentjens & de Bruyn), the uncommitted "~50/0.13%" claim is gone, and the mean's alias
   exposure is handled by two committed instruments — the |RM|<300 variant (means move ≤0.2)
   and the alias-immune sign fractions (inner-north 0.589 ± 0.028, inner-south 0.348 ± 0.031;
   3.2σ/4.9σ from 0.5) — both reproducing the pattern.
4. **"Match the literature" downgraded** to "consistent in sign and order of magnitude with the
   published maps — no published value of this exact binned ratio exists to match more
   precisely." The csc|b| ~11 comparison stays as the only quantitative external anchor.
5. **The catalogue is committed**: `data/rmsky_taylor2009.csv.gz` (l, b, RM, e_RM; 470 kB) with
   `vizier_table`/`fetched_utc` in the JSON; e_RM is now fetched at all.
6. **The correlated fixture exists and the test can fail**: the suite imports
   `rmstructure.synthetic_rm_screen` and asserts the block jackknife exceeds the i.i.d.
   bootstrap on a field with a known coherence scale — the failure the deterministic csc|b|
   fixture is structurally unable to expose (stated in Methods).
7. **The Aitoff longitude labels are fixed** (relabelled with (−x) mod 360, so the printed
   number is the true l; the caption says which way l increases).

MINORs: the latitude profile carries per-bin jackknife errors (JSON + error bars in the
figure); robustness variants committed and each free to move the answer (5°/70° edges → 7.5,
|RM|<300 → 5.2, e_RM-below-median → 6.3 — the ratio is bin-definition-dependent through the
csc geometry and the paper says so); dedup measured (2 pairs within 5″); the λ² clause now
says the helper anchors the stack but is not applied to the pre-fitted catalogue values;
display precision matched (5.4 ± 0.6); `\rmTruth` no longer ships in the macro file. Macro
namespacing (finding 13) deliberately deferred: both guards defend the file and the referee
marked it low-urgency.
