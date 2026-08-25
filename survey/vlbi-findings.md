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
which lack a boosted core and serve as **steady negative controls**. CSOs are established radio flux
calibrators (Taylor & Peck 2003 measured ~0.7% rms over 8 months for a clean sample), but two of ours
are imperfect: **OQ 208 is documented as atypically variable** (40–60% in components; Wu et al. 2013)
and **2021+614 shows long-term cm variability** (Taylor et al. 2000), so we treat the floor's
sensitivity to them explicitly below; 0108+388 and NGC 3894 are the cleaner controls. X-band light
curves span 1995–2022 with 4–123 epochs each.

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

The 4 CSO controls are **also** formally "variable" by the absolute χ² test — even the steadiest, NGC
3894, has η = 21, p ≈ 0. For a genuinely steady source this apparent scatter is *not* intrinsic: it is
amplitude calibration + (u,v)-coverage / resolved-structure differences between the heterogeneous
geodetic sessions (the baseline set, and hence the (u,v) sampling, changes session to session). This is
the signal that **the assumed 5% error is far too small** — the controls imply an effective per-session
scatter of ~19%, so absolute η/χ² is **uninformative as a discriminant**: against a 5%-too-small error
*every* source, steady or not, comes out "significant". (For OQ 208 specifically, η = 80 partly reflects
genuine variability — Wu et al. 2013 — so it is the least clean control.)

**What survives is the amplitude V, benchmarked against the controls.** The CSOs set an empirical floor
**V ≈ 0.19** (4-control median 0.193; control V range 0.183–0.250); above it sit **13 of the 14**
non-control AGN. The floor is **robust to the imperfect controls**: dropping OQ 208 gives 0.186, and
dropping both OQ 208 and 2021+614 gives 0.185 — the same 13/14 pass either way (the marginal source,
CTA 102 at V = 0.193, stays above). The boosted blazars have median V = 0.32 vs the controls' 0.19
(roughly twice as variable in amplitude) — but this ratio describes *this curated sample*, not the
blazar population, and the floor rests on only 4 controls. The one non-control below the floor is **3C
273** (V = 0.13, only 4 epochs): with so few epochs and differing (u,v) coverage, plus its extended
S/X jet, its core variability is not captured — an honest small-N edge case, not a failure.

## Honest assessment

- **The method works; the absolute significance does not.** The recover-a-known is clean (known
  blazars rank top; CSOs are the least variable), but η/χ² is uninformative here because the per-session
  flux error is not 5% — it is ~19% (the CSO floor). **V relative to a steady-control floor is the
  trustworthy discriminant**, not absolute η. The tool reports both and the figure draws the floor.
- **VLBI total flux ≠ intrinsic flux.** `Fl_int` folds (u,v) coverage and resolved structure; the CSO
  controls quantify how large that structural floor is (V ≈ 0.19). Any "variability" below it is noise.
- **η/V are X-band only;** S-band (sparser, less precisely calibrated) is used only for the index.
- **α$_{SX}$ mixes spectrum and resolution.** S (2.3 GHz) and X (8.4 GHz) correlated fluxes sample
  different angular scales (resolution differs ~3.6×), so α$_{SX}$ reflects both the intrinsic spectrum
  and differential resolution; small for core-dominated blazars, but interpret partially-resolved
  objects (3C 84, 3C 120) cautiously.
- **Imperfect controls / small N.** The floor is the median of 4 controls (range 0.183–0.250), two of
  which (OQ 208, 2021+614) are themselves mildly variable; the result is robust to dropping them but the
  floor is not tightly constrained.
- **A curated validation set, not a survey.** 18 hand-picked sources demonstrate the pipeline and the
  recover-a-known; a blind catalogue would need a large, magnitude-limited input list and the
  control-floor calibration applied per declination/epoch-count bin.
- **Reproducible:** `python -m jansky_research.vlbi --online` regenerates `results/vlbi_candidates.csv`,
  the metrics JSON, and the η–V figure from the public Astrogeo files.

## Referee round on the style conversion (2026-08-24)

**Restyle findings, fixed.** The emphasis audit passed (all nine `\emph{}` spans survive as plain
text, none deleted) and the caveat audit did not: deletions occurred *adjacent* to three of them and
each removed a qualification.

1. **MAJOR.** "We stress that $f=0.05$ is only a common starting assumption, \emph{not} a value
   prescribed for these heterogeneous networks" had become "We note ... not ...". Two stress markers
   came off the same caveat in one edit. This is the paper's only advance warning that the 5%
   per-point error is an assumption, in a paper whose central result is that it fails by ~4x
   (control median `v_floor` = 0.193 against the assumed 0.05). Restored.
2. The abstract lost "The central methodological point is that", so "much larger than a naive
   few-percent floor" stood as a free-floating claim about VLBI in general rather than something
   this paper measured.
3. A comma made the research question's condition non-restrictive: "can one recover the known
   variability of well-studied AGN**,** with the dominant VLBI systematic honestly handled?"
   presupposes the honest handling instead of asking whether recovery survives it.
4. "not a discovery claim" was demoted to the third item of a comma series; promoted to its own
   sentence, which leaves it stronger than it was before the conversion.

## Two pre-existing defects, fixed here

**1. "The most variable source is OJ 287" ranked by the statistic the paper declares unusable.**
`\viTopName` is `argmax(eta)`. The paper's own thesis is that absolute eta is not a usable
discriminant and that V relative to the control floor is. On V, BL Lac (J2202+4216) leads at 0.540
and OJ 287 (J0854+2006) is second at 0.480 --- so "the other archetypal blazars (BL Lac, ...) follow"
was false on the preferred statistic. Reworded to "the most significantly variable", with the
V-ranking stated.

**2. The 13/14 recovery count turns on a tie at the floor.** J2232+1143 has V = 0.193, identical to
`v_floor` = 0.193 at the quoted precision, and is scored `above_floor = 1`. Now stated as twelve
unambiguous plus one at the floor.

**Still open (needs a decision, not a re-run):** `\viInjected`, `\viCompleteness` and `\viPurity` are
`--` in the macro file but unused in `main.tex`, so the Methods sentence "recovered at high
completeness and purity" is an unquantified claim. Either cite those macros or drop the adjectives.

## Full referee round (2026-08-25): MAJOR REVISION, 15 findings, one BLOCKER

The negative half (absolute η is unusable on Astrogeo fluxes) is correct and well-supported;
the macro chain is clean and the figure matches the CSV row-for-row. The positive half — the
control floor — is quoted as a point estimate whose every property the referee measured from
the committed CSV.

**BLOCKER: `pecktaylor2001` is by Fassnacht & Taylor** (astro-ph/0106001 → AJ 122, 1661, DOI
10.1086/322112); Peck was never an author, the compiled .bbl renders "Taylor & Peck (2001)",
and the findings file adds a third attribution ("Taylor & Peck 2003"). It is the ONLY support
for the sentence that licenses the design (CSOs as ~0.7% calibrators) — and that figure is an
8-month campaign, supporting "steady on ~year timescales", not decades.

**MAJORs (all measured from results/vlbi_candidates.csv):**
- The floor is the control MEDIAN, so it flags a steady source 50% of the time by
  construction — two of the four controls sit above their own floor; at the control MAXIMUM
  (0.250) the headline count is 10 of 14, and the paper never states the selection function.
- The "robust to imperfect controls" check only ever dropped the two HIGHEST controls (both
  perturbations can only add sources); the full drop-one jackknife spans floor 0.186–0.201 and
  count 13→12; over the control range the count runs 13→10; with V's sampling scatter only 4
  of the 13 are above the floor at z>3 (CTA 102 sits at z = −0.01).
- V is confounded with epoch count: Spearman ρ = 0.74 (p = 0.004); OLS V = 0.074 +
  0.170·log10(n) spans the sample's whole dynamic range; the floor is calibrated on n = 4–53
  and applied to n = 4–123; an epoch-matched floor gives 12/14. Filed as future work, belongs
  in the limitations.
- The blazar-vs-control amplitude comparison cannot fail: median_v_variable is the median of
  sources SELECTED by v > floor where floor IS median_v_control — guaranteed for any input;
  the honest unselected non-control median is 0.308; \viFloor and \viMedVctrl are the same
  quantity presented as two numbers; README's "~1.7×" repeats it uncaveated.
- The offline fixture validates a DIFFERENT selector (the log η–log V outlier cut that returns
  zero candidates); the control-floor block is pragma-no-cover, so the paper's method has no
  measured completeness/purity — and Methods claims "recovered at high completeness and
  purity" citing macros that are "--" in the committed real file (flagged in this file
  previously, still unfixed).
- There is no "known" in the recover-a-known: no published amplitude is compared for any
  source; komossa2023 (MOMO VI: OJ 287 variability — the obvious anchor) sits UNCITED in
  refs.bib.
- \viNabove = 13 is not reproducible from the committed evidence: at the CSV's three-decimal
  precision CTA 102 (V = 0.193) TIES the floor (0.1935) — the Discussion admits "twelve
  unambiguous plus one tie" but the abstract states 13/14 flat.

**MINOR/NIT:** "Every number is written by the pipeline" is false (≃0.186, ≃0.185, ~19%, "four
times", 1995–2022 — the last in NO committed artifact); V is not noise-debiased though
vlass.debiased_modulation_index is computed and discarded (on controls it directly measures
how wrong f = 0.05 is); no source table mapping common names to J2000 rows (the tie source is
never named); the Kiehlmann ≲10% criterion is unverified and whether the four controls are in
the bona-fide catalogue was never checked; mode-dependent macros un-namespaced (\viN = 18 vs
400); figure x-label carries a literal "vs.\"; "p ≈ 0" overstates 3C 273's 1.6e-4.

**Status: fixes pending.** The single change: a floor that carries its own uncertainty
(jackknife range 0.186–0.201, control span 0.183–0.250), the median-threshold selection
function stated, and the headline count quoted as a range (13→10 across the span; 12 under an
epoch-matched floor) — every input is already in the committed CSV.
