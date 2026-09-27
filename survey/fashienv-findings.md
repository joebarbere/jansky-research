# Findings — the environment-split FASHI HI mass function (plan 45, DR1 first leg)

`jansky_research.fashienv` splits the FASHI HI mass function (HIMF) by large-scale environment
for the first time — void/wall (Douglass+2023 VoidFinder) and group/field (Tempel+2017) —
recovering the ALFALFA void-HIMF suppression from FAST data.

## GATE 0 (2026-07-08)

- **DR2 is NOT public yet.** arXiv:2606.31539 (156,411 sources) is out, but the catalogue links
  on zcp521.github.io/fashi 404 ("available upon publication"; release ~Aug 2026). Verified the
  404 directly. → **DR1 first leg**, the same DR1-while-DR2-embargoed pattern `rmstructure` used;
  DR2 swap is a one-line source change.
- **FASHI DR1 is on VizieR: `J/other/SCPMA/67.19511/table2`, 41,741 sources** (the ".19511" is
  the SCPMA article number, NOT the row count — the agent's 41,741 was right; my first pass
  mis-read it as 19,511 and I corrected it). Columns: RA/Dec, cz, z, W50, Ssum (flux), Dist,
  logMass — everything the HIMF needs, precomputed.
- **Novelty confirmed**: the DR2 paper and DR1 paper both publish only a GLOBAL HIMF; the sole
  environment study (arXiv:2510.22902) used 230 group galaxies. No DR2 environment-split HIMF /
  deficiency / void study exists.
- **Cross-match catalogues resolve** (all VizieR): Tempel+2017 groups `J/A+A/602/A100`
  (table2 = groups with R200/M200), Douglass+2023 voids `J/ApJS/265/7` (table1 = VoidFinder
  spheres, Planck2018 cosmology, Mpc/h coords — verified the frame matches standard RA/Dec).
- **Plan citation corrected**: Lim+2017 is MNRAS 470, 2982 (not ApJ 854, 62); not used in the
  end (Tempel groups suffice).

## Scope corrections vs plan 45 (both forced by what FASHI lacks)

- **"Gas fraction at fixed M*" DROPPED**: FASHI carries no stellar masses. The literature void
  gas-fraction excess is also weak/dwarf-only (Kreckel+2012), so void-vs-wall HIMF is the
  cleaner statement — exactly what Moorman+2014 / Jones+2018 did.
- **"HI-deficiency vs clustercentric radius" DROPPED**: classical deficiency needs optical
  diameters/types (absent), and the raw median-HI-vs-R/R200 of DETECTED sources is
  selection-biased (stripped galaxies drop out of a flux-limited sample). Replaced with the
  cleaner group-member vs field HIMF split.

## Recover-a-known (offline, in CI)

- Injected two Schechter HIMFs (void: logM*=9.70, α=−1.45; wall: 9.95, −1.25) into a
  flux-limited mock; the 1/Vmax + Schechter fit recovers both knees/slopes within 0.25 dex and
  the correct knee-offset sign. Single-Schechter round-trip also passes.

## Result (FASHI DR1, SDSS-cap overlap)

| environment | logM* | α | N |
|---|---|---|---|
| global (all DR1) | 9.94 | −1.73 | 41,741 |
| void | 9.70 ± 0.08 | −1.6 | 2,538 |
| wall | 9.95 ± 0.05 | −1.7 | (cap) |
| group member | 10.08 ± 0.05 | −1.74 | 6,119 |
| field | 9.89 ± 0.05 | −1.68 | (cap) |

(void/wall use the Planck2018-matched void geometry; group/field use Tempel R200 membership)

- **Void knee offset = −0.256 ± 0.087 dex = 2.93σ** (void suppressed vs wall) — consistent in
  sign, somewhat larger, than Moorman+2014's ALFALFA void HIMF (−0.14 dex, ~2σ). **The headline:
  an independent FAST-based measurement.** Void faint-end marginally flatter (Moorman-like), not
  steeper (Jones-like); within errors at DR1 size — report the sign only.
- **Group knee offset = +0.19 ± 0.073 dex = 2.64σ** (members *higher* than field) — statistically
  the STRONGER offset, but NOT the headline because it's survivor-biased: a flux-limited HI
  survey detects the group galaxies that RETAINED HI (the gas-rich survivors, group knee 10.08 >
  global 9.94); stripped members drop out. So the group knee is an UPPER ENVELOPE and the true
  stripping effect is more negative. The void–wall comparison avoids this (voids don't strip).
- **Distance-model sensitivity (GATE-2 catch):** the void offset depends on the comoving-distance
  frame. Using the EdS (q0=+0.5) relation dilutes it to ~1σ; using the Douglass catalogue's own
  Planck2018 cosmology (q0=−0.527, so the frames match) gives the −0.256 dex / 2.9σ above. The
  matched-cosmology geometry is the correct choice and is what we report.

## Honest caveats (GATE-2 material)

- **Absolute α (~−1.73) is steeper than the published FASHI global (~−1.3)**: the simple 1/Vmax
  uses a single flux limit, not FASHI's full completeness function. The RELATIVE knee offsets
  (same estimator both bins) are robust; the absolute slope is not, and we do NOT claim to
  reproduce the global HIMF.
- **Footprint**: void/group cross-match is the SDSS-cap subset (~27k classifiable of 41,741),
  not the full survey — FASHI's southern/high-z sky has no SDSS optical catalogue.
- **Single void-finder**: VoidFinder (Douglass table1 spheres) only; the V²/VIDE/REVOLVER
  galaxy-membership tables (keyed by NSAID) need an NSA cross-match — a stated follow-on for the
  "report per void-finder" robustness the plan wants.

## Reproduce

`uv run python -m jansky_research.fashienv --out .` (VizieR fetches, ~min). Offline CI leg:
`--offline`. DR2 swap: point `fetch_fashi_dr1` at the DR2 table when it publishes (~Aug 2026).

## Referee round (2026-08-12) — the variance was fine; the biases are not

**A test that could have failed, and didn't.** The referee expected a void jackknife to kill
the 2.93σ offset, on the grounds that the quoted error is Poisson counting noise within one
realisation and cannot see void-to-void variance. Measured by deleting each of the 186
occupied voids in turn and refitting: **jackknife error 0.039 dex**, *smaller* than the fit
error of 0.087, so the offset strengthens to 6.5σ on that axis. A bootstrap over voids agrees
(0.031). Sample variance across voids is **not** what limits this measurement. Both numbers
are now committed; the expectation was wrong and the paper says so.

**What does limit it is bias, which no resampling can see.** Two, both acting in the direction
of the reported signal, and both now stated where the reader meets them:

- *1/V_max assumes uniform density within the accessible volume* — exactly false for a sample
  selected *by* density. For the void bin the accessible volume is the union of void spheres,
  whose fraction falls to zero at the near and far edges of the box while staying finite for
  the wall bin; since V_max is mass-dependent, so is the distortion, and it suppresses φ at the
  high-mass end of the void bin. The offset is now framed as an upper bound on any true
  suppression, with SWML named as the standard fix and explicitly not attempted.
- *The "wall" bin is a Cartesian bounding box of a wedge-shaped survey.* Now quantified in the
  evidence file, which never carried it: **n_wall = 24,175 — 58% of all of DR1** — with a knee
  of 9.954 against the all-sky 9.944, i.e. **0.010 dex apart, a fifth of its own error**. The
  comparison is closer to *void versus everything* than to void versus wall.

**Numbers that were hard-typed and wrong.** The paper quoted the void knee error as ±0.08
against a committed 0.071, and the wall as ±0.05; combined those give 0.094, which does not
reproduce the quoted 0.087 (√(0.071²+0.05²) = 0.0868 does). Both are macros now. The claim
that "every number above is pipeline-generated" is narrowed to the data-derived ones, since
literature values, the flux limit and catalogue sizes are hard-typed.

**Independence overstated.** FAST against Arecibo is an independent HI measurement, but the
*environment* half of both comes from the same SDSS DR7 parent and largely the same voids —
a deeper HI sample over much the same volume, not an independent realisation of large-scale
structure.

Open and reported, not fixed: the single void-finder headline (plan 45 says never to quote
one); beam confusion, which plan 45 required and which reproduces the sign of both offsets;
the injection validation, which places galaxies uniformly in comoving volume and so switches
off the one systematic that matters; the α–M* covariance, discarded by keeping only
`sqrt(diag(pcov))`; `curve_fit` without `absolute_sigma=True`; `n_bins` silently dropped from
every metrics dict by an `isinstance(v, float)` filter; the EdS distance-frame variant quoted
as "~1σ" with no committed number; and `FASHI_FLUX_LIMIT = 0.30`, which sets every V_max and
appears nowhere in the paper.

## DR2 leg (2026-09-26) — option B weighting; the paper's headline number changes

DR2 is public (not on VizieR): Table 2 CSV from the CSTCloud share linked at
zcp521.github.io/fashi.html (`fetch_fashi_dr2`; two POST calls, signed URL). 156,411 sources;
DR1 is ~95% contained (position + velocity match). DR2's `distance` convention differs from DR1's
(matched sources ~2.5% farther), so absolute masses are not interchangeable across releases.
286 sources with z <= 0 are excluded (no comoving position).

**Weighting changed (owner decision, "option B"):** each galaxy weighted by 1/(C * Vmax) with
DR2's own per-source completeness C and Vmax (over 19,482 deg^2), on the DR2 HIMF sample
C >= 0.5 ("above the 50% flux completeness limit"; 110,520 sources vs the paper's ">109,000" —
DR2 does not publish the W20/f_sigma cut columns). This replaces the single 0.30 Jy km/s flux
cut the referee flagged. The old weighting on the same DR2 data is kept as `optA_*`, so the
method change and the sample change are visible separately.

| | DR1 committed (flux cut) | DR2 A (flux cut) | DR2 B (1/(C Vmax)) |
|---|---|---|---|
| N | 41,741 | 155,983 | 110,520 |
| global log M* / alpha | 9.944 / -1.732 | 10.048 / -1.780 | **9.907 +/- 0.021 / -1.328 +/- 0.032** |
| void - wall knee | -0.256 +/- 0.087 (2.9 sig) | -0.252 +/- 0.092 (2.7) | **-0.154 +/- 0.052 (3.0)** |
| void jackknife err | 0.039 | — | 0.041 (409 occupied voids) |
| group - field knee | +0.192 +/- 0.073 (2.6) | +0.208 +/- 0.068 (3.1) | **+0.141 +/- 0.039 (3.6)** |

1. **Recover-a-known passes:** B's global HIMF reproduces the DR2 paper's own
   (log M* = 9.89 +/- 0.02, alpha = -1.31 +/- 0.02). The single-cut weighting's alpha ~ -1.75 was
   the bias the paper already disclosed; it persists on DR2 (A), so it was the method.
2. **The environment effects survive the corrected weighting but shrink:** void offset by 40%,
   group offset by 27%. A barely moves vs DR1, so the change is the weighting, not N. Part of
   the published -0.256 was a weighting artefact. B's -0.154 sits inside the Moorman+2014
   0.1-0.2 dex range.
3. Still not addressed by B: the void-volume Vmax bias the referee raised (the per-source Vmax is
   survey-wide, not restricted to the void/wall volume), and the SDSS-cap footprint (the
   classifiable region is the void catalogue's bounding box).

Evidence: `results/fashienv_dr2_preview.json` (its own file until the paper is revised; the
committed `fashienv_metrics.json` and macros still hold the DR1 numbers under the DR1 prose).
**Next:** revise `papers/fashienv/` to DR2 + option B, then a presenter/referee round — the
headline number changes, so this is a revision, not a data refresh.

### Revision (2026-09-26, same day): DR2 + option B is now the paper; the void offset is fragile

The paper now reports DR2 with the 1/(C Vmax) weighting; `results/fashienv_metrics.json` is the
DR2 run (the preview file is retired; the DR1 file is in git history). Every comparison the prose
makes is pipeline-generated: the old weighting on DR2 (`optA_*`) and on DR1 (`dr1_optA_*`), the
DR1<->DR2 match (94.7% matched; DR2 distances larger by a median factor 1.025), the share of the
old offset due to the weighting (40%), and the wall bin's share of the sample (61%). Hand-typed
DR1-era numbers (58%, 0.010 dex, a void bootstrap of 0.031) were removed or replaced.

**The Einstein-de Sitter check, recomputed under option B, removes the void signal:** placing
galaxies with q0 = 0.5 instead of the Douglass catalogue's Planck q0 = -0.527 (a few per cent in
distance at this depth) takes the void-wall offset from -0.154 (3.0 sigma) to **-0.021 (0.43
sigma)**. The DR1 paper had reported "dilutes to ~1 sigma" in a subordinate clause. The matched
cosmology is the right one, but a signal that a few-per-cent distance change erases rests on
galaxies near void boundaries -- a fragility the fit error and the jackknife (both of which hold
the geometry fixed) cannot see. The abstract now calls the void offset a tentative upper bound,
not a detection, and names this test. Group-field (+0.141, 3.6 sigma) is unaffected by it but
remains survivor-biased.

## Referee round on the DR2 revision (2026-09-26): MAJOR REVISION, 12 findings — response

The referee checked every macro against the JSON (all agree) and found the problems were
interpretation, not bookkeeping. Two of them were my errors from earlier the same day, and are
recorded as such.

**Retracted (mine):**
- *"The EdS check shows the void offset is fragile."* Wrong. For the distant, knee-mass galaxies
  the EdS relation moves positions by 0.4-1.2 hole radii; measured, it changes the void status
  of 7% of void galaxies at z < 0.02 and of more galaxies than the void bin contains at z > 0.06
  (126%). It scrambles membership, so it would remove a real signal too. It is now reported as
  exactly that, not as evidence.
- *"The single-limit weighting gives -0.252 on the same data" / "40% of the offset is the
  weighting."* Not the same data: option A had run on all 155,983 sources, not the C >= 0.5
  sample, and the 40% compared DR1-A with DR2-B (catalogue, sample and weighting at once). Now
  measured on the same sample with a paired void jackknife: the weighting alone moves the
  offset by +0.078 +/- 0.018 (4.3 sigma; -0.228 -> -0.150).
- *"Upper bound."* The paper gave the bounding-box bias opposite signs in two places, and the
  1/Vmax bias's sign was asserted. Dropped.

**New computations (all in `results/fashienv_metrics.json`):**
- **Full VoidFinder voids** (table2, all holes): 17,289 void galaxies vs 6,399 in maximal spheres;
  offset -0.150 +/- 0.043 vs -0.157 -- the classification barely matters.
- **Random-void null** (the referee's key recommendation): every void moved rigidly (distance,
  holes, size kept) to a random position inside the Tempel footprint, 100 placements, same
  weights. Mean **-0.091 +/- 0.018 (s.d.)**: the estimator and geometry produce ~60% of the
  measured offset, in the direction of the signal (the sign is now measured, not argued). No
  placement reached the measured value (one-sided p = 0.0099). **Excess from the real voids:
  -0.059 dex**, 3.3 sigma against the placement scatter, 1.8 sigma against the void jackknife
  (0.033). Tentative, not a detection.
- **Group frame fix:** Tempel `zcmb` converted to heliocentric with the Planck dipole (was
  Dist.c x H0=70 vs the catalogue's 67.8, in the wrong frame). Group offset +0.165 +/- 0.043
  (3.85 sigma).
- **Fit quality committed:** reduced chi2 19 (void), 26 (wall), 41 (global). curve_fit scales
  the errors by it; the knee values depend on the Schechter form.

**Claim changes:** survivor bias withdrawn as the group explanation (1/Vmax counts a stripped
galaxy at its current mass; it leaves the sample only far below the knee); the group offset is
now unexplained, with beam confusion (FAST ~3'), the field bin's out-of-footprint sky and
distance-dependent group finding listed. "Recover the published HIMF" is now "reproduces, as a
bookkeeping check". The mock-validation sentence says what the mock tests (the sign, with no
void geometry). The jackknife is read as "no excess void-to-void variance", not as a second
significance. The "every number is pipeline-generated" sentence is qualified, and the 142 other
dropped sources are stated. `fashi_dr2` now carries the Crossref-verified journal DOI
(10.1007/s11433-026-3072-1, Sci. China PMA 69, 129811).

**Headline now:** void-wall knee -0.150 +/- 0.043 dex, of which -0.091 is reproduced by randomly
placed voids; the real-void excess -0.059 dex is tentative. Group-field +0.165 +/- 0.043,
cause not established.

**Still open:** a density-insensitive estimator (SWML); a true footprint mask; a random-placement
null for the groups; a beam-confusion test; whether FASHI v_opt and the Douglass redshifts share
a velocity frame.

## Second referee round on the DR2 revision (2026-09-26): MAJOR REVISION

Bookkeeping now clean (every macro and every derived ratio re-derived from the JSON and agrees).
Of the 12 first-round findings: 9 RESOLVED, 3 PARTLY (#5 group reasoning; #10 `assign_groups`
still uses H0=70 against Tempel's h=0.678 for R200/distance, ~3%; #12 `fashi_dr2` fourth author
is Hong Guo per Crossref, not "W.-K. Guo" -- my "Crossref-verified" checked DOI/volume/article
but not every author). The headline now rests on the random-void null, and the referee argues
it has not been shown fair:

1. **The null's scatter (0.018) is below every noise estimate for the real void bin** (fit 0.043,
   jackknife 0.033) -- the expected symptom of placements landing on mean-density structure and
   filling a better-populated void bin. Against the real bin's own noise the excess is 1.4-1.8
   sigma, not 3.3. Needs per-placement occupancy committed and a jackknife on a few placements.
2. **The null is not signal-free or geometry-matched:** placed voids can overlap the real voids
   (leaking signal into the null), overlap each other, spill holes outside the footprint (only
   the centre is tested), and the classifiable box stays the real voids'. "Alone produce 60%"
   is stronger than shown. Needs a constrained null + an overlap regression.
3. **EdS (-0.063) sits at the null's 93rd percentile** -- a coherent few-Mpc shift moves the
   statistic by the size of the excess, a variation the placement scatter cannot see.
4-11 (minor): void fit uses 14 bins vs the wall's 18; beam confusion and the unresolved
   velocity frame missing from the void section ("attributable to the real voids" -> "not
   reproduced by random placements"); the DR1 offset moved -0.256 -> -0.135 with full voids
   (undisclosed) and "DR1 and DR2 agree" is wrong (0.095 dex apart); the group-section detection-
   limit argument is wrong (the limit reaches the knee at 380 Mpc); p = 0.0099 is just the
   0-of-100 floor (run >= 1000); run the null under the single-limit weighting too; ALFALFA
   comparison should be raw-vs-raw; caption should mention the null.
NIT: stale `weighting_shift_pct` (40.0) still in the JSON as a `_merge` carry-over.

**Single change recommended:** make the null demonstrably fair (per-placement occupancy,
real-void overlap, spill-over; constrained placements; >= 1000 draws), and until then lead the
abstract with the conservative 1.4-1.8 sigma.

## Null upgrade (2026-09-26, run `fashienv-round3`): the environment split does not survive

Round-2 referee asked for a demonstrably fair null. Built and run (1,000 unconstrained + 1,000
constrained void placements, per-placement diagnostics, a 10-placement jackknife on each, a
200-placement group null, a common-bin fit, DR1 under both void definitions). A disjoint
("no overlap") void null proved infeasible: 687 of 1,163 voids cannot be re-placed disjointly,
because the voids fill too much of the volume; overlap is measured instead.

- **Occupancy:** random void bins hold ~28,800 galaxies vs the real 17,289 -- placements land on
  mean-density sky, as the referee predicted.
- **But noise is not the explanation:** a delete-one-void jackknife on the placements gives
  0.028-0.042 (median ~0.034), the same as the real bin (0.033). The placement-to-placement
  scatter (0.017-0.018) is simply smaller than within-configuration void-to-void variance.
- **The null is not a clean bias estimate.** Null offset vs the fraction of null-void galaxies that
  are real-void members: slope +0.35 (constrained) / +0.42, intercept at zero overlap -0.20 /
  -0.22 -- *more negative than the measurement*. The geometric bias depends on where the placed
  voids sit, so no single null mean can be subtracted.
- **The null-subtracted excess depends on the weighting:** -0.063 under 1/(C Vmax), -0.160 under
  the single flux limit (null means -0.087 vs -0.068). The referee's can-fail test, failed.
- **Groups:** randomly placed groups give +0.126 +/- 0.056; 44 of 200 reach the measured +0.166
  (p ~ 0.22). The group offset is reproduced by geometry and estimator alone.
- Common mass bins: -0.134 +/- 0.043 (14 bins) vs -0.150; bin coverage is not the issue.
- DR1: -0.256 with maximal spheres, -0.135 with full voids -- to be disclosed.
- EdS sits at the 92nd percentile of the null.

**Conclusion:** the raw void-wall offset (-0.15) is robust as a *number*, but 1/Vmax
environment splits in this survey geometry produce offsets of the same size with no
environmental signal, that bias is environment-dependent, and the residual is
weighting-dependent. **No environmental dependence of the HIMF is claimed.** The group offset
is consistent with random placement. The defensible paper is a methodological caution: in a
FASHI-like geometry, random-placement nulls reproduce knee offsets as large as published void
effects (ALFALFA -0.14), so environment-split 1/Vmax HIMFs need such a null. Reframe pending
owner decision; the paper text still carries the previous framing.

## Third referee round, on the methodological-caution rewrite (2026-09-26): MAJOR (text-level)

Reframing endorsed; bookkeeping clean (every macro vs JSON, both figures); "worth publishing as a
short RNAAS-scale caution once #1-5 are fixed". Nearly all fixes are claim strength, not new
computation. Two are my misreadings of my own regression, recorded as such:

1. **I misread the overlap regression.** real_overlap_frac spans only 0.280-0.367 (5-95%:
   0.304-0.347); R^2 = 0.07; across the sampled range the fitted line moves 0.03 dex, not "~0.1
   dex". The zero-overlap intercept (-0.201) that I set beside the measurement is an
   extrapolation 3.2 range-widths out, and the same line predicts +0.15 at overlap 1 against a
   measured -0.15 -- the linear model fails at both ends.
2. **"Overlap" is a proxy for occupancy and redshift, with the wrong sign for leakage.** Offset
   vs occupancy + median z: R^2 = 0.200; adding overlap: 0.201, and its coefficient falls 0.351 ->
   0.056. A positive slope is the opposite of real-void signal leaking into the null.
3. **The Discussion's mechanism is contradicted by the rows:** within the null, emptier
   placements give LESS negative offsets; the real voids (17k members) lie outside the null's
   occupancy range (25-32k), so the null cannot say what a density-matched split would give.
4. "No environmental dependence survives" overclaims -- the void offset is beyond all 1,000
   placements under BOTH weightings; the supported statement is "cannot attribute". The intro
   also said the tests showed the DR1 -0.256 not attributable, but DR1 was never null-tested.
5. "Arbitrary volumes" / "no environmental signal" blame the null offsets on the estimator
   alone, contradicting Sec. 4.2 ("not signal-free"); the volumes keep the catalogues' distances,
   which is the portable point.
6-14 (minor/nit): Moorman+2014 may have used a density-insensitive estimator (2DSWML; to verify)
   so the ALFALFA juxtaposition needs narrowing; the single-limit weighting fails the global-HIMF
   check (alpha -1.73, chi2 190) so it should not carry half the argument; the group null
   records no occupancy/redshift diagnostics; 687/1163 is order- and budget-dependent and
   contradicts the blanket reproducibility sentence; "beyond all 1000" is conditional on one
   galaxy sample; Fig. 1 caption should state fitted bins; "126%" undefined; fashi_groups year is
   2026 not 2025 (Crossref); stale \feRealWeightingShiftPct macro.

**Recommended single change:** rebuild the "no unique baseline" argument on what the rows show
(no density match; offsets track occupancy and redshift; the preferred-weighting residual is
-0.063, 1.5-1.9 sigma against the real bin's noise) and conclude "we cannot attribute".

### Response (run `fashienv-round4`; every earlier number reproduced exactly)

All fourteen findings addressed in text. Three new diagnostics were computed rather than asserted:

- **Regression, recomputed in the pipeline** (`void_null_constrained.B_regression`): overlap
  0.280-0.367; R^2 overlap alone 0.072, occupancy + z 0.199, all three 0.200; overlap coefficient
  under controls 0.056; the fitted overlap line moves **0.031 dex** over the sampled range. The
  zero-overlap intercept is gone from the paper.
- **Group-null diagnostics:** placed group regions hold a median **1,949** galaxies at z 0.0225,
  against **18,696** at z 0.031 for the real groups, with a median 0.43 real-group overlap. So
  44/200 is a weak test, not evidence of no group effect. The paper now says "neither confirms
  nor excludes" (this phrase is the triage LOW "overclaim" hit, which is a false positive).
- **No-overlap trial in the pipeline** (seed 67, 500 tries): **688/1163** voids unplaced (the
  one-off test earlier gave 687, a different draw order). Reported as a single seeded trial, not a
  reproducible constant.

Text: "cannot attribute" replaces "no environmental dependence survives"; the volumes keep the
catalogues' radial selection; Moorman+2014 used 2DSWML (verified), so the ALFALFA contrast is
narrowed; DR1 -0.256 never null-tested; the single-limit weighting is demoted (global alpha -1.73,
reduced chi2 190); Fig. 1 caption states the fitted bins (void points below log M ~ 8 lie below the
curve, checked on the rendered figure); fashi_groups year 2026; stale WeightingShiftPct removed.

## Fourth referee round (2026-09-26): MINOR revision

Bookkeeping re-verified in full: every macro against the JSON, including the derived residuals
and sigmas, both figures, and four citations via Crossref. 13 of 14 third-round findings are
resolved; #3 is partial and #5 is fixed in the body but not the abstract. Nothing needs new data.

1. **(major) The recommendation is scoped too broadly.** Both environments are weighted by the
   survey-wide Vmax, never a Vmax restricted to the environment's own volume. This round-1 item
   was never closed (see the DR2 leg above). The null offset may be mostly that known
   misapplication. Scope the claim to survey-wide-Vmax splits, name the environment-restricted
   Vmax as the direct remedy, and list it as next work. Optional: compute the void-volume
   fraction vs distance with `void_membership_holes` on randoms, and refit.
2. **(major) The null's one directional hint is omitted.** corr(offset, n_in_void) = -0.42:
   emptier placements give LESS negative offsets, and the real voids (17k) are emptier than
   every placement. Extrapolating would enlarge the residual, not remove it. State the sign
   alongside the refusal to extrapolate. "Tracks" at R^2 = 0.2 should read "correlates weakly".
3. **(major, venue) Not RNAAS-sized:** about 1,590 body words and two figures against ~1,000
   and one.
4. Hedging has accumulated. Cut to the four-sentence contribution. Trim the EdS test, the
   regression detail, 688/1163 and the duplicated DR1 history. Drop the single-limit weighting
   from the abstract (it fails the global check) but keep the 0.078 weighting dependence once.
5. The sigma denominators are the offset's fit and jackknife errors, not "the void bin's"
   noise. Say why they are used rather than the null std.
6. The abstract gives the null offsets a single cause (the estimator); the Discussion gives
   two (estimator plus real structure).
7. "Redshift distribution is preserved" sits next to the member median-z mismatch (0.0225 vs
   0.031). Clarify.
8. Three-decimal false precision ("C >= 0.500", "51.300%", "91.700 percentile").
9. The Fig. 2 caption doesn't name the weighting.
10. `README.md:95,195` and the stale DR1 arXiv tarball still carry the -0.26 dex claim.

### Response: environment-restricted Vmax (run `fashienv-envvmax`; `env_vmax` block)

Round-4 #1 was run, not just scoped. Each side of each split was reweighted by
C * Vmax * g_E(dmax). Here g_E(D) is the fraction of the volume within D that environment E
occupies, from 2x10^6 randoms uniform in volume over the DR2 footprint, classified by the
galaxies' own rules (`env_vmax_offset`, `environment_volume_fraction`, `survey_randoms`). A
unit test shows the method's point: a distance-only "environment" with one HIMF gives >0.1 dex
under survey-wide Vmax and <0.05 dex restricted. The nulls reuse the same seeds, so the 200 void
and 200 group placements are the committed ones. The survey-wide offsets reproduce the committed
rows with **max |diff| = 0**.

| | survey-wide Vmax: measured / null | env-restricted Vmax: measured / null |
|---|---|---|
| void - wall | -0.150 / -0.087 +/- 0.017 | **-0.098 +/- 0.034 / -0.046 +/- 0.016**, 0/200 reach |
| group - field | +0.166 / +0.126 +/- 0.056 | **+0.088 +/- 0.029 / -0.025 +/- 0.049**, 3/200 reach |

- **The restriction removes the group null offset and halves the void one.** Groups' share of
  the classifiable volume falls from 3% within 50 Mpc to 1% beyond 200 Mpc, because Tempel groups
  are found in a flux-limited sample. That trend was the group bias. The void share is 30-44%
  with no trend, so the half of the void bias that remains is not a Vmax-bookkeeping effect.
  The referee's hypothesis ("the null offset may be mostly the misapplication") holds for groups
  and half-holds for voids.
- **The void residual barely moves:** -0.063 before, -0.052 now (1.5 sigma against the fit
  error). Still not attributable.
- **The group offset changes status:** within the null before (44/200), now in its tail
  (3/200, p ~ 0.02; +0.113 dex excess). Not claimed, because the null regions hold a tenth of
  the real groups' galaxies. FAST-beam blending also grows with density and would raise the
  group knee, which is the sign seen. That makes it the first systematic to test.
- The paper is rebuilt around this as an RNAAS-length note: 868 body words, one figure (both
  weightings on the same 200 placements).
  - #2: added the occupancy direction (corr -0.42; emptier gives less negative; not
    extrapolated), and "tracks" became "correlates weakly".
  - #4: the hedging is cut.
  - #5: denominators are named, with the reason.
  - #6: the abstract now has the two-cause wording.
  - #7: the redshift sentence is gone.
  - #8: percent macros are formatted.
  - #9: the caption names both weightings.
  - #10: the README rows are updated. The gitignored DR1 arXiv tarball is stale and must be
    rebuilt (`make arxiv`) before any submission.
- Fig. 1 (void/wall HIMF) is no longer in the paper; `figures/fashienv.pdf` is still generated.

## Fifth referee round (2026-09-26): MINOR revision

The env-restricted Vmax implementation was verified wherever it could be checked:
- The DR2 Vmax inverts exactly to FASHI's comoving distance, so there is no lower distance limit.
- Galaxies and randoms use the same classifiable box and group-cylinder rules.
- All 46 macros match the JSON; the figure matches its rows; Crossref confirms both FASHI bib
  entries.

9 of 10 round-4 findings are resolved; #5 is partial. The problems are interpretive, and all
lean the same way:

1. **(major) "Without a trend" is false.** The void share of classifiable volume per shell is
   0.30 within 50 Mpc and ~0.42 beyond. Cumulative g_void/g_wall goes 0.44 -> 0.72. That trend
   is WHY the restriction moved the void numbers (constant g cannot move a knee, per the
   docstring). **My inference in the round-4 response was therefore wrong:** "no trend, so the
   remaining void bias is not Vmax bookkeeping" rested on a false premise.
2. **(major) "Necessary but not sufficient" / "does not remove the void bias" assume the
   remaining null offset is bias.** The null contains real structure, so part of it may be
   signal. "Necessary" is not shown either: the null-subtracted void residual was -0.063 before
   the restriction and -0.052 after (a 0.011 change). The Discussion's opening is single-cause
   again (round-4 #6 regressed).
3. **(major) The group tail status is created by the weighting change.** The restriction moved
   the group null by -0.151 but the measurement by only -0.078. Under survey-wide weighting the
   residual was +0.040 (44/200); under the restriction it is +0.113 (3/200), so two-thirds of
   the excess is the differential response. The env null mean is -0.025 +/- 0.003, not zero:
   "removes" should read "reduces to". The paired rows do not explain the differential response
   (r = 0.00 with occupancy). The real groups lie outside the null on z, occupancy and overlap.
   The blending caveat has the right sign, and would also depress the void knee.
4. Use one significance yardstick for both environments. In quadrature: void 1.4 sigma, group
   2.0 sigma. Two tests make p ~ 0.02 into ~ 0.04.
5. The distance frames are mismatched: FASHI comoving vs h70 has median ratio 0.951, and 10.4%
   of galaxies have dmax below their own h70 distance. The effect is estimated at << 0.01 dex
   but not measured. Rerun with dmax rescaled into the randoms' frame, in a worktree.
6. The occupancy correlation -0.42 comes from the survey-wide rows. On the paired 200 rows it is
   -0.39 survey-wide and -0.22 restricted. Group diagnostics are collinear and give no direction.
7-9 (nits): write "first 200 of the same seeded placements"; the runner writes the results
   JSON directly with no `--out`; "41.492" should be "41"; the figure legend overlaps the void
   histograms.

### Response (run `fashienv-envvmax2`, frame-corrected; pairing checks still exactly 0)

- **#5, the distance frame: measured, and small.** Each galaxy's dmax is now rescaled into the
  randoms' h70 frame (median ratio 1.052, 5-95% 0.958-1.156; 10% of galaxies had dmax below
  their own h70 distance uncorrected). The measured offsets move by <= 0.002 dex; uncorrected
  values are kept in `void_frame_uncorrected`/`group_frame_uncorrected`. The null means moved
  by <= 0.006, which mixes the correction with a fresh random draw (the randoms' outer radius
  changed). The previous run's summaries are kept in `env_vmax.previous`.
  - Void: measured -0.092 +/- 0.035, null -0.040 +/- 0.017, 0/200 reach it.
  - Group: measured +0.089 +/- 0.032, null -0.019 +/- 0.053 (SE 0.004), 4/200 reach it.
- **#1 accepted.** The text now says the void share rises from 30% within 50 Mpc to 42%
  overall. The retracted inference is gone.
- **#2 accepted.** "Necessary but not sufficient" and "does not remove the void bias" are gone.
  The Discussion now says the restriction removes the part of the null shift caused by the
  environment's volume share changing with distance. What remains is the estimator's response
  plus real structure, which the null cannot separate. The void residual barely depends on the
  weighting (-0.063 vs -0.052).
- **#3 accepted.** The abstract and results say the group residual grows from +0.040 (44/200)
  to +0.108 because the restriction moves the group null about twice as far as the measurement
  (-0.145 vs -0.077). "A group signal that appears under only one weighting is not yet a
  measurement." "Removes" is now "moves the null means to".
- **#4.** One yardstick for both (`excess_sigma_quadrature`): void 1.34 sigma, group 1.75
  sigma; rank p ~ 0.02, 0.05 for two environments.
- **#6.** Correlations are computed on the paired rows under both weightings (-0.22 restricted,
  -0.39 survey-wide). The group diagnostics are stated as collinear, with no direction.
- **#7-9.** "First 200 of the same seeded placements"; the runner has `--out`. Nit 8 was
  already handled: `write_results` merges through `preserve_live_results`. Reduced chi2 is now
  an integer; the figure legend has headroom. The blending caveat now names both signs (it
  raises the group knee and lowers the void knee).
- About 960 prose words plus the caption; one figure; triage and lint clean.

## Post-round-5 tests (2026-09-27): occupancy/redshift-matched null + beam blending (`robustness`)

Both were named as next steps in the paper. Run `fashienv-robust2`; code `label_shuffle_null`,
`confusion_counts`, `robustness_leg`; runner `scripts/fashienv_robustness.py --out`.

**Shuffle null.** Environment labels are given to the real number of members per dz = 0.0025
bin, drawn at random from the classifiable sample, 1,000 replicates. This matches occupancy and
redshift exactly, which random placement could not, but has no spatial coherence. The real
members are included at the base rate (void 23%, group 25%, against 32% overlap for the
placements), which dilutes the null toward the signal.

| | survey Vmax: meas / null / excess | env Vmax: meas / null / excess |
|---|---|---|
| void | -0.150 / -0.043 +/- 0.014 / **-0.107 (2.4 sigma, 0/1000)** | -0.092 / +0.027 +/- 0.014 / **-0.119 (3.2 sigma, 0/1000)** |
| group | +0.166 / +0.078 +/- 0.012 / **+0.088 (2.0 sigma, 0/1000)** | +0.089 / -0.006 +/- 0.012 / **+0.095 (2.8 sigma, 0/1000)** |

(sigma = excess / quadrature of the fit error and the null spread.)

- Under this null both residuals are larger than under random placement, and neither depends on
  the weighting. The group residual's weighting-dependence (round 5, #3) was a property of the
  random-placement null, not of the data.
- The two nulls disagree about the baseline. Survey-wide void means: -0.087 (placement) vs
  -0.043 (shuffle). The placement offsets correlate with median z (r = 0.31) and occupancy, both
  of which the placements get wrong and the shuffle matches. Placements carry spatial coherence
  (structure at average density); the shuffle does not. Nothing here decides which baseline is
  "right". They bracket the void residual at -0.05 to -0.12 and the group residual at
  +0.04 to +0.11.

**Blending**, within the SDSS footprint (55,893 classifiable sources). A source is "confused" if
2 or more Tempel SDSS galaxies (r < 17.77) lie within 2.9' and the velocity window.

| | confused in / out | env offset all -> unconfused | env shuffle excess, unconfused |
|---|---|---|---|
| group, W50/2+100 window | **25.5% / 1.6%** | +0.123 -> **+0.038** | +0.064 (1.8 sigma) |
| group, fixed +/-300 km/s | 29.1% / 1.6% | +0.123 -> **+0.044** | +0.079 (2.2 sigma) |
| void, W50/2+100 window | 6.8% / 10.8% | -0.088 -> -0.078 | -0.100 (2.6 sigma) |
| void, fixed +/-300 km/s | 7.5% / 12.2% | -0.088 -> -0.080 | -0.102 (2.6 sigma) |

- **The mass-independent window gives the same answer.** The drop in the group offset is not
  caused by W50 scaling with mass.
- **But the flag is mass-selective under either window:** median log M is 9.9 for confused
  members against 9.47 for unconfused ones. A source needs an SDSS counterpart plus a bright
  neighbour, both likelier for distant, massive, optically bright galaxies, and the beam covers
  more physical area at larger distance. Removing confused sources therefore removes massive
  group members whether or not their HI is blended. **Excluding them cannot separate blending
  from "massive galaxies have more neighbours".** The removal is a selection, the dr20radio
  censoring lesson again. What the test does establish: a quarter of group members are
  candidates for blending, and without them most of the group offset goes (+0.12 -> +0.04).
- **The void offset is insensitive to the exclusion** (-0.088 -> -0.078/-0.080) and stays 2.6
  sigma beyond the shuffle null.

**Where this leaves the claims (not yet in the paper):**
- **Void.** 2.4-3.2 sigma against the occupancy/redshift-matched null, robust to blending
  exclusion; 1.3-1.5 sigma against coherent random placement. The honest statement names both
  baselines and does not choose between them.
- **Group.** Not attributable. It is carried by the quarter of members that are blending
  candidates, and exclusion cannot tell blending from selection. Resolving it needs
  interferometric HI or a forward model of FAST-beam confusion.

### Reframing (2026-09-27)

The paper is rebuilt around the two nulls and the blending test. It is RNAAS-length: about 920
prose words, a 158-word abstract and one figure. The figure now adds the shuffle nulls as
mean +/- 1 sigma bands.
- **Void:** "reported as 0.05-0.12 dex beyond the two nulls rather than as a detection". The
  argument: the shuffle matches occupancy and redshift but not coherence; placement keeps
  coherence but not occupancy. A void is defined by being emptier than any random placement, so
  no placement can match both.
- **Group:** "carried largely by the quarter of members with a neighbour inside the beam";
  exclusion cannot separate blending from mass selection, so it is not attributed.
  Interferometric HI or a confusion forward model would decide it.
- **Recommendation:** for single-dish environment-split 1/Vmax HIMFs, use both nulls plus a
  confusion test, with each null's occupancy reported.
- **Removed from the text** (still in the JSON): the frame-ratio detail and the per-shell
  volume shares. The weighting-dependence of the group residual against placement now survives
  only as the two placement excesses quoted side by side.

## Sixth referee round, on the reframing (2026-09-27): MAJOR revision

The conclusions are pitched at about the right strength. Macros, figure and citations are
correct (Crossref: moorman2014, douglass2023, tempel2017, fashi_groups). The evidence under the
headline has three holes and one mis-stated bracket.

1. **(major) A quarter of wall/field galaxies have no environment classification.**
   "Classifiable" is a padded box, not the SDSS sky: 74,024 in the box against 55,893 in the
   SDSS cells. All 18,131 outside galaxies (24.5%) are labelled wall AND field. The headline
   offsets and both nulls use the box pool; the blending test uses the footprint pool. The
   offsets differ: group +0.089 (box) vs +0.123 (footprint); void survey -0.150 vs -0.171. The
   paper quotes both as "the restricted group offset". The shuffle excess before (+0.095, box)
   and after (+0.064, footprint) confused-source removal differs in two ways at once. This
   round-1 item was never closed. **Fix:** move everything to the in-footprint pool.
2. **(major) The shuffle is conditioned on redshift, not depth.** DR2's rms spans 0.48-2.45
   mJy (5-95%), and the mass floor at fixed z follows it. The docstring claim "mass is the only
   thing that can differ at fixed z" is false. Depth differences between void regions and the
   rest would inflate the 3.2 sigma end, and could explain the gap between the nulls (placements
   keep sky position; shuffles do not). **Fix:** stratify by (z, rms) or (z, sky cell).
3. **(major) The blending flag largely restates the Tempel group definition.** Two SDSS galaxies
   within 2.9' (~0.1 Mpc) and ~300 km/s are almost necessarily FoF-linked, so 25.5% vs 1.6% is
   largely by construction. Exclusion also removes close pairs, so it cannot separate blending
   from mass selection OR from close-pair environment. **Sharper tests with existing data:**
   W50 and W20/W50 residuals of flagged vs unflagged sources at fixed mass and z (blending
   broadens lines); dependence on beam radius; HI mass at fixed optical luminosity. A random
   exclusion matched in z and mass would be tautological.
4. **(major) The bracket is quoted without significances, and its low end (1.3 sigma) is
   consistent with zero.** "Bracketed by the two nulls" also misstates what is bracketed: it is
   the excess, not the offset. **Fix:** "an excess of 0.05 dex (1.3 sigma, placements) to 0.12
   dex (3.2 sigma, shuffles)".
5. "Real members pull each null towards the measurement" is wrong for the shuffle, which puts
   members in and out at the same rate. It holds for placements only (0.325 in vs 0.175 out).
6. The post-exclusion sigmas use the full-sample fit error, not the clean sample's, so they are
   overstated.
7. The round-4 occupancy sign regressed out of the text. Emptier placements give less negative
   offsets, so the placement end of the bracket is the conservative one; restore "weakly".
8. "Reproduce most of it" holds for one weighting only (58% survey-wide, 44% restricted); say
   "about half".
9. "Does not depend on beam confusion" overstates a test that is weak for voids (7-11% flagged;
   blind to faint companions). Say "insensitive to excluding sources with an SDSS-bright
   companion".
10. The figure's shuffle bands carry no measured-offset error, so it reads as ~8 sigma.
11. The shuffle unit test cannot fail. Add a depth-confound case.
12. Nits: ConfIn/ConfOut precision mismatch; 2.5 vs 2.6; "emptier than a typical randomly placed
    volume"; the title still says one null.

On "no placement can match occupancy": the argument is sound for rigid placements. The feasible
improvement is a (z, depth)-conditioned shuffle, not a better placement.

### Response (run `fashienv-round6`: the full real leg, everything on the in-footprint pool)

- **#1 fixed at the source.** "Classifiable" is now the padded box AND the SDSS footprint (1 deg
  cells holding Tempel galaxies) for galaxies and randoms alike (`in_sdss_footprint`). The pool
  is 55,893; all nulls, weightings and the blending test use it.
  - Headline offsets move: void -0.150 -> **-0.171 +/- 0.047**; group +0.166 -> **+0.200 +/-
    0.043**. Env-restricted: -0.088 +/- 0.038 and +0.123 +/- 0.034. The wall is 35.2% of the
    sample.
  - `test_comparison_bin_size_is_committed` had asserted the wall knee matched the global one to
    0.05 dex. That was true only because the box padded the wall with unclassifiable galaxies, so
    the test locked the defect in. It now asserts pool consistency.
- **#2.** The shuffle is stratified by (z bin x rms tercile; edges 0.81, 1.19 mJy); z-only is
  kept for comparison. Depth conditioning lowers every excess by 0.01-0.02 dex, as the referee
  predicted. A unit test now shows a pure depth confound fooling a z-only shuffle and not a
  (z, depth) one, but only when the Vmax is depth-blind: a Vmax that knows local depth already
  removes the confound. That is what the null can and cannot see.
- **#3 confirmed and tested.**
  - 99.4% of flagged group members' window galaxies are in the same Tempel group
    (`confused_same_group`), so the flag marks compact group cores.
  - Line-width test (`linewidth_residual`, log W50 on log M + z fitted to unflagged sources):
    flagged group members are broader by **0.019 +/- 0.004 dex** at fixed mass and z. The excess
    **grows with search radius**: 0.009 / 0.019 / 0.031 at 1.5 / 2.9 / 4.5 arcmin. Beam
    blending predicts the opposite (a neighbour outside the beam cannot broaden the profile).
    HI deficiency in dense cores (broader lines at fixed HI mass) predicts this direction.
  - So the round-6 draft's "carried by blending candidates" leaned the wrong way. The paper now
    says the group excess sits in compact cores and the data cannot separate blending from the
    core environment.
- **#4.** Excesses are quoted with sigmas throughout. With the restricted weighting, void 1.1
  (placements) to 2.4 sigma (shuffles); group 2.2 to 3.2 sigma. The void is weaker than in
  round 6 (the depth and footprint fixes each took some); **the group is now the more significant
  offset.**
- **#5.** Real-member dilution is stated for placements only: 33% in the placed voids against
  30% of the sample.
- **#6.** Clean-sample fit errors are used. The group without flagged sources is +0.038 +/-
  0.035 (1.4 sigma over the shuffles); the void without them is -0.078 (2.3 sigma).
- **#7.** The occupancy sign and "weakly" are restored: correlation -0.26 restricted, -0.45
  survey-wide. The placement end is the conservative one.
- **#8-12.**
  - "About half" is dropped: the text gives the excesses directly.
  - Void blending is no longer claimed.
  - The figure has measured fit-error bars and the legend below the panels.
  - A depth-confound test is added.
  - The title is now "...Against Two Null Tests"; percent precision is unified; "a typical"
    randomly placed volume.
  - Knee offsets print at three decimals.
- About 900 prose words plus the caption; 133-word abstract; one figure; triage and prose lint
  clean.

## Seventh referee round (2026-09-27): MAJOR revision

The round-6 fixes are all present (a single 55,893 pool, clean-sample errors, bracket sigmas,
placement-only dilution, occupancy sign, figure error bars). Every macro matches the JSON, the
sigmas were recomputed by hand, the figure matches, and the headline claims are pitched correctly.
The major findings are in the reasoning added in round 6:

1. **(major) "Grows with search radius, which beam blending does not predict" is not
   supported.** It refutes only a fixed-kernel, main-lobe blending model. Alternatives that
   can produce a rising trend:
   (a) DR2 integrates flux over an `ell_maj` aperture (3.5-5.2' in the first rows), not the
       beam, so neighbours at 2.9-4.5' can be inside the measurement;
   (b) the |dv| <= W50/2+100 window selects broader sources at fixed mass. Estimated
       0.012-0.046 dex for residual scatter 0.10-0.20, the size of the effect, plausibly
       growing with radius. The width test never used the fixed +/-300 flags;
   (c) the 4.5' flag is a density proxy (faint companions Tempel misses);
   (d) the samples are nested and there are no errors on the differences.
   Voids show the same trend (0.009/0.015/0.023) and flagged walls are broader (0.031, 7.4
   sigma) than flagged group members, so the effect is not specific to dense cores. **Tests:**
   the radius series with the fixed window; an annulus class (4.5' & ~2.9'); an aperture-aware
   flag (beyond ell_maj/2 + 1.45'); median log M and z per radius set.
2. **(major) HI deficiency has the wrong sign for the offset.** It lowers HI mass, so it
   predicts a LOWER group knee; the measured offset is +0.123/+0.200. It can explain the widths
   but not the sign. The environmental reading that fits the sign is segregation of massive
   galaxies into cores; otherwise only blending explains a raised knee. **My round-6 sentence
   was wrong in sign.**
3. **(major) "Depth" (`rms`) is aperture-integrated spectral noise, not depth.**
   rms/rms_beam = 1.47 x (ell_maj/2.9), so at fixed z it tracks source size and mass, and
   blending. Stratifying on it partly conditions on the outcome. `rms_beam` ("per-source
   detection sensitivity") is the depth. The terciles are coarse, and the unit test cannot
   catch this. **Fix:** stratify on rms_beam, with finer bins; say whether C/Vmax already
   include sensitivity.
4. **(major) The group placement null holds a tenth of the members** (1,909 vs 18,683, z 0.0225
   vs 0.031). Its spread (0.056) is 4x the shuffle's, so the 2.2 sigma end is mostly noise, and
   the "conservative end" argument holds for voids only. Report the occupancy and median z; say
   the group placement null is weak by construction.
5. The exclusion is a mass truncation (flagged median log M 9.92 vs 9.47); restore that caveat.
6. "Compact group cores" is untested; `confused_same_group` checks any two window galaxies, not
   counterpart plus neighbour. Use `clustercentric_radius` or write "close pairs within groups".
7. Residual model: compare in matched (M, z) bins; refit on the 4.5'-unflagged set.
8. Reduce the recommendation to what was shown: two nulls with occupancy, plus a confusion test
   with a fixed, outcome-independent window.
9. Footprint edge cells admit unsearched sky; check 0.5-degree cells or a minimum N per cell.
10-12 (nits): the 2.9' FWHM used as a radius reaches the 6% response point; "one is normally
    the counterpart"; the results paragraph chains 16 numbers, so move the survey-wide ones
    to the JSON.

### Response (run `fashienv-round7`: the robustness leg rebuilt; headline JSON unchanged)

- **#1 settled, and the referee was right.** Line widths use a fixed +/-300 km/s window,
  compared within (log M, z) cells (`matched_median_diff`) against isolated sources. Classes:
  inner (neighbour inside ell_maj/2 + 1.45'), ring (neighbour only in a 3' ring beyond it;
  cannot be blended), isolated.
  - Group members: inner **-0.006 +/- 0.006 dex**, ring +0.010 +/- 0.005.
  - Fixed-window radius series: -0.008 / -0.005 / +0.001 at 1.5 / 2.9 / 4.5'. Flat.
  - **The round-6 0.019 dex broadening and its rise with radius were produced by the
    W50-dependent window** (finding 1b), not by anything physical. The "blending does not
    predict" sentence is retracted.
  - Equal widths at fixed mass also cannot exclude blending, because a blend's mass and width
    grow together. The paper says so.
- **#2 accepted; the sign was wrong.** The HI-deficiency sentence is gone. The Discussion names
  the two readings that fit a POSITIVE knee: blending, and massive galaxies concentrated towards
  group centres.
- **#3 fixed.** The shuffle is stratified on `rms_beam` (per-source detection sensitivity)
  quintiles, edges 0.40 / 0.47 / 0.55 / 0.66 mJy; terciles give the same result (converged).
  - Restricted weighting: void -0.103 (2.5 sigma), group +0.128 (3.5 sigma).
  - z-only: -0.110 / +0.139.
  - With the right column, depth conditioning barely moves anything. The earlier `rms` strata
    had moved the nulls because they tracked source size.
  - New unit test: a stratifier that tracks mass absorbs a real offset; an independent one does
    not.
- **#4.** The group placement null's occupancy is stated (1,909 at z 0.022 vs 18,683 at 0.031),
  "weak by construction", and the conservative-end argument is limited to voids.
- **#5.** The exclusion is stated as a mass truncation: "a drop that removing the most massive
  members would cause whatever its origin".
- **#6.** Flagged group members sit at median R/R200 **0.35 vs 0.52**, so they are nearer group
  centres.
- **#7.** Cell-matched comparisons are the headline; the residual-fit versions are in the JSON.
- **#8.** The recommendation is narrowed: two nulls with occupancy, plus confusion flags built
  with a line-width-independent window.
- **#9.** A strict footprint (at least 5 Tempel galaxies per cell; pool 55,651) moves no offset
  by more than 0.001 dex.
- **#10-12.** "FWHM used as a radius"; "normally one of them its counterpart"; the survey-wide
  numbers are reduced to one clause.

## Eighth referee round (2026-09-27): MINOR revision

"The paper claims neither offset as environmental, and its evidence supports that." All 38
macros match the JSON, every sigma was recomputed by hand, the figure matches, and four DOIs
check out on Crossref. What remains is in the group half, and all of it changes what the
abstract should say:

1. **The group placement sigma should leave the abstract.** Its null is occupancy-biased, not
   just noisy: the restricted-weighting offset rises 0.021 dex per 100 members (corr +0.45),
   and the null mean moves -0.098 -> +0.054 across its own occupancy range. At the real 18,683
   members an occupancy-matched null would likely sit closer to the measurement. Describe it as
   uninformative, with the trend.
2. **The line-width classes overlap, and there is no counterpart matching.** "ring" and
   "isolated" share 449 group members (1,110 field). A source with no Tempel counterpart and
   one galaxy in its aperture is called isolated although it is exactly a blending candidate.
   Match the counterpart explicitly; give no-counterpart sources their own class; make the
   classes disjoint.
3. **Residual mass mismatch inside the 0.25-dex cells** (0.02-0.04 dex x slope 0.41) is as
   large as the quoted +/-0.006. It biases inner toward broader, so the corrected difference is
   more negative, which is what a small-offset blend predicts (summed mass, individual widths).
   Report the within-cell dM and dz; try 0.1-dex cells.
4. **The group interpretation is incomplete.**
   (a) R/R200 0.35 vs 0.52 is a selection effect under both readings; drop "nearer centres"
       from the abstract or say it follows from the selection.
   (b) The data support discriminators before interferometry: HI mass at fixed optical
       luminosity (needs Tempel r-band); the inner class split by neighbour velocity offset,
       plus W20/W50; knee vs R/R200 for ring + isolated members only.
   (c) Missing candidates: optically dark HI (tidal or intragroup gas); fit shape (group alpha
       -1.13 vs field -1.32, chi2nu 15/13), so show the offset with a common alpha.
   (d) Check whether Jones+2018 found the knee RISING with density. If so it is the precedent
       for a positive offset and the wrong citation for stripping.
5. **The headline exclusion uses the W50-dependent window the paper recommends against.** Use
   fixed300 throughout (29.1%, +0.044 +/- 0.035, 2.0 sigma). Carry the mass-truncation caveat
   into the abstract.
6. **Say why depth conditioning barely matters.** C and Vmax already track rms_beam at fixed z
   and M (corr -0.4 to -0.55), and the weights carry most of the depth dependence. rms_beam is
   an acceptable proxy.
7. **The strict-footprint check can hardly fail.** It removes 0.4% of the pool. Use 50% of the
   median cell count, or 0.5-degree cells.
8-10 (nits):
   - The void placement sentence mixes the 1,000-placement and 200-placement nulls.
   - "For 99.4% the neighbour is in the same group" overstates the any-two test.
   - Void RR macros render as "--" (unused).
   - The abstract chains four sigmas; results paragraph 1 chains about 12 numbers.
