# Changelog

All notable changes to `jansky-research` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/) as codified in [`VERSIONING.md`](VERSIONING.md).

Every PR adds an entry to `## [Unreleased]`. `scripts/next_version.py` reads that section to
recommend the next version number.


## [Unreleased]

### Changed
- **Review round 6 verdicts recorded: `peaked`, `stokesv`, `swaves`, `windwaves` -- four major
  revisions, eight blockers.** The batch was chosen by inherited-defect risk and every
  inheritance hypothesis paid out. peaked (3 blockers): the "100% rising" validation is a flux
  cut that cannot fail (measured against the sample's real TGSS fluxes: ~53%, a third FALLS);
  the candidate list exists in no committed file and its vetting step in no code; the global
  25 mJy TGSS limit vs the field's own local noise moves the count 6 -> 1-3. stokesv (3
  blockers): two of nine "circular detections" have |V|/I of 568% and 135% (the I and V cutouts
  are paired by independent unordered CASDA queries with no same-observation constraint); the
  fifteen rows exist only inside a PDF. swaves (1): the headline speed has no uncertainty and
  came from the live non-converged fit (slope/mask/R^2 from different iterations); the
  "fundamental" number is hand-halved; the title distance is the band edge. windwaves (1): the
  self-declared dominant systematic (CME density) is never propagated and both abstract
  framings flip inside its stated range (0.08-0.17 c, 10-20 R_sun; harmonic and density are
  exactly degenerate). Verdicts and full finding lists in each survey/<slice>-findings.md;
  fixes to follow per slice.

### Changed
- **Review round 5 verdicts recorded (presenter + referee): `southern`, `hi`, `torchfdmt`,
  `solarbursts` -- four major revisions, five blockers.** southern: GLEAM-X DR2 *does* publish
  per-band errors (the assumed 10% proportional error is ~10x too small for the faint majority)
  and the alpha_lo purity gate is anchored on the noisiest sub-band -- three of the figure's own
  six candidates flip under a one-sub-band change; hi: the load-bearing ~7 km/s citation is a
  CO-vs-HI cross-tracer offset, not the threshold-vs-fit measurement the paper says it is, and
  "flat" is never quantified; torchfdmt: the committed figure plots the raw track sum (maximum
  at DM ~118, not the quoted 56.59) and the butterfly S/N 6.0 carries no trials factor against
  an expected noise maximum of ~5.2, while the 29x benchmark ratio splices two invocations
  keeping the CPU number that flatters it; solarbursts: the figure plots the unclipped fit
  (R^2 0.41) under the R^2 0.811 caption, and the robust fit is unconverged -- the headline
  speed spans 0.111-0.147 c over a hard-coded iteration count. Verdicts and full finding lists
  in each survey/<slice>-findings.md; fixes to follow per slice.

### Fixed
- **`peaked` revision (round-6 referee, three blockers): the validation became real, the limit
  became local, and the candidate list became evidence.** `validate_hfp` now uses each Dallacasa
  source's measured TGSS flux (83/100 detected) with the VLASS QL correction: 53% rising / 32%
  falling / 38% still rising at 3 GHz -- the old 100% was a flux cut the sample could not fail,
  and the paper says so. The TGSS limit is the field's own 7-sigma floor (29.4 mJy from the ADR1
  local-rms column) with the count sweep committed (8/6/1/0 at 20/25/29.4/39.9 mJy); the
  headline is 1 peaked + 1 GHz-peaked (the latter appearing exactly as predicted once the VLASS
  flux-scale correction was applied), and the old six-candidate list is described as a statement
  about the assumed limit. results/peaked_candidates.csv publishes every rising candidate with
  index errors and in-pipeline SIMBAD/NED vetting; all three Callingham bins are committed
  (0/81, 3/26, 1/6); macros are namespaced with '--' validation defaults (closing the
  zero-overwrite hazard); the fixture gains extended fakes; the hand-authored arxiv.yaml is
  rewritten from the committed metrics.

- **`windwaves` + `swaves` revision (round-6 referees, one blocker each): the estimator is
  rebuilt on the independent unit, and the dominant systematic is propagated.** The shared
  `beam_speed` now fits ONE point per distinct time sample (column median) with a
  leave-one-sample-out jackknife -- en route, the converged per-point clip was measured to
  COLLAPSE on the STEREO ridge (2 surviving samples, 0.031 c vs 0.160 c), the leverage
  pathology both referees circled, seen live. windwaves: 0.102 +/- 0.011 c with the committed
  (mode x density) grid spanning 0.054-0.215 c and 6-19.5 R_sun, the f_p^2 ~ n degeneracy
  stated, and the abstract claiming neither distance regime exclusively; the residuals'
  band-monotonic structure is reported as unexplained (the matched-cadence fixture disproves
  the old band-compression story). swaves: 0.160 +/- 0.014 c at R^2 = 0.986 over 11 samples
  with a tight estimator bracket -- the track is genuinely straight at sample level -- the
  hand-halved fundamental replaced by grid macros (0.079 c / 41 R_sun), and the title carries
  the 0.2-0.4 AU range. Both slices: ridge CSVs, burst epochs, pad/snr provenance, namespaced
  macros, figures that draw the exact quoted fit, matched-cadence fixtures, and
  krupar2015/reiner2015 fixed against Crossref.

- **`southern` revision (round-5 referee, 18 findings, two blockers): re-run with the
  catalogue's own errors, the method survives and the over-density vanishes.** GLEAM-X DR2
  *does* publish per-sub-band uncertainties; `fetch_gleamx` now uses them, fits admit only
  >=3-sigma detections (closing the positivity-censoring bias), and the rising-side gate is a
  fitted in-band index required significantly rising/flat instead of a two-point ratio anchored
  on the noisiest sub-band. The committed cascade runs 160 naive -> 49 -> 41 candidates
  (was 90), median nu_pk 253 MHz, median reduced chi^2 1.45; the density ratio vs RadioSED II
  is 1.21 (the old ~2.5x over-selection was the assumed-error artifact, and the paper says so).
  Callingham validation with honest denominators: 31/50 covered recovered at 0.091 dex (90%
  within 2x; unconditional 0.095 dex, 76% of all tested), per-bin tested/recovered committed so
  "recovery climbs with nu_pk" is rate-backed. The full 1,545-source catalogue is committed
  (positions, classes, nu_pk, alpha +/- err, chi^2), crossmatch multiplicity measured (0), the
  gap fraction (0.49 interpolated turnovers) stated in the abstract, the synthetic field can
  fail on both cuts (noise floors + flattening contaminant + extended fakes), band counts and
  racsmid2024/kerrison2025 citations fixed, and the arXiv package rebuilt clean.

- **`solarbursts` revision (round-5 referee, 18 findings, one blocker): the fit converges, the
  figure draws the fit it captions, and the event selection became falsifiable.**
  `_robust_linfit` gains a converge mode (mask fixed point; slope, mask and R^2 from the same
  final fit; the legacy 3-iteration default is kept for the three other slices whose committed
  evidence used it). The converged headline is 0.1173 c (R^2 0.935, 55/66 channels, fit band
  32.4-62.4 MHz) with a committed analysis-choice spread 0.1173-0.1332 c, replacing a point
  value whose undeclared iteration-count ambiguity spanned 0.111-0.147. `run_candidates` commits
  all four candidates plus the X6.9 storm control at the headline parameterization (JSON +
  macros): the rejects give R^2 0.48/0.054/0.342 -- one with a wrong-direction positive
  drift -- and no longer all land in the canonical band, so the recover-a-known can fail.
  Band/drift bound to the used channel set; provenance (pad_s, snr, clip sigma) in the results
  and the ridge CSV; six-figure km/s rounded; harmonic assumption stated with the F-H check;
  "directly consistent" downgraded to a range comparison; stale findings sections marked
  superseded; arXiv package rebuilt.

- **`torchfdmt` revision (round-5 referee, 15 findings, two blockers): the figure now plots the
  statistic that produced its own caption's number, and the significance carries a trials
  accounting.** The butterfly panel plotted the raw track sum (maximum at DM ~118); it now
  plots the per-row z-scored peak that `best()` maximises, with a test pinning figure == quoted
  value. A 200-rep per-channel circular-shift null is committed (median 5.23, p99 6.03, max
  6.26 -- the referee's Gumbel estimate to the decimal): the observed 6.0 has p = 0.030
  against it, so the paper now leads with the positional coincidence (2.9 of 1,903 DM trials,
  p = 0.0035) and the boxcar S/N 14. The benchmark was re-measured in a single invocation with
  `--bench-devices` decoupled from the science device and the hardware string emitted by torch
  introspection: brute speedup is 24x (the spliced row's 29x had retained the older, slower CPU
  time). En route, a guard bug: `preserve_live_macros` treated torchfdmt's mixed source string
  as synthetic and silently discarded real reruns' macro updates; its rule now matches
  `preserve_live_results` (mixed counts as real), with tests. Sclocco et al. 2016 cited and the
  universal CUDA claim narrowed; boundary condition, oracle timing, and pulse position/width
  committed; stale findings table corrected in place; arXiv package rebuilt.

- **`hi` revision (round-5 referee, 13 findings): the citation became a measurement.** The paper
  attributed its 9% normalization excess to a ~7 km/s threshold-vs-fit offset cited from
  McClure-Griffiths & Dickey 2016 -- whose number is actually a Clemens-1985-CO-vs-HI
  cross-tracer comparison. `terminal_velocity_edge` (error-function edge fit) now runs on the
  same spectra: the measured threshold-vs-edge offset is 13.8 km/s and the edge-fit flat mean is
  242.9 +/- 8.0 km/s, within 2.9% of the Reid 2019 V0 -- the excess is demonstrated in-house to
  be the estimator. Flatness is quantified (slope 2.8 +/- 1.8 km/s/kpc; Keplerian contrast 176
  vs 263 at the outermost point); the synthetic fixture's edge widened to a realistic 5 km/s so
  the wing overshoot is reproduced offline (13.6 km/s) instead of hidden; threshold sweep and
  drop-a-point sensitivities committed; the +/-6 relabelled as scatter with the SEM quoted; a
  pipeline-written per-longitude table added; macros namespaced (hiSyn*/hiReal*); figure shows
  both estimators, the Keplerian curve, and the excluded bar points; arXiv package rebuilt.

- **`rmstructure` revision (referee round of 2026-08-24, two blockers): the sample now matches
  the release and the headline error was understated eleven-fold.** The "S/N >= 8" sample
  (333,173 rows) was the goodRM sample with tile-overlap duplicates; `load_spice_racs_dr2` now
  dedups to one row per `cat_id` (246,508 -- the release's published post-dedup count) and
  commits the full cut cascade with the S/N column named. The enhancement ratio's error is now
  a leave-one-block-out jackknife over 601 10-degree sky blocks (11.0 +/- 1.1, display
  precision matched) instead of the i.i.d. source bootstrap (+/-0.10) the slice's own recorded
  lesson warns about, with a test pinning jackknife > bootstrap on the correlated fixture. The
  quality-flag claim is measured on committed runs of both variants (unflagged high-|b| break
  0.87 deg vs 2.29 deg flagged; the old 0.5-vs-3.70 numbers were pre-dedup and uncommitted).
  The 30-seed ensemble ratios are committed (`results/rmstructure_synthetic.json`); ladder rows
  carry plateau/sigma errors and pair accounting; the floor subtraction quotes its one input
  (12.3 +/- 0.21) and its latitude-dependent leverage; provenance macros are namespaced so an
  offline rebuild cannot invert the file's marker; the fixture's scope (pole cut 15 vs 60 deg,
  band-median vs peak boost) is stated; thomson2023 is cited; the stale arXiv package with the
  dangling-citation abstract is rebuilt clean under the 1,920-char limit.
- **`vlass` revision (referee round of 2026-08-24, one blocker): the completeness curve
  measured its own error model, and the archival rejection is now committed evidence.**
  `injection_recovery` scaled the whole per-epoch error with the injected flare
  (`e[k] *= fac`), capping the flare epoch's chi^2 and manufacturing an apparent 52%
  saturation attributed to a "three-epoch ceiling" the committed thresholds refute; it now
  rescales only the flux-proportional systematic term, injects from the census's own usable
  population (including two-epoch light curves) against the census's committed thresholds,
  extends the factor grid to 50x, and repeats over ten seeds so the 50% factor carries a
  Monte-Carlo error. The prose-only archival rejection of the second image-confirmed
  candidate is replaced by a coded stage (`archival_vet`/`run_archival`: Epoch-4 forced
  photometry + NVSS/FIRST/TGSS/AllWISE counterparts, committed to
  `results/vlass_archival.json` and an `archival_survives` CSV column); per-epoch forced
  rms/offsets are committed so every confirmation and rejection is auditable; the census
  count closure (3,722 Epoch-1-only sources) and Epoch-1 anchoring are stated; the stale
  "few x 10^-4" fraction is replaced by pipeline-derived floor/corrected fractions; and the
  paper drops "necessary and sufficient", fixes the gordon2021 title (right DOI, wrong
  title), and quotes the candidates' catalogue-vs-forced amplitude collapse from committed
  numbers.
- **`glitchpop` revision (all 19 referee findings): the title flip was Monte-Carlo noise, and
  the census is pinned.** At 2x10^5 bootstrap draws with per-pulsar seeds and valid
  (k+1)/(B+1) p-values, the headline flip (J2229+6114) dissolves along with a second marginal
  quasi-periodic call --- n_flipped = 0, QP count 10 -> 8, while the aggregate excess stands
  (Poisson-binomial p = 7.9e-5 under empirically measured count-dependent false rates). The
  title becomes "the Post-2018 Increment" and the paper's finding is that no classification
  change is secure at current sample sizes. The catalogue snapshot is committed (728 rows;
  a fourth distinct count in four retrievals) with a from_csv analysis path and stated
  retrieval date. "~184 post-2018 glitches" was vintage arithmetic: the measured increment is
  **89**, with 21 retroactive pre-2019 additions counted via the previously discarded is_new
  flag. The validation can now fail (injection surface: completeness 0.16-0.44 in the
  borderline cv/n regime the old cv-0.12-only test never visited); the pooled false-positive
  rate is disaggregated by outcome; the clustered "lower bound" is withdrawn to
  chance-consistent (expected 1.55 false, p = 0.80); the hidden second cut is stated; the
  gap-factor sweep is committed; the out-of-sample Howitt anchor (B1338-62) is promoted into
  the paper; and a real calibration anomaly is surfaced and stated (census p-values non-uniform,
  KS p = 1e-4 --- mild population-wide regularity or a missing monitoring dead-time in the null).

### Fixed
- **`typeii` revision (all 15 referee findings): the contaminant that can fail was injected,
  and it failed.** A 30-seed ensemble measures the detector per contaminant class: clean against
  fast type III (0/480) and narrowband RFI (0/240), but **false-positive at ~38%** (91/240)
  against coherent slow-drift window-filling background --- the class the real band contains and
  the old single-seed purity 1.0 never tested. "This is not a detector defect" is retracted; the
  census's 83% window saturation and the ensemble's failure mode now corroborate each other, and
  the null is stronger for it. **The match-rate deficit is retracted as a coverage artifact**:
  with coverage-aware diagnostics and a fresh CDAW fetch, the covered match rate is 0.68, above
  the 0.612 chance rate; the committed `association_is_background_like` no longer keys on it.
  The speed test is sharpened with the matched fast-fraction (0.184 vs background 0.057, KS
  p ~ 0): a real but modest enrichment, the flare--CME confound's signature. Also: the
  **arxiv.yaml blocker fixed** ("48 days, zero failures" --- the original clobber's number
  surviving in the hand-authored submission abstract --- is now 768); the synthetic leg writes
  its own committed results file (the paper's provenance claim had been false for every
  synthetic number); the harmonic-cut sweep committed with denominators ("confirming" downgraded
  --- 8 matched CMEs, trend reversing at the next cut); 332 detections = 320 distinct structures
  after adjacent-window dedup, stated; the RSTN comparison scoped as cycle-phase inflated; the
  committed figure's inverted frequency axis and null-as-zero bars fixed.

### Added
- **Review round 4: full presenter/referee rounds on rmstructure, typeii, glitchpop and vlass.**
  All four returned **major revision** (six blockers between them); every verdict and complete
  finding list is in the slice's findings file, fixes pending:
  - *rmstructure*: the headline 11.17 +/- 0.10 uses the i.i.d. bootstrap the slice's own
    recorded lesson warns about, on a spatially correlated, ~25%-duplicated sample (honest error
    likely ~6x larger; fix is a spatial block jackknife); the sample definition conflicts with
    the release's published counts; the abstract's closing claim rests on an uncommitted run.
  - *typeii*: the tracked arxiv.yaml still says "48 days, zero failures" --- the original
    clobber's number, in the submission metadata, outside the macro system's protection; the
    CME match-rate deficit REVERSES (0.676 > chance 0.619) once the 57 detections beyond
    CDAW's coverage end are excluded, so the committed association_is_background_like boolean
    is wrong as computed; the synthetic leg has no committed results file at all.
  - *glitchpop*: the title flip has margin two bootstrap replicates (0.15 MC-SE of its own
    null), the abstract counts glitches and pulsars on different samples, no catalogue snapshot
    is committed (three scrapes, three counts), and the out-of-sample Howitt+2018 anchor the
    validation needs (B1338-62) exists in the census and is never used.
  - *vlass*: the 52% completeness saturation is an artifact of `e[k] *= fac` in the injection's
    error model (the paper's "three-epoch ceiling" explanation is refuted by its own committed
    v_threshold), and the second candidate's rejection --- the step that yields "one genuine
    variable" --- exists only in prose while the committed CSV still marks it confirmed.
- README corrections from the round: the glitchpop rows now match the committed census
  (31 classifiable: 20/10/1 --- the old row omitted the clustered class entirely), and the
  rmstructure row no longer says the DR2 contrast "awaits the public file".

### Fixed
- **`rmdipole` revision (all 16 referee findings): the null becomes a limit, and every
  diagnostic resolved in the paper's favour.** The scramble null's amplitude percentiles are
  now committed per leg (the old run stripped them), so the headline reads "RA-projected
  power-dipole amplitudes |p|/m > 0.354 excluded at 95%; nothing below that constrained"
  instead of an unscoped "isotropic". The clip cannot manufacture the result (an injected 0.3
  dipole survives the 0.99 clip at p=0.001; quantile swept 0.95--0.995, p stable at
  0.78--0.93); the injection is unbiased (0.3006 +/- 0.0091 over 20 realizations --- the old
  single-seed "partial-sky bias" explanation was wrong and the paper now says so); a dipole
  painted onto the real heavy-tailed field is detectable at p=0.001 (with an uncalibrated
  amplitude scale, so conclusions are stated as detection/exclusion); the clipped-out tail's
  own apex is (25 deg, +65 deg), at the survey's Dec edge, not the full-sample apex previously
  imputed to it; and the >=5-neighbour mask is measured inert. Bib: four now-published
  preprints completed (incl. the DR2 data citation) and three missing DOIs added; the
  kinematic expectation restated in the fitted statistic; apex range corrected to 80--153 deg.
  All original leg values unchanged.

### Fixed
- **`frbwait` revision (all 15 referee findings): the population claim is retracted, and the
  experiment that killed it is the paper's new spine.** A transit-censored injection curve
  (pure Poisson processes observed through CHIME's 15-minute daily window, fitted exactly as
  the census fits real sources) recovers k = 0.78--0.91 at typical census rates and **k = 0.47
  at 3.2 bursts/hr** --- so the census median (0.831, CI 0.73--0.98) is what Poisson arrivals
  look like through the transit comb, and "sub-Poissonian clustering is the population norm"
  is gone. The three flagged sources still sit below the curve at their own rates (some
  clustering beyond selection remains indicated); one source is significantly
  super-Poissonian and is no longer hidden by the one-sided flag. A group-preserving scramble
  (whole transits move rigidly) settles the detections: the anchor survives at the floor while
  the two other "significant" peaks collapse to p = 0.046 and 0.021 --- their significance was
  within-transit multiplicity. FRB20220912A's peak is flagged as railed at the period-grid
  edge. The anchor's duty is convention-matched (the full-containment arc, 6.3 d, is wider
  than the published 5-d window; the 90% arc, 3.5 d, sits inside it --- the old "consistent
  once conventions are matched" was wrong as written). Median-k CI and sign test, run
  configuration, per-row declinations, grid-edge flags and the bias curve are all committed;
  cat2's bib entry gains its published volume/page; the placeholder-asserting test is replaced
  by a two-step guard test. No previously published value moved.

### Fixed
- **`junodam` revision (all 16 referee findings): the null is quantified, and it contains
  structure.** The blocker --- a one-sample test that existed nowhere in committed evidence ---
  is closed by computing and committing it: two-sided sign test p = 0.453, geometric-mean
  monthly contrast 0.92, **95% CI 0.53--1.59**, so the paper now states that "does not reject
  unity" equally does not reject a 1.6x enhancement, and the title softens from "Do Not
  Coherently Organise" to "Are Not Detected to Organise". The new evidence surfaced two real
  things: **Io-B alone is enhanced (1.67** vs 0.92/0.79/1.28 for A/C/D), the per-region signal
  the union contrast averages away; and the box-shift scan is flat at every rigid Io-phase
  offset (0.75--1.18), so no convention error rescues the boxes and the machinery-validated
  scan (it peaks at zero on injections) makes the weak organisation a measurement. Also: the
  aggregate contrast now carries a day-block bootstrap error (**1.12 +/- 0.18**; the 44,294
  active bins are 6,784 episodes, so per-bin errors would overstate N by ~170x); the injection
  calibration extends to the near-unity regime (injected 1.25 recovers 1.16, the ~30%
  contraction of the estimator stated with its direction); the noise-promoting upward-rescale
  null is replaced as primary by a censored census that cannot promote noise (near/far 2.71,
  same-detector pair 330.4 -> 2.71 quoted in the abstract instead of mixed detectors); the
  single v01 month is measured, not assumed, benign (v02-only CI 0.48--1.33); the dataset
  citation is repointed to the v02 DOI with the correct author list; and the per-month table,
  per-quartile denominators and raw-p90 duty cycles are committed. One previously published
  value moved: near_far_corrected 2.2 -> 2.25 (the referee's rounding finding).

### Fixed
- **`skr` revision (all 21 referee findings): the 1.39 is retracted, and the retraction is the
  result.** Three defensible constructions of the 1/r^2 sensitivity null give three answers on
  the same data --- 1.39 (rescale total flux: moves the range-independent noise floor with the
  signal, biased toward a collapse), 1.68 (rescale excess upward: amplifies far-range noise,
  biased toward a reversal; the flat-rate control comes out at 0.26), and 3.58 (the adopted
  common-sensitivity census, which only ever scales excesses down and passes both injection
  controls). The adopted estimator's own numbers then show the decomposition is not measurable:
  orbit jackknife +/-2.75, rule sweep 1.14--34.95, per-day background 0.59. The paper now quotes
  the raw trend (**3.33 +/- 0.93**, leave-one-periapsis-out over 10 passes) and no residual, and
  says plainly that the old "the null removes essentially the whole effect" was an artifact of
  rescaling the floor. Also: the rotation anchor restated as a location match (0.05%/0.06%
  within a 0.08 h resolution; peak power NOT significant under a day-block permutation, p=0.325;
  the analytic ls_fap dropped); the ~10.34 h broad-band peak committed as evidence; the stale
  junodam-priority sentence deleted; provan2019's author list and pages corrected per Crossref;
  the fabricated ye2016 "Local Time and Latitude Dependence" entry (a subtitle and page range
  that exist nowhere) removed and the duplicate-DOI pair collapsed to the one real paper;
  gurnett2009's title corrected; the detection rule, band, reference range, bin edges (periapsis
  reaches ~1 R_S, not the "~8" the paper claimed), corrected per-quartile duty cycles, date
  range and NaN count all committed; the latitude confound's sign stated; the figure added to
  the paper; the DOY-bucket hardcode in the fetcher fixed. Two new controls that can fail are in
  the test suite, and both rejected null models measurably fail one of them.

### Changed
- **README rewritten in a plain register and brought current.** Stale facts fixed: "Forty
  slices" is forty-five; "What's next" recommended work that was executed weeks ago and now
  describes the actual current work (review campaign, publishing queue, station); the arXiv
  queue reads atlas3i -> dr20radio -> lptv -> innerrc; the three RNAAS notes are listed; the
  skills/agents/plans listings match the tree; the missing lptduty row is added. The Reviewed
  column now records the 2026-08 referee-driven fixes paper by paper. Register: the 22 longest
  table cells cut from paragraphs to a line, "honest/honestly" reduced from ~30 uses to one,
  bold-pair contrasts and stock phrasing removed.

### Added
- **Review round 3: full presenter/referee rounds on rmdipole, frbwait, skr and junodam** ---
  the four highest-risk unreviewed papers. All four returned **major revision** (three with a
  blocker); every verdict is recorded in the slice's findings file, fixes pending:
  - *rmdipole*: the isotropy null has no stated sensitivity, and the evidence needed to state
    one (the scramble-null amplitude distribution) is stripped before commit. The tail clip is
    never tested against a genuine dipole; the injection control is single-seed and Gaussian
    where the real field is heavy-tailed.
  - *frbwait*: the three "clustered" sources are exactly the three highest-rate sources
    (p = 0.0022 by chance), which is what the paper's own disclosed censoring bias predicts;
    two significantly super-Poissonian sources are never mentioned; the offline validation runs
    ~25x below the worst-case rate and cannot fail.
  - *skr*: the sibling-census comparison is contradicted by the repo's own committed junodam
    evidence; the headline 1.39x has no uncertainty; the 1/r^2 null rescales the noise floor
    along with the signal, in the direction that manufactures the collapse.
  - *junodam*: the abstract's "one-sample test against unity does not reject" exists nowhere in
    committed evidence, and the 95% CI (0.49-1.64) admits a 1.6x enhancement; the dataset DOI is
    the obsoleted v01 with a wrong author list; the frame convention whose failure mode IS the
    reported null has no test that could catch it.

### Fixed
- **junodam's committed macros were carrying a synthetic run's values** --- the abstract's opening
  counts rendered the 28-day fixture (161,280 bins / 12,729 active) instead of the real 210-day
  census (1,209,600 / 44,294), with `\jdSource` reading "synthetic orbit". Introduced by the
  2026-08-23 commit that fixed a different instance of the same clobber class; caught by the
  presenter round. Macros regenerated from the committed real JSON, and the four remaining
  mode-dependent names (`\jdNbins`, `\jdNact`, `\jdOccIo`, `\jdOccOut`) namespaced
  `jdReal*`/`jdSyn*`; the namespace audit is clean and the referee verified the fix, including
  that the committed figure is the real run.

### Added
- **`triangulate`: the miss-distance sweep is measured, and the cut is not load-bearing.** The real
  2013-05-15 leg was re-fetched from SPDF and the threshold swept over 15/30/60/100 R_sun via the
  new tested `miss_sweep` (pure filtering on a track built with the cut open). Tightening 60 -> 15
  keeps only 12 of 38 channels and moves the correlation 0.989 -> 0.977 and the ratio 2.18 -> 2.20;
  no kept channel misses by more than 60, so the 100-row is identical. Stated in the paper with
  `\triSweep*` macros, and the per-channel arrays (cut open) are committed as
  `results/triangulate_channels.csv`, so the choice is auditable and any future sweep is offline.
- **`rfitrend`: the validation can now fail.** A third fixture arm injects local RFI into the
  lines' *flanking* channels only --- the systematic the paper names as decisive and the
  line-vs-adjacent difference cannot cancel, which the old arms (common-mode gain, broadband
  bursts) cancelled algebraically. Measured: the same injected rising line (+0.2404/yr clean)
  recovers at **-0.1907/yr** under flank contamination, a bias of -0.4311/yr that **flips the
  recovered sign** --- the ALMATY mechanism reproduced end-to-end. The paper's validation paragraph
  now reports this (`\rfSynFlankSlope`/`\rfSynFlankBias`) instead of disclaiming it.
- **`solarbursts`: the systematics grid is emitted, and the stale-grid mystery is solved.**
  `speed_grid` computes all three harmonic x fold points from the same fitted ridge as the
  headline, so the middle point IS `\sbSpeedC` by construction; the raw 66-channel ridge is
  committed as `results/solarbursts_ridge.csv`. The discrepancy the findings file recorded was a
  *parameterization*, not code drift: `pad_s=10.0` reproduces the committed evidence byte-for-byte,
  while `--recover` had `pad_s=5.0` pinned --- so the one command documented as regenerating the
  committed result produced different numbers (r2 0.897 vs 0.811, speed 0.1368 vs 0.1347).
  `--recover` now pins the committed parameterization. Grid: 0.0813 / 0.1347 / 0.2448 c.
- **`ecallisto_census`: the committed real leg is cited.** `\ecsReal*` macros (both namespaces on
  every run, per the accumulates-values-not-names rule) carry the 168-day illustrative ingest
  (5 events on 2 days, r = 0.28), and the limitations section now presents it as the data-volume
  point made concrete rather than leaving committed real evidence unmentioned.

### Fixed
- Two hand-typed numbers in `solarbursts` were stale against the committed grid and are now
  macro-backed: the abstract's bracket ("0.09 c ... 0.27 c" -> `\sbGridFundOne`/`\sbGridHarmFour`)
  and the Discussion's outer radius ("3.6 R_sun" -> `\sbGridHarmFourRhi` = 4.01, which
  extrapolates *beyond* the Newkirk model's constrained range, not "near its edge"). And the claim
  "All three lie within the canonical 0.1--0.5 c band" was false --- the fundamental point (0.0813,
  previously 0.086) is below it --- and now says so.

All three real re-runs (SPDF, e-Callisto) are purely additive against the committed evidence:
no previously published value moved.

### Fixed
- **`pulsarspec`'s null now has a sensitivity, which turns it from an impression into a limit.**
  "Millisecond pulsars are not significantly flatter" rested on two rounded means; `run()` computed
  both subsample dispersions and discarded them, so nothing in the committed evidence quantified
  "indistinguishable". New tested helper `compare_subsamples` returns the difference, its Welch
  standard error, the observed significance and **`resolvable`** --- the smallest offset the two
  subsamples could have distinguished from zero at 2 sigma. Measured on the real catalogue: the
  means differ by **0.018 +/- 0.122 (0.15 sigma)**, and the sample resolves an offset only above
  **0.24 in alpha**. The abstract and Results now state that bound: offsets above 0.24 are excluded,
  smaller ones are not constrained. This is the `frblens` lesson applied to a different shape of
  null --- a non-detection divides by what the measurement could have seen.
- The same slice's offline fixture injects a millisecond--normal offset of **+0.2**, which is
  *below* the real sample's 0.24 resolution, and **no test ever asserted it was recovered** while
  the Methods section told readers the flatter-millisecond sub-population was validated.
  `test_offline_recovers_the_injected_msp_offset` now asserts it, at the fixture's size where it
  is detectable. The pre-existing `test_run_offline` window was ~28 sigma wide on an exact
  algebraic transform and could not fail.
- `pulsarspec` also gains a denominator (473 of **2536**, i.e. 18.7%, not "the whole catalogue") and
  states the *direction* of its selection bias: requiring detection in both bands biases a
  two-frequency index **flat**, since a steep source near the S400 limit has already dropped below
  the S1400 limit. Same mechanism recorded for `dr20radio`.
- **`ppdot`'s Crab validation now touches the real data path.** `fetch_atnf_ppdot` requested the
  `PSRJ` column and then discarded it, so the run could not identify any named pulsar and the
  abstract's two Crab numbers were typed from a textbook, present in no committed artifact. The
  column is kept, and the new tested `named_pulsar_derived` reads the Crab out of the analysed
  sample: **B = 3.8x10^12 G, tau = 1257 yr, P = 0.0333924 s**. Reading the row exercises the `P1`
  column units, the row parsing and the positive-Pdot cut; evaluating closed forms at hand-entered
  literature values, as the old unit test did, could fail only on a typo in two-term algebra.
- `ppdot`'s "Of the ~3500 entries" was a recollection and wrong: the table has **2536** rows, the
  same figure `pulsarspec` obtains from the identical call, so the two slices no longer disagree
  about the size of one table. Its hand-typed arithmetic ("a factor of ~80,000", "~1.8% below our
  line") is now derived in `run()` and macro-backed, so it cannot go stale after a re-run, and the
  false universal "Every number above is written by the pipeline" is scoped to catalogue-derived
  numbers, naming the two literature figures that are not.

Both real re-runs were purely additive: every previously published value is unchanged and no macro
changed value.

- **A test in the suite was making a live VizieR call**, and it failed CI here by timing out.
  `stokesv.forced_photometry_recover` fetched the radio-star catalogue *before* checking for
  `CASDA_USERNAME`, so `test_real_path_does_not_write_the_synthetic_figure_before_the_real_leg` ---
  which asserts the real path raises on a missing credential --- only reached its `RuntimeError`
  after a network round trip. It passed whenever VizieR was fast, which is why it had never been
  noticed. The credential is now checked first, which is also the right behaviour (a missing
  credential dooms the run either way, so there is no reason to query the archive). That test file
  now runs in under two seconds; it was the only test exercising a real path without a monkeypatch.

### Changed
- **Traditional-style conversion complete: all 45 papers.** Batch 8 converts `ppdot`, `pulsarspec`,
  `rfitrend`, `stacking` and `vlbi`. Every batch's referee round returned revisions (eight for
  eight), and all of this batch's findings are fixed here. The restyle MAJORs were again scope and
  claim-strength drift no linter sees: `stacking` moved its limitations sentence from *this work's*
  limitations to the method class while the list still opened "the result depends on ..."; `vlbi`
  softened "**We stress** that $f=0.05$ is only a common starting assumption" to "We note", the
  paper's only advance warning that the assumption its central result overturns is an assumption;
  `rfitrend` promoted a signpost to "the **decisive** check" twelve lines above the paragraph calling
  the same test "an underpowered global-vs-local discriminator"; `ppdot` collapsed an appositive so
  the sqrt(2) convention caveat modified a phantom list item rather than the field formula it is the
  only caveat on.
- **The em-dash-to-comma appositive collapse is now confirmed as this campaign's characteristic
  defect**, appearing in five papers across two batches. A comma re-attaches an appositive to the
  nearest noun: in `stacking` the ~2x excess acquired a generic plural subject instead of staying
  tethered to two quoted bin values; in `pulsarspec` "squarely in the literature range" moved inside
  the parenthetical series and so appeared to endorse the *scatter*, which the Discussion explicitly
  declines to endorse.

### Fixed
- `vlbi` called OJ 287 "the most variable source" by ranking on eta --- the statistic the paper's own
  thesis declares unusable as a discriminant. On V, the statistic it argues for, BL Lac leads (0.540)
  and OJ 287 is second (0.480), so "the other archetypal blazars (BL Lac, ...) follow" was false on
  the preferred measure. Also states that one of the 13 sources above the floor sits *at* it to the
  quoted precision, so the count is twelve unambiguous plus one tie.
- `rfitrend`'s abstract attributed the trend to the Starlink count three sentences before claiming no
  Starlink attribution; `\rfRealNMonths` = 286 was attached to three stations but counts two (286 =
  161 + 125; three would be 448); and the config-stability screen's attrition was invisible although
  `n_months_raw` was in the committed metrics all along. It is uneven and it matters: HUMAIN keeps
  161 of 174 and GLASGOW 162 of 163, but **ALMATY keeps 125 of 174 (28% dropped)**, and ALMATY is the
  station whose falling slope carries the sign disagreement the null rests on. Now emitted as
  `\rfReal<ST>NMonthsRaw` and stated.
- `stacking`'s magnitude paragraph pointed at the wrong figure panel (the caption puts magnitude in
  the middle; both paragraphs said "right"), and its "the monotonic ordering across the three bins is
  clean" is now stated at its measured strength: adjacent steps are 1.00 and 0.60 sigma, end-to-end
  1.60 sigma, so the ordering is monotonic in the central values only.

### Documented
- Two pre-existing MAJORs are recorded in the findings files rather than patched, because each needs
  a real re-run (VizieR confirmed reachable, so both are tractable next):
  **`pulsarspec`'s null has no sensitivity.** "Millisecond pulsars are not significantly flatter"
  rests on two rounded means with no dispersion, no N for the normal sample, no SE and no test
  statistic in the metrics --- the code computes the dispersions and discards them. Reconstructed, the
  difference is 0.02 +/- 0.12 (**0.17 sigma**) and the smallest resolvable offset is **~0.24 in
  alpha**, which is the paper's real result and is nowhere stated. The offline fixture injects
  **+0.2**, *smaller than the real sample's 2-sigma threshold*: the analysis could not have detected
  the offset its own validation builds in.
  **`ppdot`'s Crab validation never touches the real data path.** Both quoted numbers are hand-typed
  and in no committed artifact, and `fetch_atnf_ppdot` requests `PSRJ` and then discards it, so the
  run cannot identify the Crab at all. Its "Of the ~3500 entries" is also wrong: the table returns
  **2536** rows, the number `pulsarspec` records for the same call.

### Fixed
- **`triangulate` cited the wrong one of two same-year Krupar companions.** The entry was
  internally consistent --- DOI `10.1007/s11207-014-0522-x` matched its own volume 289, issue 8,
  pages 3121--3135 --- so a coordinate check passes it. But it is the *flux-density* paper, and the
  citation supports an apparent-source-size claim. The ~60 deg figure is in the goniopolarimetric
  companion (`10.1007/s11207-014-0601-z`), whose 125 kHz--2 MHz band is also exactly this paper's
  `\triFlo`--`\triFhi`. When two same-year companions share a title stem, the discriminator is the
  subtitle and the claim, not the identifier. The same sentence also stated ~60 deg as a constant
  where the source reports it at the *lowest* analyzed frequencies; now qualified.
- `triangulate`'s Method said channels are kept when the ray miss distance is "below a threshold"
  and never gave it. It is 60 R_sun, against inferred distances of 15.3--106.1 R_sun, so at the
  high-frequency end it admits rays missing by four times the distance being measured. Stated, with
  its permissiveness noted. Sweeping it needs a real-leg re-run: the committed results file keeps
  only summary scalars, not the per-channel arrays the headline is computed from.

### Fixed
- **`ecallisto_pipeline`'s only quantitative result was drawn from the wrong run.** Found by the
  referee round on the style conversion. `_write_macros` emitted seven un-namespaced names from both
  legs and `run()` wrote both legs to one results file, so the real 2011-09-14 archive day (which
  ran last) owned the macros while **all twelve macro uses in the paper describe the synthetic
  day**. The abstract typeset as *"a real burst at **0** stations ... confirms exactly **0** event
  and rejects 8"*, and the Method attributed the real day's `8 of 1512` and `-5.941 MHz/s` to a
  synthetic run over ten stations. `preserve_live_macros` was powerless: both legs wrote real values
  under one name. Macros are now `\ecSyn*`/`\ecReal*`, each leg writes its own results file
  (`ecallisto_synthetic_metrics.json`, allowlisted like `vgpra_synthetic_*`), and the paper cites
  the synthetic namespace. Correct values: 4 stations, 1 event, 3 rejected, 7 of 10 flagged.
  The committed real day is now cited once instead of sitting unread in `results/`.
- Two lessons worth keeping from that fix. **`preserve_live_macros` accumulates values, not names**
  --- it rewrites the lines the new run emits and drops any existing macro the new text never
  mentions, so a writer emitting only its own namespace *deletes* the other leg's numbers; each leg
  must emit both namespaces, filling its own and leaving the other as `--`. And **the un-namespaced
  test was the lock**: `test_run_offline_writes_artifacts` asserted `\ecNbursts` existed, so it
  passed because the defect was present and would have failed on the fix. Third instance in this
  repo of a test encoding "X must not exist" where X was the fix.

### Changed
- Papers restyled to traditional pre-LLM register (batch 7: `driftsearch`, `ecallisto_census`,
  `ecallisto_pipeline`, `frbperiod`, `offsets`) --- 35 of 45 papers converted. Gates clean on every
  file before review; **all five referee rounds returned revisions** (seven batches for seven).
  The MAJORs were all scope and claim-strength drift that no linter can see: `ecallisto_census`
  narrowed *"all event counts in **this paper** are synthetic"* to "in this validation" and turned
  *"This paper ... does not report a measured solar-cycle correlation"* into agentless passive,
  which reads as a claim about the literature and is false on the paper's own citation;
  `driftsearch` lost an inferential connective, so one detector on one file became *"A detector
  validated on injected tones ... fails on the real Voyager-1 data"* (any such detector), and
  softened "does not generalize **as-is**" to "as configured", offering an escape hatch the same
  paragraph rules out; `frbperiod` retitled *"A null, as expected"* to a declarative heading, in a
  paper that quotes no sensitivity whatsoever for that source; `ecallisto_pipeline` lost its
  limitations head sentence, so limitation (i) read as a description of the method, and promoted a
  LaTeX comment ("none typed by hand") into a body claim the paper falsifies one paragraph up.
- `driftsearch`'s title drops "Honest". The conversion removed all six body instances, leaving it an
  orphan; in journal register it is a claim about the authors, not the result, and `\shorttitle`
  already omitted it.
- Recurring mechanical signature this batch: **collapsing an em-dash appositive to a comma
  re-attaches it to the nearest noun.** In `offsets` it made the 24.0x excess a property of the
  Rayleigh expectation rather than the ratio of the two fractions; in `frbperiod` it let "on 19
  bursts" qualify the agreement claim rather than the periodogram. Grep for it in future batches.

### Fixed
- **Three wrong span claims in the solar-burst papers, two of them in abstracts.** Each was a
  round-number restatement of a macro pair that nobody re-derived after the runs settled.
  `triangulate` called its 0.125--1.975 MHz band "two decades of frequency" in four places; it is
  **1.20** decades (its distance range, 15.3--106.1 R_sun, is 0.84). The text now cites
  `\triFlo`/`\triFhi`, so the span cannot drift from the run again. `type3synthesis` claimed
  "more than three decades" for 0.125--78.94 MHz, which is **2.80**; the same paper counted "four
  single events" and "four regimes" where its own committed table has three, and its abstract
  claimed validation of a coronal density model the triangulated band never touches --- that leg
  is entirely interplanetary, so only Leblanc is tested, now stated in the limitations too.
- `solarbursts` printed `$0.137\,c$` for the harmonic/1x grid point three lines below
  `\sbSpeedC` = **0.1347** for the same point; that value now cites the macro. The two flanking
  grid values (0.086, 0.272) come from the same superseded run and are **not** silently patched:
  `exciter_speed` needs the raw ridge and only its summary is committed, so they cannot be
  recomputed from evidence. `survey/solarbursts-findings.md` tabulates the whole stale-vs-committed
  set (R^2 0.90 vs 0.811, drift -3.3 vs -2.55 MHz/s) and records the fix as emitting the grid on
  the next real run. Inventing two numbers to match would have been the worse repair.

### Changed
- Papers restyled to traditional pre-LLM register (batch 6: `solarbursts`, `swaves`, `windwaves`,
  `triangulate`, `type3synthesis`, `rmsky`, `sourcecounts`). Gates clean on every file before
  review; **all three referee rounds still returned revisions** (six batches for six), and all
  twenty-seven conversion defects are fixed here.
  The MAJORs: `sourcecounts`'s abstract stopped predicating "not a new measurement" on *this work*
  and attached it to the published counts it agrees with; `windwaves` turned "a recover-a-known"
  (the class of exercise) into "recovering a known result" two lines after saying its speed came
  out **below** published values, collapsed an appositive so the fitted speed attached to the
  *outer point* rather than the whole track, and downgraded "The honest limitations are
  substantial" to "Several limitations apply" in the paper whose dominant systematic exceeds a
  factor of two; `solarbursts`'s limitations list lost its head sentence, so a factor-of-3
  systematic now elaborates a positive clause.
  **One mechanical transformation had a signature**: `--- but they ---` promoted to a full stop in
  both `rmsky` and `sourcecounts`, leaving "They" pointing at the positive clause in each. A
  referee flagged the *transformation* rather than the instance, which made it greppable; it was
  confined to those two.


### Fixed
- **The two slices no guard could protect are now namespaced.** The audit's own output pointed
  at them: of the eight slices emitting no `*Source` macro, `vlass` and `type3synthesis` were the
  two with real numeric hazards, so `preserve_live_macros` could not even tell their run modes
  apart. `vlass`'s census counts are now `vlassReal*`/`vlassSyn*` (the live clobber was
  `\vlassNconfirmed` 2 -> 0 and `\vlassArea` 703 -> 0) and `type3synthesis`'s ladder is
  `synReal*`/`synSyn*` (corona speed 0.1347 real against 0.3002 offline). Both now emit their
  provenance, and `type3synthesis` records a `source` in its metrics --- it had none, which is why
  it was the single results file a forced offline rebuild could still change. Verified: a direct
  offline run leaves both papers' macros byte-identical, and the audit reports 0 for each.
- **The audit had two bugs of its own, both found by using it.** `re.IGNORECASE` made its
  `Syn|Real` marker match the `syn` *prefix* of every `type3synthesis` macro, silently exempting
  the whole slice --- a detector reporting a clean bill for a slice it never looked at. And it
  flagged the `*Source` macro itself, which is *supposed* to differ between modes since it is what
  the merge guard reads. Now case-sensitive, requiring an upper-case component boundary, with
  provenance macros exempt.


### Added
- **`scripts/audit_macro_namespaces.py`** --- the repo-wide grep for mode-dependent macros that
  are not namespaced. It runs every slice's offline leg into a throwaway directory and diffs the
  macros against the committed ones, so it measures the hazard rather than pattern-matching for
  it. **34 of 42 slices carry at least one**; about twenty change a *numeric* value under a
  shared name (`\hiVflat` 257 -> 231, `\rmRatio` 5.4 -> 8.4, `\vlassNconfirmed` 2 -> 0). Exit
  code 1 when any is found.

### Fixed
- **`preserve_live_macros` was defeated by `make figures`, silently, since it was written.**
  Snakemake **deletes a rule's declared output before running the job**, so the DAG removed each
  `generated/macros.tex` and the merge found no file to merge with. A forced run in the repo root
  **overwrote 511 real macro values across 39 papers** and blanked real values to `--`
  (`\ptRealNFit` 136 -> `--`), which the guard's own docstring says cannot happen. Direct
  invocation was protected throughout, which is why every previous test of the guard passed.
  The DAG now builds into `build/figures` and writes nothing into `papers/` or `results/` --- what
  a smoke build should always have done. Verified: a forced full `make figures` now leaves every
  committed macro file, results JSON, figure and CSV untouched, while still producing all 43
  artifacts under `build/`.
- **`preserve_live_macros` now also refuses a synthetic-over-real overwrite**, not just
  synthetic-over-real-with-a-placeholder. It reads the provenance every macro file already
  carries (the `*Source` macro) and, when a synthetic run meets a real file, keeps the existing
  values while still adding any macro the file lacks. This neutralises all 34 audited slices at
  once; per-slice namespacing remains the clean fix and is still worth doing.


### Fixed
- **`pte2`'s abstract asserted the opposite of its own body.** It listed "three facts undercut a
  giant-pulse interpretation" and gave as (i) that the excess does not correlate with spin-down
  luminosity, "contrary to the giant-pulse--Edot expectation" --- quoting the *raw* Spearman
  (-0.027), which the body establishes is count-confounded. Controlling for pulse count gives
  rho = **+0.098**, in the **predicted** direction, and the flagged sources' median log Edot sits
  **0.45 dex above** the rest, also predicted; neither is significant. The abstract now states
  the confounding, the partial correlation and the median offset, and concludes that the data
  neither support nor contradict the expectation --- which is what the body always said. No
  number changed.
- **`junodam` denied "an ~180x intrinsic rise" four lines after rendering `\jdNearFarRaw` =
  196.2.** The 180 was hand-typed and stale; now macro-backed.
- **`junodam`'s recover-a-known was not in committed evidence and its recovered value was
  wrong.** The paper hand-typed "injected 8.75, recovered 7.2" while `expected_contrast` was
  `null` and `\jdExpContrast` rendered `--`. Re-running the offline leg into a tmpdir gives
  8.75 -> **6.95**. Both are macros now, and the paper states that the injected contrast is
  several times the contrast measured on the real orbits, so the round trip cannot test recovery
  near unity.
- **`junodam` had a live macro clobber.** `io_contrast`/`expected_contrast` mean different things
  in the two run modes and shared one macro name, so an offline rebuild would have written the
  synthetic recovery into the macro the paper uses for the real measurement --- the `\tiiNEvents`
  shape, which `preserve_live_macros` cannot arbitrate because both runs write real values. Now
  `\jdRealContrast`/`\jdSynContrast`/`\jdSynExpContrast`, with a test that reproduces the clobber.
  Two existing tests asserted the un-namespaced names, i.e. they were pinning the defect.
- `junodam` now cites `louis2021jgr`, which sat in `refs.bib` uncited while measuring the
  latitudinal beaming this paper invokes as its preferred interpretation.


### Changed
- Papers restyled to traditional pre-LLM register (batch 5: `pte2`, `junodam`, `rmdipole`, `hi`,
  `glitchpop`, `vlass`). Gates clean on every file before review; **all five referee rounds still
  returned revisions** (five batches for five), and all twenty-two conversion defects are fixed
  here.
  **The first defect in this campaign where a style edit changed a FACT rather than a claim's
  strength**: folding a run-in label into a topic sentence, the `glitchpop` conversion wrote
  "the full **post-2018** catalogue" where the code (`glitchpop.py:283`) reclassifies on "the
  full catalogue". That describes a disjoint-epoch comparison never run, and it makes both
  derived quantities incoherent -- `newly classifiable` is a union operation and the single flip
  is recorded as `n_pre: 8 -> n_now: 9`, i.e. cumulative. No number, macro or citation moved, so
  no mechanical gate could see it.
  Also: `pte2`'s "Honest bottom line." became `\subsection{Summary}`, promoting to the paper's
  designated summary a paragraph quoting only the naive 19% incidence -- the endpoint its own
  referee round established is not the answer; `hi`'s closing paragraph lost the word
  *validation*, so its contribution class appeared nowhere after Results; `junodam` lost the only
  heading naming a limitation, in a paper whose one positive claim is explicitly an upper bound;
  and `vlass` swapped "pipeline" for "selection" in a sentence whose predicate is false of the
  selection alone.

### Fixed
- `glitchpop`'s abstract had lost "among these" (the denominator: the 1.55 expected false
  positives are 31 x 0.05) together with the two words distinguishing the *quasi-periodic
  fraction* from *individual members*, in a single edit.
- `rmdipole`'s abstract had split its significance claim from the caveat that kills it, so the
  formally significant p = 0.001 dipole could be excerpted without the clause establishing that
  the rejection is carried entirely by the top 1% of |residual|.


### Added
- **`report.preserve_live_results` / `report.write_results`** --- the results-JSON counterpart of
  `preserve_live_macros`, wired into **all 49** `results/*.json` write sites across 46 modules.
  A results file may now only gain information: a real result is never replaced by a synthetic
  one, a partial real re-run retains the other leg's keys rather than dropping them (the
  `torchfdmt` GPU-column case), and any key carried across runs is listed under `_merge` so a
  spliced row cannot look like one invocation. Sixteen tests.
  **Verified against the documented incident**: `typeii.run(".", offline=True)` --- the command
  CLAUDE.md records as deleting 3429 lines and flipping `is_real` True->False --- now leaves
  `results/typeii_metrics.json` byte-identical, while the offline run still writes its `tiiSyn*`
  macros. A *forced* full `make figures` in the repo root changed **25** committed results files
  before and **1** after; that one carries no provenance marker at all.
- `make guard-real` now reports how many results files carry no `source`/`is_real` marker (76
  today), so the gap neither it nor `preserve_live_results` can cover is visible rather than
  silent.

### Fixed
- The first cut of the results guard protected **3 files out of 25** because it keyed on
  `is_real`, which most slices never set. It now uses the `source` field that
  `guard_real_results.py` already treats as authoritative, and counts a *mixed* source
  ("synthetic recover-a-known + real RACS-mid epoch pair") as real --- otherwise
  `stokesv_discovery` was left unprotected by the very string that documents its real leg.
- Seven `results/*.json` writers outside the `*_metrics.json` naming convention were missed by
  the first pass, including `hi.py`'s `rotation_curve.json`. A forced offline rebuild overwrote
  it with synthetic output and **`make guard-real` failed** --- caught only because the forced
  rebuild was run as a test rather than assumed to be unnecessary.


### Fixed
- **`torchfdmt` quoted "~24x" for a ratio its own adjacent macros make 29x.** The abstract and
  benchmark section rendered "~24x (44.12 -> 1.5 s)"; 44.12/1.5 = 29.4. The 24 traces to a
  superseded 36.1 s CPU timing that appears nowhere in `results/`, in a file whose header claims
  every number is `\input` from the pipeline. Fixed at the source: `_write_macros` now
  **derives** `\spBruteSpeedup` and `\spFdmtSpeedup` from the committed timings, so no ratio in
  this paper can be hand-typed or drift again; a CPU-only run emits `--`, which
  `preserve_live_macros` will not write over a real value. Three regression tests, one checking
  the macro on disk against the JSON on disk.
- **`torchfdmt`'s real-data DM recovery was quoted as a bare "0.3%".** The FDMT butterfly indexes
  rows by whole samples of dispersive delay, so the DM is quantised rather than fitted: at this
  file's band and sampling one row is 0.0627 pc/cm^3, making the 0.18 offset **2.9 trials**. The
  abstract's "within delay quantisation" frame belonged to the synthetic leg and was inheriting
  to the real one. The real leg now records `real_dm_step_pc`, the paper says the DM is
  grid-quantised, and it quotes `\spRealSnr` (6.0) -- generated and used nowhere, despite being
  the number that says how well the peak can be localised.
- **`torchfdmt`'s benchmark row is a splice of two invocations.** `git log -p` shows a CPU-only
  run with GPU values patched in later while `device` stayed `"cpu"`; the code's own device
  logic means a `--device cuda` run writes both columns *and* sets `device: "cuda"`, so no single
  run can produce the committed row. The conclusion survives on same-run numbers (3.3x), so this
  is provenance: the paper now says the columns come from separate invocations, and
  `bench_devices` is recorded going forward.
- **`results/singlepulse_metrics.json` labelled its synthetic block with the real file's name.**
  `source` sat directly above `recovered_dm: 56.63`, a synthetic value (the real one is
  `real_recovered_dm: 56.59`). `source` now names both legs and states the key convention.


### Changed
- Papers restyled to traditional pre-LLM register (batch 4: `typeii`, `vgpra` main+note,
  `spectra` main+note, `frbstats`, `peaked`, `torchfdmt` --- nine files). Gates clean on every
  file before review; **all five referee rounds still returned revisions**, and all twenty-eight
  conversion defects are fixed here.
  **The sharpest case in the campaign so far: `typeii`'s "100%" changed meaning without changing
  value.** At HEAD the literal sat inside a disavowal ("This curve --- not a headline `100%' ---
  is the honest characterization"); after conversion it was an assertion ("recovered at 100%
  completeness"). The numeric multiset is identical, so the guard passed --- and the disavowal
  was not stylistic, it was forced by GATE-2 finding R2 because the 1.0/1.0 is an easy-synthetic
  ceiling. Also: `vgpra` lost "The result is a controlled null" and the word *controlled* fell to
  zero occurrences while its own note still used it; `peaked` lost the only *comparative* form of
  a robustness claim that survives twice in unscoped form; `spectra` and `frbstats` each lost the
  single skim-level carrier of a negative result.

### Fixed
- **A style rule that was wrong.** The `peaked` conversion stripped `\emph{}` from the mission
  name *Fermi* to stay under an `\emph`-per-kw threshold, but the style guide's own reference
  says italics are *for* mission names. Fixed with `\textit{Fermi}`. When the lint and the guide
  disagree the guide wins; when the guide and a science gate disagree the science gate wins.
- `papers/peaked/arxiv.yaml` re-derived after the conversion. Its overridden arXiv metadata
  abstract still asserted "far more robust ... and alpha high is TGSS-independent" where the
  paper now says "more robust ... because alpha high is TGSS-independent"; a metadata override
  that no longer matches its paper is worse than none.


### Fixed
- **`torchdsp` claimed its science ran on the GPU; the committed evidence says CPU.** Settled
  from git history: the commit that created `benchmark_device` (602e0ca, 2026-08-10) states
  that one field labelling both legs meant "a real GPU benchmark could only be stored by
  mislabelling the CPU-run CHIME science as GPU". `results/torchdsp_metrics.json` agrees
  (`device: cpu`, `benchmark_device: cuda`). Exactly four numbers are GPU-backed, all wall-clock
  timings. The paper's "run entirely on the ROCm GPU", "validated end-to-end" and "runs
  identically everywhere" are corrected to what the evidence covers, and the same wrong
  attribution is fixed in `README.md` (2 rows) and `survey/torchdsp-findings.md` (3 places).
  Cross-device numerical agreement is not measured and the paper now says so.
- **`stokesv` printed a reproduction command that destroyed committed evidence.** `--offline`
  defaulted to False and `--out` to `"."`, so the command in the paper ran the *real* path from
  the repo root. `--out` is now required, and the synthetic figure is no longer drawn on the
  real path at all -- previously a real leg that raised partway (a missing CASDA credential is
  enough) left synthetic output standing where the committed two-panel real figure had been.
  Two regression tests.
- **`southern`'s macros were unmerged AND un-namespaced**, the combination that ships a wrong
  number rather than a hole. One `make figures` would have replaced 1545 real matches with the
  synthetic field's count and written `\soCallTried{0}` over a real 50. Counts are now
  `soSyn*`/`soReal*`, the real-only Callingham validation defaults to a `--` placeholder rather
  than 0, and the writer merges via `preserve_live_macros`. Verified empirically: an offline
  rebuild now fills `\soSynNmatched{40}` and leaves `\soRealNmatched{1545}` untouched.
  The existing `test_run_offline` had asserted the *un-namespaced* names, i.e. it was pinning
  the defect -- another instance of a test locking a bug in.
- **Audited the `preserve_live_macros` invariant across every slice, not just the ones in
  hand.** Four more were missing it: `atlas3i` (three macro writers, and it heads the
  submission queue), `driftsearch`, `peaked`, `stacking`. All fixed. **42 of 42 slices now call
  it**, against 5 of 42 at the 2026-08-12 audit. The CLAUDE.md lesson stands: a rule stated in
  that file is not a rule the repo follows until something greps for it.


### Changed
- Papers restyled to traditional pre-LLM journal register (batch 3: `rmstructure`, `torchdsp`,
  `stokesv`, `southern`, `lpt`, `skr`). Lint 3/2/3/3/3/2 HIGH -> 0 on all six; diff-guard clean;
  triage unchanged. **All six referee rounds still returned revisions**, and all twenty-nine
  conversion defects are fixed here. The severe ones were single deleted qualifiers that invert
  a headline: `skr`'s abstract lost **"raw"** from "a raw near/far ratio of 3.33" (the paper
  exists to show that ratio is a 1/r^2 artefact collapsing to 1.39); `lpt`'s abstract lost
  "own", "already" and **"per-value provenance"** -- the contribution named in its own title --
  from its precedence concession to Rea et al.; `rmstructure` had its quality-flag result
  demoted from "a lesson promoted to the main result" to a methods aside across four edits;
  `southern` *gained* an unsupported "directly" on a fitted extremum and had "The recovery"
  become "The recovery **fraction**", a quantity the pipeline cannot produce; `torchdsp` lost
  "as of this writing" from a novelty claim its own Introduction contradicts, and fused a
  provenance sentence so synthetic numbers inherited real-archive attribution.

### Added
- `traditional-style` skill rule: **run-in labels often carry ranking.** Dissolving a
  `\textbf{Label.}` opener into a topic sentence, or promoting it to a `\subsection{}`, changes
  how the reader weights what follows -- batch 3 turned "An honest null:" into a heading naming
  an attempted measurement as though it were a claimed one, and lost a "(a near-null)" marker so
  a paragraph read as a detection for three sentences.


### Fixed
- **`frblens` quoted its own retracted limit as its first conclusion.** Discussion item (i)
  read "not lensing-contaminated at the $\gtrsim$9\% level" as a hard-coded literal --- that
  is `\flRealLimitNaive` = 0.091, the count-based value this slice's referee round retracted
  in favour of 0.368, and which the Results paragraph calls "wrong by a factor of four, in the
  optimistic direction". Now macro-backed so it cannot drift again.
- **`frblens`'s printed equation was not reproducible.** `\flRealEpsSum` (8.3) and
  `\flRealEpsMean` (0.24) were computed over 34 sources while the search ran on 33 and the
  limit divides by the restricted sum, so the page showed `2.996/8.3 = 0.3683` where that
  quotient is 0.3624. The extra source is FRB20210601A: `injection_efficiency` accepts any
  train with >= 5 bursts, the search also requires a span > 2 d. Only the macros were wrong ---
  the limit was always computed on the searched subset --- so `_write_macros` now restricts to
  `m["rows"]` and the macros were regenerated from the committed metrics with no re-run.
  8.3 -> 8.1, 0.24 -> 0.25; 2.996/8.1333 = 0.3684. Two regression tests, one of which checks
  the committed evidence so page and JSON cannot drift apart again.
- **`stokesv_discovery`'s `\svdRealDet` counted neither systems nor rows.** The paper said it
  "counts systems rather than rows" while the macro is 2: the four >= 5-sigma rows are the two
  CNS5 entries of one binary in a single beam, so the chain is 4 rows -> 2 entries -> 1 system.
  The sentence now states that chain and the estimator records what it counts.
- **`stokesv_discovery`'s census is a 12-arcsec peak search and the Method claimed otherwise.**
  `search_arcsec` defaults to 12 ("brightest Stokes-I pixel within that radius"); only
  `<= 0` is forced. Re-verified in the committed evidence: **all 54 quiescent rows have I > 0**
  at median I/sigma_I 2.19, where a fixed-pixel measurement would go negative about half the
  time (p = 2^-54). The Method now describes the search, the two remaining "forced" method
  claims are gone, and Honest limits states that the quiescent limits bound a beam-scale
  maximum rather than the star. The forced re-measurement remains outstanding.
- **`wdpulsar` quoted a Stokes-V limit as the depth of the Stokes-I null** in the abstract and
  conclusion (\wdRealMedianVLimit is `median_v_limit_mjy`); section 3 always labelled it
  correctly. Both now say Stokes-V.
- **`frbwait`'s table caption said "top 20 of 15"**; the table has all 15 rows above the
  completeness cut and is not truncated.

### Changed
- Papers restyled to traditional pre-LLM journal register via the `traditional-style` skill
  (batch 2: `svsbi`, `stokesv_discovery`, `frblens`, `wdpulsar` main.tex, `frbwait` --- the
  first three previously refereed, so their science was frozen). Lint 3/3/2/2/2 HIGH -> 0 on
  all five; diff-guard clean; triage unchanged. **All five referee rounds still returned
  revisions**, and all sixteen conversion defects are fixed here: `frbwait`'s Conclusions had
  lost the clause scoping "clustering is ubiquitous" to its own census (3 of 15 sources),
  `wdpulsar`'s abstract had collapsed a two-of-five control record into an unqualified
  re-detection and had recast a VLASS catalogue cross-match as a photometry validation,
  `svsbi` had promoted a `though` to a full stop so a power caveat stopped bounding a
  calibration pass with a 0.004 margin, and `frblens` had widened a novelty claim by dropping
  "observational".

### Fixed
- **The diff-guard can be satisfied by relocation.** Shortening `frblens`'s abstract, the
  style pass moved `$\pm$230\,s` into the Discussion, preserving the numeric multiset by
  construction so the guard passed --- and the number landed in a clause stating the Roemer
  term is "~30x our tolerance ... reaching +/-230 s there", where the ratio is 46, and which
  the Method attributes to a different delay. The abstract meanwhile lost the reason
  barycentring is mandatory. Both restored; the skill now forbids relocating a numeric literal
  across sections. Multiset equality is not content equality.

## [1.8.0] - 2026-08-22

### Added
- `lptxray` (plan 93): uniform X-ray cross-match of the 16 LPTs and the 56 white-dwarf-pulsar
  candidates against 5XMM-DR15, 2RXS, eROSITA-DE eRASS1 and CSC 2.1.1, with the
  chance-coincidence rate **measured** two independent ways (field density, and rigid
  position-shift trials computed offline from one cached cone per position; both ~1%).
  Among the candidates, **accretion strongly predicts X-ray brightness -- 20/21 accreting
  (95.2%) against 3/35 others (8.6%), ratio 11.1** -- and survives a common-footprint cut
  (14/14 vs 3/24) and an all-sky 2RXS-only cut (18/21 vs 0/35). Set against `wdpulsar`'s
  radio photometry at matched depth (0.818 vs 0.825 mJy), **accretion does not predict radio
  loudness: 0/19 accreting and 0/30 others**, though this constrains persistent emission only
  and cannot exclude LPT-like pulses. Three X-ray associations absent from Pelisoli's
  compilation, including the confirmed WD pulsar J1912-4410. Evidence:
  `results/lptxray_metrics.json`, `survey/lptxray-findings.md`.
- `lptxray.classify_coverage`: splits archival pointings into targeted and serendipitous by
  aimpoint offset rather than target name -- all 24 pointings near ASKAP J1935+2148 are
  observations of SGR 1935+2154, whose designation shares the same RA digits. Census: 12 of
  16 LPTs have archival X-ray coverage, 8 have a dedicated pointing.

### Fixed
- Plan 93's LPT leg, **cut at GATE 0 on measured evidence**: all three LPTs with a published
  X-ray detection are absent from every serendipitous source catalogue while each has pointed
  coverage taken *after* those catalogues were built (XMM 0953011101, XMM 0973390301, Chandra
  26681/26682/29265/29266). Catalogue recall on the LPT sample is **1/3**, against 16/16 on
  the candidates, so a catalogue-derived LPT detection fraction would measure catalogue
  latency rather than the sources. `lptxray.catalogue_recall` enforces this as a guard, and
  the two legs are never differenced. The plan's headline was also overstated: only one of
  the 16 LPTs is a confirmed accretor.

### Added
- `lptlin` (plan 92): linear polarization at LPT pulse positions from VAST Stokes Q/U.
  One detection -- ASKAP J1832-0911 at L/I = 10.9% (33 sigma, EVPA 98 deg) -- plus three
  non-detections. **The 10.9% is depolarization-dominated, not the source's linear fraction**:
  the published RM of +89 rad/m^2 leaves only ~9% of the intrinsic linear polarization in a
  band-averaged image, so this is a strong lower limit consistent with an intrinsically ~90%
  polarized source (matching the published 92+/-3% total). Ricean debiasing, a leakage veto,
  a `total_polarization` bound check and a `faraday_depolarization` calculator. Evidence:
  `results/lptlin_metrics.json`, `survey/lptlin-findings.md`.

### Fixed
- Plan 92's premise: an earlier feasibility check generalised from a query at **one** position
  to conclude that only EMU carries Q/U at LPT positions. Querying all sixteen shows 14 have
  Q/U imaging and VAST produces it too. The plan's scope reduction was wrong and is corrected
  in place.

### Added
- `lptspec` real run and **negative result** (plan 91, now closed): Taylor-term cutouts were
  staged for the three pulses GATE 0 cleared, and the method fails. alpha comes out -1.24,
  -2.99 and -17.30; the last is the taylor.1 image's *global minimum* at the source pixel
  (-49 sigma), and a recover-a-known over the same image puts every steady source at
  |T1/T0| ~ 1-2 against our 17. MFS fits one constant-flux power-law source across the whole
  synthesis, so an LPT pulse present for only part of it cannot be represented and the
  deconvolution absorbs the mismatch into taylor.1 -- invalid rather than merely biased, and
  not fixable by choosing brighter pulses. `results/lptspec_metrics.json` is retained as
  evidence of the failure and stamped "do not quote the alpha values".

### Added
- `lptspec` GATE 0 (plan 91): `src/jansky_research/lptspec.py` + `scripts/lptspec_gate0.py`
  decide, without fetching an image, whether an ASKAP Taylor-term in-band spectral index is
  recoverable for any LPT pulse. The short in-band lever arm makes the Taylor-1 image ~10.7x
  noisier than Taylor-0, so `sigma_alpha ~ 10.7/(S/N)`; folding in the published MT-MFS
  penalty (Rashid et al. 2024, arXiv:2405.18978) gives **3 of 7 pulses usable** at
  `sigma_alpha <= 0.3` — including ASKAP J174508.9-505149, the source whose reported
  frequency drift motivated the plan. Novelty pass clear: no LPT paper reports in-band alpha
  from Taylor terms. Evidence: `results/lptspec_gate0.json`, `survey/lptspec-findings.md`.

### Added
- `papers/lptduty/rnaas.tex`: the plan-90 note (RNAAS register, macros generated from the
  committed JSONs by `lptduty.write_paper_assets`, prose-lint and triage clean). Leads with
  the ephemeris finding: the measurement is limited by what discovery papers publish, not by
  sensitivity.

### Fixed
- `scripts/triage_papers.py` checked only `paper/main.tex`, so **no RNAAS note had ever been
  triaged** -- including wdpulsar's, at the head of the submission queue. It now checks every
  file carrying a `\documentclass`. All five previously-unchecked notes come back clean.

### Fixed
- `lptduty` GATE-2: the Poisson upper limit used `2.996 + k`, correct only at k=0. The exact
  one-sided 95% limit is `0.5*chi2.ppf(0.95, 2k+2)` (2.996, 4.744, 6.296 for k=0,1,2), so the
  three sources *with* a detection had limits 19-26% too tight -- and the write-up quoted only
  the unaffected zero-detection limits, so it was invisible. Corrected to 0.0521, 0.0668,
  0.1431 (and 0.2965 in the phase leg); `poisson_upper_95` now pins the exact values.
  Also: the detection criterion now applies `lptv`'s leakage veto (|V| > 0.006|I|) as its
  docstring already claimed (no committed number changes), and GATE 0 records the family-wise
  alpha it actually tests against.

### Fixed
- `data/lpt_sample.csv`: five errors found by the plan-90 ephemeris audit and verified
  against the arXiv record before changing anything -- three wrong discovery arXiv ids
  (GCRT J1745-3009, GLEAM-X J162759.5, GPM J1839-10 all cited unrelated papers), ASKAP
  J1832-0911's period derivative (9.0e-12 -> 9.8e-10, two orders of magnitude) and its
  period (2656.2554 -> 2656.247, the value its own cited paper gives). The three arXiv ids
  had propagated into `papers/lptv/refs.bib` notes; triage's Crossref DOI check could not
  see them because the titles and DOIs were correct. No committed `lptv` number changes.

### Added
- `lptduty` increment 3 (`scripts/lptduty_phase.py`, `results/lptduty_phase.json`):
  phase-resolved separation of f_active from the in-period duty cycle. Validated by
  reproducing exactly the pulse phases (0.489, 0.951) that `lptv` published independently.
  Seven of ten sources are excluded for having no published reference epoch at all; of the
  three attempted, one supports the single-window model (f_active ~ 0.06).

### Added
- `lptduty` GATE 0: novelty pass (clear -- no class-wide LPT activity-fraction constraint
  exists; the field review says the selection function is unmodelled) and an aliasing test
  (`scripts/lptduty_gate0.py`, `results/lptduty_gate0.json`). The aliasing test invalidates
  two of ten constraints: GPM J1839-10's snapshots are not uniform in pulse phase
  (Kuiper V 0.275, Bonferroni p 0.012), and ASKAP J142431.2-612611's period is quoted too
  imprecisely for the test to run. Both verdicts are stamped per source into
  `results/lptduty_metrics.json` (`constraint_valid`, `phase_sampling.verdict`).

### Changed
- **Airflow 3.** The container base moves 2.9.3 -> 3.3.1, the `airflow` extra to
  `>=3.3,<4.0`, DAG imports to the Task SDK (`airflow.sdk`), and compose from `webserver`
  to `api-server` (removed in 3.0) with an explicit `EXECUTION_API_SERVER_URL` and
  SimpleAuthManager (`airflow users create` belonged to the removed FAB auth manager).
  Verified by running the DAG: `airflow dags test ecallisto_ingest 2011-09-14` completes
  with all 7 tasks successful under 3.3.1.

### Fixed
- DAG `ecallisto_ingest`: `scan_station` now returns JSON-safe values. Airflow 3 serialises
  XCom as **strict** JSON, so a NaN in the returned row raises "Out of range float values
  are not JSON compliant" and fails the task; Airflow 2's encoder allowed it. The DAG parses
  cleanly either way, so this only appears at runtime on stations whose scan yields a NaN --
  found by running the upgrade rather than by CI, which never exercises this stack.

### Added
- `lptduty` (plan 90, increments 1--2): constrains how often each LPT is on, from the 647
  measured VAST snapshots `lptv` already committed -- no new data. Reports
  `p = f_active x (w + T) / P` (whose factors are not separately identifiable without an
  ephemeris) with an efficiency-weighted denominator rather than an epoch count, and an
  efficiency floor so shallow epochs cannot sum into sensitivity that does not exist.
  Three sources give p = 0.011--0.046; seven give 95% limits of 0.029--0.107 at 5 mJy.
  Evidence: `results/lptduty_metrics.json`, `survey/lptduty-findings.md`. GATE 0 and the
  aliasing check are open -- not paper-ready.

### Added
- Plans 91--93, the follow-ups Rose et al. 2026 suggests, each gated on a feasibility check
  run before the plan was written rather than after: `lptspec` (in-band spectra of LPT
  pulses from the ASKAP Taylor-term images -- 98 `taylor.1` products confirmed present at an
  LPT position, so no sub-band re-imaging is needed); `lptlin` (linear polarization --
  **scope reduced**: a CASDA query showed RACS carries no q/u at these positions, only EMU
  does, so the intended all-sky census does not exist to be done); `lptxray` (does accretion
  predict radio loudness -- cross-match the LPT catalogue and the Pelisoli candidates against
  public X-ray catalogues, with a measured chance rate).

### Fixed
- `lptv` refs: `rose2026` upgraded from an arXiv e-print to the published record
  (Nature Astronomy 10, 1166--1178, `10.1038/s41550-026-02882-x`), verified against
  Crossref. It is the discovery/classification paper for ASKAP J174508.9-505149, one of
  the census's three RACS detections, and `lptv` is in the submission queue.

### Added
- `plans/90-lpt-duty-cycle.md`: quantify how often LPTs are on, from the 966 committed VAST
  snapshots — Rose et al. state the switch-off behaviour qualitatively for one source and
  nobody has measured it across the class. Plan records up front that only the product
  (active fraction x in-period duty cycle) is identifiable without phase information, and
  that the denominator is efficiency-weighted exposure rather than epoch count.

### Fixed
- `frbstats`: annotate the repeater/one-off property dicts as `dict[str, np.ndarray]`.
  A comprehension over literal keys infers `dict[Literal[...], ...]`, which dict
  invariance rejects against `compare_populations`' parameter — an error mypy >=2.3
  reports and 1.11 did not, and the only thing blocking the dev-tooling bump (#139).

### Fixed
- CI: the tectonic package bundle is now cached (`actions/cache` on `~/.cache/Tectonic`,
  rolling key) in both `paper.yml` and `release.yml`. Every run previously re-downloaded
  the whole bundle from `relay.fullyjustified.net`, which 429s under load — it failed the
  `paper` and `release` workflows on 2026-08-19 (the v1.7.0 tag push, so the release had
  to be created by hand) and flaked twice the day before. After one successful run the
  compile step needs no network. Restore and save are split so the save runs on
  failure too — a 429-interrupted run still downloads part of the bundle, and keeping
  that partial cache lets successive runs converge instead of each starting cold.
  The cached path is `~/.cache/tectonic` (lowercase; `~/.cache/Tectonic` is also listed
  for older builds) — the first attempt named only the capitalised form, which cached
  nothing while the step still reported success, because `actions/cache` logs a missing
  path as a warning rather than an error.

## [1.7.0] - 2026-08-19

### Changed
- Papers restyled to traditional pre-LLM journal register via the `traditional-style` skill
  (batch 1: atlas3i main + RNAAS note, dr20radio, lptv, innerrc — the submission queue).
  Prose only: every macro, citation, and numeric literal is multiset-identical to its
  pre-conversion state (`prose_lint.py --diff-guard`). Each conversion then went through a
  `paper-referee` round, all four returning *minor revision* for scope drift the guard
  cannot see (deleted qualifiers, an appositive that re-pointed a claim at the wrong
  subject, negotiated verb strengths, ranking markers) — all restored; see
  `survey/stylecorpus-findings.md`.
- `wdpulsar` figure: the left panel is now real data (forced I/sigma over the 442 usable
  candidate epochs against the unit normal), replacing the synthetic injection fixture —
  the referee's deferred suggestion. The synthetic panel remains as the offline fallback.

### Fixed
- `CHANGELOG.md` structure: entries from seven PRs had accumulated *inside* the intro
  sentence rather than under `## [Unreleased]`, so `scripts/next_version.py` reported
  "nothing to release". Entries relocated and regrouped; the release recipe works again.

### Added
- `traditional-style` RNAAS extension (plan 89 follow-on): genre baseline from 391
  pre-LLM arXiv-deposited notes (`results/stylecorpus_rnaas.json`; the main corpus's
  refereed filter had excluded all 1,035), style-guide §8, `prose_lint.py --file/--genre`,
  and the wdpulsar RNAAS note converted through the full gate sequence (referee round
  restored two qualifiers the multiset guard cannot see; blind A/B passed).
- `traditional-style` skill + `style-editor` agent (plan 89, Stages 2–4): 2,000-paper
  pre-LLM corpus acquired (manifest committed), fingerprint engine + corpus-vs-us delta
  (`results/stylecorpus_{fingerprints,selfscan,manifest}.json`), corpus-derived style
  guide with before/after conversions, and `prose_lint.py` (corpus-percentile lint +
  prose-only `--diff-guard`).
- `stylecorpus` (plan 89, Stage 1): scoping of the pre-LLM radio-astronomy literature
  (1933–2021) for the `traditional-style` skill — `src/jansky_research/stylecorpus.py`
  (era strata, ADS/arXiv query builders, stratified size sampler, bootstrap total-PDF-size
  estimator; network runners no-cover), `scripts/style_corpus_scoping.py` (resumable CLI),
  committed evidence in `results/stylecorpus_scoping.json` + `survey/stylecorpus-findings.md`.
- `papers/wdpulsar/rnaas.tex`: standalone RNAAS note (refereed; ~630 words) — the null
  census condensed with the control's true 2/5 efficiency attached to the f<6.1% bound.
- `lptv.summarize_vast` + `fold_phase` (tested): VAST-sweep reduction with per-detection
  pulse phases on the published J1839-0756 ephemeris; paper gains a VAST-extension section
  reporting two pulses absent from the discovery paper's observation table (one at
  interpulse phase 0.489±0.004), with unreleased-SBID bookkeeping (107 epochs never
  publicly released, verified via obscore) and the honest 2026 phase-scrambling caveat.
- `scripts/lptv_vast_real.py`: VAST-collection extension of the lptv forced-V sweep
  (12-min epochs, ~966 across 10 covered sources; full-precision epoch MJDs and per-row
  durations, fixing the round-2 timestamp-quantization finding).

### Fixed
- `wdpulsar` main.tex: two peak-search-era numbers corrected against the forced CSV
  (max candidate significance is 2.9σ not 3.3–3.9σ; brightest measurement 3.7 mJy not
  4.3 mJy) — found by the RNAAS note's referee pass.
- Triage: a deliberately-committed `*synthetic*` metrics file with `is_real: false` no
  longer flags no-evidence when the real sibling exists (the vgpra two-file design);
  16.11 h allowlisted as a Voyager literature constant. vgpra "confirms"/"proves" and
  wdpulsar hard-typed 7.4/"confirms" fixed; zero HIGH/MED across all 46 papers again.
- `lptv` round-3 referee (VAST section): 2026-epoch chain recomputed self-consistently in
  code (20 epochs, p=0.53); phase errors include exposure smearing + anchor (±0.016);
  interpulse claim downgraded to candidate on the published V/I tension; novelty pass run
  and recorded with all section evidence in `results/lptv_vast_adjudication.json`; sweep
  end date corrected (2026 July, not August); leakage veto |I| guard; boundary tests;
  three bib entries completed/corrected (lee2025, mcsweeney2025, hurleywalker2023).
- `lptv` detection provenance: all three census detections adjudicated against their
  discovery papers' e-prints — each is an epoch the discovery paper already reported
  (J1745's SBID 20398 discovery-search epoch; J1651's single archival detection of
  2024-11-21; J1839's discovery observation). Paper reframed: zero new detections, three
  validated re-measurements (`results/lptv_detection_provenance.json`).
- `lptv` round-2 referee pass (major): rescoped negative recomputed from the forced CSV
  (median I/sigma = -0.14, no |V|/I constraint on non-detections — the 1.9 / 49–114% figures
  traced to the superseded peak-search data); J165130.3 promotion re-based on locked-pixel
  flux recovery; J183950.5 bright epoch identified as the published discovery observation
  (SBID 57929, lee2025) and reframed as a guaranteed (non-blind) known-pulse recovery with
  committed adjudication evidence (`results/lptv_j183950_adjudication.json`); ASKAP
  J1832-0911 catalogue row's erroneous name corrected (coordinate always matched the VLBA
  position); stale two-detection text, figure caption, epoch-composition, limit-count, and
  phase-coverage errors fixed.

## [1.6.0] - 2026-08-13

### Added
- `scripts/triage_papers.py` — a mechanical first-pass review over all 46 papers, checking
  only defect classes found in real papers this cycle: cited `--` placeholders, un-namespaced
  mode-dependent macros, hard-typed numbers duplicating macros, DOIs whose Crossref metadata
  contradicts the entry, `author = {others}`, missing evidence, and retracted verbs. Ships
  with an evidence-alias map for pre-convention filenames and a reasoned allowlist for
  adjudicated coincidences. Final state: zero HIGH, zero MED across all 46 papers.
- `svsbi`: `prior_sensitivity` — the same census re-inferred under a widened prior box over
  multiple torch seeds (torch was previously unseeded, so NPE training was irreproducible).
  Verdict per parameter is committed: `log_Lbreak` is prior-driven (+0.96 dex when the box
  widens, against 0.02–0.14 seed scatter), `f_beam` is prior-located, only the slope median
  is data-driven. The prior box itself is now in the evidence file.
- `stokesv.measure_circular_pol`: a genuinely forced mode (`search_arcsec <= 0`) that reads
  the pixel at the propagated position. The default 12″ mode is a peak search — on blank sky
  a noise-maximum statistic (54/54 quiescent census targets had I > 0, p = 2⁻⁵⁴ for a true
  fixed-pixel measurement) — and the docstring now says which is which. A test pins the
  distinction: forced photometry on noise goes negative about half the time.
- `frblens`: `efficiency_per_source` — injection efficiency for every searched train at the
  census's own detection threshold, and `lensed_fraction_limit(..., eps_sum=)` divides by
  Σε rather than the source count.
- `fashienv`: `void_jackknife_offset` — delete-one-void jackknife on the knee offset, plus
  `n_wall`/`n_field` committed (the comparison bin held 58% of DR1 and its size was nowhere
  in the evidence).

### Changed
- **`svsbi`'s headline is retracted to a lower limit**: the one detection does not "pin"
  log L\* (the posterior piled against the prior wall and its median tracked the box); the
  census supports log L\* ≳ 13.9. `f_beam` is stated as the product of emitter fraction,
  beaming and duty cycle — the model cannot separate them. `\svbNTargets`/`\svbSource` are
  namespaced (an offline rebuild previously wrote the 400-star synthetic parent size under
  the macro the abstract used for the 38-target real census).
- **`stokesv_discovery`'s 10σ inter-epoch decline is demoted to marginal**: I and V fell
  together (V/I constant to 3.4%), the flux-scale signature; 5.2σ/3.5σ at a 3%/5% per-epoch
  scale term. The title's "Variability Recovery" is now "Recovery". One system, not two rows:
  BL and UV Cet share a beam and the photometry was bit-identical.
- **`frblens`'s limit is 4× weaker and now honest**: f < 0.37, dividing by measured Σε = 8.3
  (mean efficiency 0.24; four searched sources have ε = 0 and constrain nothing). The
  abstract's lensed-one-off framing is corrected to lensed repeaters — the ≥5-burst and
  M_max ≥ 2 cuts exclude the one-off channel entirely (30/33 sources have M_max = 1).
- **`fashienv`'s offset is reframed as an upper bound**: the void jackknife *strengthens* it
  (0.039 dex vs the 0.087 fit error — sample variance is not the limiter), so the paper now
  names the two biases that are, both acting in the signal's direction: the 1/Vmax
  uniform-density assumption, and a bounding-box comparison bin 0.010 dex from the all-sky fit.
- 33 modules now route `_write_macros` through `report.preserve_live_macros` (5 of 42 did,
  against the stated repo rule); a single `make figures` would previously have blanked every
  real macro in 33 papers. All 46 papers now carry the `\software{}` block citing
  `jansky-research` (4 did).

### Fixed
- 19 citation defects across 12 papers, each verified against Crossref/arXiv metadata and
  annotated in the .bib: five entries had the right DOI under a wrong title, most wrong DOIs
  were the *adjacent* identifier (constructed by pattern instead of looked up), one paper was
  cited under the wrong first author (Ye et al. 2016, not "Gurnett et al.") with the prose
  corrected to match, and five `author = {others}` entries are filled in from the arXiv API —
  including `fashi_groups`, which the fashienv novelty claim depends on, and `racslow2`,
  published since the entry was written and now carrying its PASA coordinates.
- Two arXiv-assembler bugs that reached a submission payload: `\rightarrow` had no symbol
  mapping so `$V=9.26\rightarrow7.11$` rendered as the nonexistent flux "9.267.11 mJy", and
  `\citealt` (textual) was deleted like parenthetical `\citep`, leaving "(catalogue, )".
- Nine hard-typed prose numbers replaced with the macros recording the same quantity;
  `rfitrend`'s empty stable-line list now renders "none" instead of the repo's
  missing-value marker "--".

### Added
- `dr20radio`: **the spectral index is now measured, not assumed** (`run_alpha`, new
  `results/dr20radio_alpha.json`). The survey-overlap band (-40 < dec <= +30) is covered by
  both VLASS and RACS, so a quasar detected twice gives a two-point alpha directly: 5,571 of
  174,171 band quasars are detected by both. That joint sample is truncation-biased flat
  (median -0.62), because at the RACS limit a steep source has already fallen below the VLASS
  limit at 3 GHz; above S_RACS > 6.2 mJy the truncation cannot operate for any alpha >= -1.5,
  and those 4,190 quasars give **alpha = -0.72 +/- 0.02**. The canonical -0.7 is right.
  The contrast evaluated at the measured index is 4.09% vs 2.82%, a gap of 1.27 pp, moving
  only over 1.25-1.29 across +/-1 sigma -- against 0.23-1.66 for the old 0 -> -1 sweep, which
  is now the range the measurement *excludes* rather than its uncertainty.
- `dr20radio`: `kaplan_meier_median` -- survival analysis replaces the completeness cut as the
  index estimator. A RACS detection with no VLASS counterpart is not missing data: it says
  S_VLASS < the VLASS limit, i.e. alpha is LEFT-CENSORED at
  log(S_lim/S_RACS)/log(nu_V/nu_R), and those are exactly the steep objects a cut discards.
  KM uses all 6,626 RACS-detected band quasars (5,571 measured + 1,055 censored) and is
  unbiased over the whole range instead of over the range a cut assumes: **alpha = -0.755
  +/- 0.012**, against -0.722 for the cut and -0.615 for the naive joint-detection median.
  It lands steeper, as it must, and inside the floor progression the cuts were converging to.
- `dr20radio`: the flux dependence re-tested with the censoring handled. Measured on
  detections only the bins give -0.46/-0.50/-0.52/-0.84; by Kaplan-Meier they give
  -0.61/-0.57/-0.55/-0.90. So up to 0.15 of the apparent trend in the faintest bin was
  truncation -- but the trend **survives**: the brightest sources are genuinely steeper. A
  flux LIMIT is converted at the faint end, where the value is -0.61 (`ALPHA_THRESHOLD_REGIME`),
  giving a gap of 1.11 pp against 1.32 at the sample median; both are reported.
- `dr20radio`: `fetch_racs_total_flux` -- RACS `total_flux_source` over the whole overlap band
  (1.56M sources), so the beam-resolution systematic is **measured rather than estimated**.
  Recomputing with integrated flux on both sides gives alpha = -0.699, a shift of +0.056
  flatward (the 2.5" beam resolves out flux the 25" one keeps). That is a sixth of the
  flux-range term, and replaces a ~0.1 estimate the previous round could only assert.
- `dr20radio`: three sensitivity checks on the measured index, each able to fail, and each
  of which did (`completeness_floor_sensitivity`, `per_epoch`, `flux_bins`). Lowering the
  assumed completeness floor to alpha >= -2 and -2.5 steepens the median monotonically
  (-0.722 -> -0.779 -> -0.856), so -0.722 is an **upper bound on flatness**; single VLASS
  epochs give -0.703/-0.781 against -0.722 for the max-of-epochs (a positively biased
  estimator for a flux ratio); and the median runs -0.50/-0.52/-0.84 across S_RACS bins, so
  it is not transferable across flux. The paper now quotes alpha = -0.722 +/- 0.015 (stat)
  +/- 0.36 (sys) -- the bootstrap SE is the smallest term in the budget by an order of
  magnitude -- and identifies a fourth, unquantified term (peak fluxes across a 2.5" and a
  25" beam, ~0.08 in alpha per 10% flux-ratio error).
- `dr20radio`: `VLASS_S_LIM_CONSERVATIVE_MJY` and `luminosity_matched_vlass_conservative` --
  the MIRROR of the RACS conservative variant, and the side that can actually move the ratio.
  VLASS's 1 mJy is a per-epoch *reliability* threshold while RACS's 3 mJy is a 95%
  *completeness* limit, so at the measured index the north is cut at 1.25 mJy against the
  south's 3.0. Equalising them takes the gap 1.27 -> 1.10 -> 0.91 pp, a 28% reduction: the
  contrast is materially but not wholly a consequence of how the two limits are defined.
- `dr20radio`: `luminosity_matched_per_source_alpha` tests the single-index approximation by
  giving every quasar its own index drawn from the measured distribution. This reports a bias
  rather than a variance on purpose: the realization spread is 0.016 pp and would make any
  breadth of distribution look harmless. The scatter shifts each fraction by ~0.3 pp (to
  3.83%/2.53%) but is nearly common-mode, moving the gap by 0.03 pp. The paper now states
  both -- the contrast survives the scatter, the absolute fractions do not.
- `dr20radio`: spectral-index sweep (`ALPHA_SWEEP = 0, -0.35, -0.7, -1.0`) on the
  luminosity-matched fractions in both legs, committed to `results/dr20radio_{north,south}.json`
  under `luminosity_matched_alpha`. This is the sensitivity test the published 5 mJy
  "robustness check" could not perform: because the common luminosity limit is the RACS one in
  both legs, raising the flux floor rescales north and south identically and leaves their ratio
  unchanged (1.4391 -> 1.4434), whereas alpha moves it 1.08 -> 1.61.
- `make release-check` — verifies a tagged release carries the hand-built papers asset. The
  CI-built papers are synthetic and must never be attached; the upload is deliberately manual.
- `papers/atlas3i/arxiv.yaml` — restores the hand-curated arXiv metadata (astro-ph.IM primary,
  astro-ph.EP cross, page counts, cleaned abstract) that a `make arxiv` run had silently
  overwritten with keyword-inferred guesses.

### Changed
- `dr20radio`: the north/south contrast is now reported as a **range over spectral index**
  (0.23-1.66 percentage points) rather than a single 4.06%/2.82% pair, in the abstract,
  Section 4.3, and the summary. The apparent hemispheric difference is principally an artefact
  of the cross-frequency K-correction, not a measured property of the two populations; at
  alpha = 0 it nearly vanishes. No astrophysical inference is drawn from it.
- `dr20radio`: two hard-typed figures replaced with pipeline-generated macros — the
  out-of-redshift-range census fraction (`\drOutZNorthPct`/`\drOutZSouthPct`) and the northern
  census fraction lying outside the RACS footprint (`\drNorthOutsideRacsPct`).
- `atlas3i`: seven referee findings applied to `main.tex` and `rnaas.tex`, including the
  drift-grid caveat on the EIRP limit; `\aiTotScans` added so the scan count is evidence-backed.

- `innerrc`: the sensitivity scan now excludes variants with a fitted parameter on a bound
  (`bound_contact`, `railed_variants`, `n_fitted`), commits the full per-variant parameter
  vector including `v_halo`/`h_halo`, and reports `chi2_per_n` plus the primary fit's 1-sigma
  rho_DM interval from its own covariance.
- `innerrc`: the anchor now scores the paper's Table-1 solution in chi2 and measures its
  outer-curve (R > 8 kpc) residual, committed as `paper_table1_chi2_per_n` /
  `paper_table1_outer` / `paper_table1_params`.

### Changed
- `innerrc`: **the compatibility-with-consensus claim moved off the variant scan and onto the
  fit's own covariance.** rho_DM = 0.24 with 1-sigma 0.16-0.31 (the scan's old 0.19-0.32 had
  its maximum set by a variant with v_bulge at exactly the 800 km/s bound; 6 of its 8 variants
  are railed, and the 2 interior ones span only 0.20-0.24). The anchor result is also restated:
  their published halo sits +1.33 sigma per point below their own curve beyond 8 kpc over 32
  points, chi2/N 3.18 vs 1.92, so the refit is not merely "another corner of a degeneracy".
- `innerrc`: three overclaims corrected against committed evidence -- the bar-region inner peak
  does *not* reproduce (-36 km/s over 28 bins inside 2 kpc; their 255 km/s peak at 550 pc reads
  214 km/s at 1450 pc), the E/W asymmetry replicates *qualitatively* (period 36% long, damping
  railed, fit seeded at their published answer), and the estimator calibration route is not an
  independent confirmation (algebraically the same statistic on a sub-sample).

### Fixed
- `dr20radio` `refs.bib`: `arnaudova2024` title was the Macfarlane et al. 2021 title. Corrected
  to "Exploring the radio loudness of SDSS quasars with spectral stacking" against Crossref
  (10.1093/mnras/stae233).
- `dr20radio` Limitations now states that the two legs are not matched-sky samples: 14.1% of the
  northern census lies north of the RACS footprint and is unobservable from the south at any depth.

## [1.5.0] - 2026-08-10

### Added
- **`report.preserve_live_macros`** — merges `generated/macros.tex` across run modes instead of
  overwriting it, so a run may only *add* information: a real value always beats a placeholder
  and never the reverse. Wired into `typeii`, `rmstructure`, `torchdsp` and `singlepulse`.
- **`papers/<slice>/arxiv.yaml`** — a tracked, hand-authored override the arXiv assembler merges
  over its auto-extracted metadata. `arxiv-submission/` stays a disposable build artifact, so
  `make arxiv` is idempotent and the human decisions are reviewable in a PR. Written for
  `innerrc`, `peaked`, `frblens`, `vgpra`, `rmdipole`, `pte2`, `rfitrend`, `glitchpop`, `typeii`.
- **`arxiv-submit` resolves `\citet` from `refs.bib`** into natbib textual form, handling
  comma-form and braced compound surnames and `and others`.
- **`rmstructure.N_SYNTHETIC_REALIZATIONS`** and the ensemble macros `\rmsSynRatioEns`,
  `\rmsSynRatioEnsSd`, `\rmsSynNReal` — the honest uncertainty on a recover-a-known over a
  correlated random field.
- **`benchmark_device` / `benchmark_hardware`** in `torchdsp` and `singlepulse` results, so a
  GPU timing run can be recorded without mislabelling CPU-run science as GPU.
- **GPU benchmarks measured** on this workstation's AMD Radeon RX 7600 XT (gfx1102,
  torch 2.12.1+rocm7.1): `torchfdmt` brute-force dedispersion **1.5 s** GPU vs 37.9 s CPU;
  `torchdsp` FFA **0.65 s** vs 7.27 s, SumThreshold **7.63 s** vs 3.04 s — confirming with
  numbers the paper's claim that its per-series loop is GPU-hostile.
- **`tests/test_report.py`** — including the collision case `preserve_live_macros` deliberately
  does *not* rescue, which is the argument for namespacing.
- New blocking validations in the arXiv assembler: unresolvable `\citet`, an abstract starting
  lowercase, macros with no committed value, non-S-parameter and wrong-port-count Touchstone,
  and a self-parse of the generated YAML.

### Changed
- **Every mode-dependent `typeii` macro is namespaced `tiiSyn*`/`tiiReal*`.** `\tiiSource`,
  `\tiiNEvents`, `\tiiCompleteness`, `\tiiPurity` and `\tiiComp*` previously meant different
  things in the two run modes while sharing one name, so merging could not arbitrate and an
  offline rebuild turned `\tiiNEvents` from 768 real observing days into 48 synthetic events —
  under prose reading *"…days, zero failures"*. `papers/typeii/main.tex` cites the namespace it
  actually means.
- **`rmstructure`'s recover-a-known no longer overclaims.** *"recovers an injected plane
  enhancement (4.64 ± 0.35 for an amplitude boost of 5)"* became *"responds to an injected
  low-latitude amplitude boost (3.15 ± 1.11 across 30 field realizations; the
  single-realization bootstrap, 4.64 ± 0.35, understates that scatter threefold)"*. The
  bootstrap resamples within one field realization and so measures sampling noise rather than
  realization variance; the default seed sat 1.3σ high. *"…for an amplitude boost of 5"* is gone
  because the statistic is a band-average over |b| < 10 deg of a profile 5 deg wide and was
  never going to equal the injected peak.
- `make arxiv` packages every paper and fails at the end with the list, instead of stopping at
  the first failure and silently skipping the rest.

### Fixed
- **Four papers' abstracts cited macros that had been blanked to `--`.** Not uncomputed —
  *clobbered*: each slice has two run modes producing different metrics (offline synthetic
  validation vs real census, CPU vs GPU), both write the same `generated/macros.tex`, and each
  emits the other mode's macros as `--`. Whichever ran last silently won. The abstracts cite
  **both** namespaces, so no single run can populate them.
  - `report.preserve_live_macros` merges instead of overwriting: a run may only *add*
    information, a real value always beats a placeholder and never the reverse. Wired into
    `typeii`, `rmstructure`, `torchdsp` and `singlepulse`.
  - **Merging alone was insufficient.** `typeii` left `\tiiSource`, `\tiiNEvents`,
    `\tiiCompleteness`, `\tiiPurity` and `\tiiComp*` un-namespaced, so both modes wrote real
    values and nothing could arbitrate: an offline rebuild turned `\tiiNEvents` from 768 real
    observing days into 48 synthetic events, under prose reading *"…days, zero failures"*.
    Every mode-dependent macro is now `tiiSyn*`/`tiiReal*`, and `papers/typeii/main.tex` cites
    the namespace it actually means. Caught by the science reviewer.
  - Recovered: `typeii` purity **1.0**, completeness **0.917**, curve 0.333/0.625/0.917/1.0 at
    SNR 2/2.5/3/4; `rmstructure` synthetic recovery **4.64 ± 0.35**.
  - `torchfdmt` and `torchdsp` still block — their GPU benchmarks need a ROCm run, which this
    machine cannot do. The merge means a future GPU run will no longer be wiped by a CPU one.
- **Running an offline mode in the repo root destroys the real results JSON.**
  `typeii.run(".", offline=True)` overwrote `results/typeii_metrics.json` with synthetic
  output — 3429 lines deleted, `is_real` True→False, `event_list` gone. Reproduced twice
  (once by me, once by the reviewer). `make guard-real` catches it only at packaging time.
  Documented in `CLAUDE.md`: run offline modes with `out=<tmpdir>`.
- `papers/peaked/arxiv.yaml`: restored the scope-limiting caveat *"a tooling and methodology
  contribution, not a discovery"*, which the first trim dropped. That is a disclaimer against
  over-reading the 6-candidate list, not connective tissue.
- `tests/test_report.py`: unit tests for `preserve_live_macros`, including the collision case
  it deliberately does **not** rescue (which is why namespacing was also required).

- **The two GPU benchmarks were measured**, on this workstation's AMD Radeon RX 7600 XT
  (gfx1102, torch 2.12.1+rocm7.1). `torchfdmt`: brute-force dedispersion **1.5 s** on GPU
  against 37.9 s on CPU. `torchdsp`: FFA **0.65 s** GPU against 7.27 s CPU, and SumThreshold
  **7.63 s** GPU against 3.04 s CPU — confirming, with numbers, the paper's claim that its
  per-series loop is GPU-hostile. **Every paper now packages clean.**
- `torchdsp`: split `benchmark_device` from `device`. One field was labelling both the science
  leg and the timing run, so a real GPU benchmark could only be recorded by mislabelling the
  CPU-run science as GPU. The results JSON now carries `benchmark_hardware` too.
- **`rmstructure` no longer overclaims its recover-a-known.** The abstract said it "recovers an
  injected plane enhancement (4.64 ± 0.35 for an amplitude boost of 5)". Two problems: the
  bootstrap resamples within *one* field realization and so measures sampling noise rather than
  the realization variance of a correlated random field, and the default seed sits high. Across
  30 realizations the recovered ratio is **3.15 ± 1.11** — 3.2× the quoted error, with seed 0's
  4.64 a high outlier. `run(offline=True)` now computes the ensemble
  (`N_SYNTHETIC_REALIZATIONS = 30`) and emits `\rmsSynRatioEns`/`\rmsSynRatioEnsSd`, and the
  paper reports it: *"responds to an injected low-latitude amplitude boost (3.15 ± 1.11 across
  30 field realizations; the single-realization bootstrap, 4.64 ± 0.35, understates that
  scatter threefold)"*. "Recovers … for an amplitude boost of 5" is gone: the statistic is a
  band-average over |b| < 10 deg of a profile 5 deg wide, so it was never going to equal the
  injected peak, and comparing them implied a target the measurement cannot reach.
- **`papers/typeii/refs.bib` `lawrance2024` corrected against Crossref** (DOI
  10.1007/s11207-024-02317-8): authors are Lawrance, Devi, Chandra & Miteva — "Moni-Bidin" was
  not an author — and it is article **75**, not page 58.


## [1.4.0] - 2026-08-07

### Added

- `dr20radio` paper (plan 88): `papers/dr20radio/` (AASTeX, 4 pp) — GATE-2 passed (no
  blockers); should-fixes applied: SDSS Collaboration byline, conservative 5 mJy RACS-limit
  variant (contrast robust: 3.45% vs 2.39%), north carton validation split by selecting
  survey (49% RACS-selected / 27% LOFAR-selected at 3 GHz), overlap-band luminosity-matched
  3.12% surfaced. Slices registered in the Makefile paper targets.
- `dr20radio` increment 3 (plan 88): luminosity-matched north/south contrast — counterpart
  fluxes through the match chain, rest-1.4 GHz common-limit fractions (north 4.06% vs deep
  south 2.82%, overlap 3.12%), `paper_assets` figure + macros from committed evidence.
- `dr20radio` southern leg (plan 88): **the categorical first** — 73,074 DR20 quasars south
  of −40° (pure LCO, confirmed in-data) × RACS-low DR1 → **3.95% radio-detected** at 5″
  (measured false-match 0.06%); overlap band 3.76% vs VLASS's 4.67%; carton validation
  split by selecting survey resolves the pooled 66%: racsradio×RACS 88.9% (the pipeline
  validation) vs lofarradio×RACS 17.2% (cross-frequency fading). Evidence committed.
- `dr20radio` increment 2 code (plan 88): the RACS southern leg — strip-wise resumable
  RACS-low DR1 fetcher (CASDA sync TAP, 1° Dec strips, verified live: 2,123,638 sources),
  `run_south` (deep-south categorical-first census + the −40..+30 overlap band for the
  VLASS cross-check + carton validation against the SELECTING survey), and the two-survey
  synthetic that models the spectral-fading blind spot increment 1 exposed. 9 tests.
- `dr20radio` northern leg (plan 88): 202,691 clean DR20 quasars × VLASS E2/E3 — radio-
  detected fraction 4.67% (any epoch, 2.5″, measured false-match 0.009%), with the
  classic high-z selection rise; the radio-carton VLASS match rate (31–35%) reframed
  honestly as a cross-frequency measurement (selection was at 144/888 MHz — the ~100%
  validation belongs to the RACS leg). E3 interim-list schema handled. Evidence committed.
- `dr20radio` increment 1 (plan 88): census tooling + synthetic recover-a-known — spAll
  quasar selection with the radio-carton circularity exclusion (verified to matter in the
  synthetic round-trip), sky cross-match, position-shift-measured false-match rates
  (density-scaling verified), Wilson-interval detection fractions, and the resumable
  spAll-lite fetcher. 7 offline tests.
- `plans/88-dr20-radio-census.md` — new slice plan (`dr20radio`): the first radio-counterpart
  census of the SDSS-V DR20 Black Hole Mapper quasars (VLASS north / RACS south). GATE 0
  done 2026-08-07: catalog pinned to the file+column level (spAll-lite, 177 MiB, `OBS`
  hemisphere flag), novelty verified (all prior SDSS×radio work is legacy-catalog based; the
  southern SDSS×RACS pairing was categorically impossible before DR20), the eROSITA-rights
  landmine scoped out (X-ray leg deferred), and the radio-targeted open-fiber carton
  circularity identified with its exclusion baked into the plan.

## [1.3.1] — 2026-08-05

### Fixed

- `atlas3i` pre-submission hardening (PRs #169–#171, entries added retroactively): the
  self-citation re-pinned to v1.3.0; pre-flight review of the arXiv package caught and fixed
  a missing figure (`main.tex` referenced `fig:survey` with no figure environment — the PDF
  printed "Figure ??"), a 404ing Choza et al. DOI, the since-minted Jacobson-Bell RNAAS DOI,
  and the FAST paper's DOI/initials; final reader-pass fixes: the cadence duration corrected
  to the paper's 30 min, the EIRP-agreement claim restated with its exact substitutions
  (δν_t = δν, β = 1) instead of "algebraically identical", and the figure title's literal
  backslash removed at the source.

## [1.3.0] — 2026-08-05

### Added

- `innerrc` GATE-2 PASS with fixes: Sofue & Kohno DOI added (10.1093/pasj/psaf114); the
  pre-calibration Table-2 offset made traceable (`median_dv_kms_fixed_sigma15` = −14.5,
  superseding a pre-Q4-fix −13.4); bulge claim scoped and macro-fed; sensitivity variants now
  commit bulge/disc params; E/W amplitude compared to their stated ≃15 km/s mid-disc figure;
  MG&D16 apples-to-apples clause; hidden-constraint and Davis+2025 phrasing softened.
- `innerrc` increment 3 (plan 86): the full-length paper (`papers/innerrc/main.tex` +
  `refs.bib`, compiles in the tectonic container) — anchor degeneracy + raw-HI4PI
  replication + estimator head-to-head, all numbers from `paper_macros` reading the two
  committed evidence JSONs (the committed-real-results pattern end-to-end).
- `innerrc` increment 2 (plan 86): the HI4PI real-data leg — 1,113 sightlines through both
  TVM estimators close the `hi` slice's caveat with a measured number (threshold reads
  +17.9 km/s high; per-estimator σ calibration gives 9.5 vs 26.6 km/s, bracketing the
  paper's adopted 15); their Table 2 reproduces at the ~4% level from raw survey data; the
  E/W asymmetry replicates qualitatively (period/phase match, softer amplitude). Fourth-
  quadrant |sin ℓ| bug and a CDS stall (now timeout+resume+retry) found and fixed. Five
  paper-style figures + evidence JSON committed.
- `innerrc` increment 1 (plan 86): tooling + synthetic recover-a-known + the **offline anchor
  result** — decomposing Sofue & Kohno 2025's own published RC tables (vendored from the arXiv
  source) validates their ρ_DM = 0.107 GeV/cm³ arithmetic exactly (through the eq.-25 4π
  convention, now covered by a test) while showing it is one corner of a broad disc–halo
  degeneracy: an unconstrained refit of the same curve fits marginally better at
  ρ_DM = 0.24, and all eight sensitivity-scan variants land in 0.19–0.32 — the consensus
  range. Gaussian-decomposition TVM (greedy + joint refinement) recovers injected terminal
  velocities on crowded synthetic spectra and measures the `hi` threshold estimator's
  documented high bias on the same spectra. Committed evidence: `results/innerrc_anchor.json`.

### Changed

- **Committed-real-results migration (the structural fix from the 2026-07-31 incident):**
  real-run outputs — `results/*.json|csv`, every paper's `figures/` and `generated/` macros —
  are now git-tracked evidence, reviewed in PRs, replacing the "regenerable artifacts are
  gitignored" policy that made the synthetic-clobber possible. The offline Snakemake DAG is
  demoted to CI smoke only; `make guard-real` continues to gate packaging. All artifacts
  committed in this change are guard-verified real (or by-design-disclosed synthetic).

## [1.2.0] — 2026-08-04

### Fixed

- **Integrity: packaged papers can no longer mislabel synthetic output as real data.** Caught
  2026-07-31 via the released `hi` paper, which quoted the offline synthetic-fixture rotation
  curve (231 km/s, ±0.3 scatter) under a caption claiming LAB data — `make papers-zip` depended
  on the offline `figures` DAG, which regenerates every slice's figures/macros from synthetic
  fixtures into the same files the papers `\input`. Now: `papers-zip` depends on a new
  `make guard-real` (`scripts/guard_real_results.py`) that fails the build if any
  `results/*.json` is synthetic-sourced (allowlist requires per-entry justification, currently
  empty); the CI release workflow no longer attaches its offline build to Releases (renamed to
  a clearly-labelled `papers-synthetic-smoke` artifact with a MANIFEST warning) — distributable
  zips are built locally from real-sourced results and uploaded manually. The v1.1.0 release
  asset's `hi/main.pdf` was regenerated from the real LAB leg (257 ± 6 km/s) and replaced, with
  the correction noted in the zip's MANIFEST.

### Added

- `atlas3i` full 1–12 GHz completion: S/C/X sweeps done (60 unique node-cadences total,
  1,938,860 raw hits → 294 on/off survivors → **0 confirmed** — the paper's full-band null
  reproduces). `survey_summary`/`survey_macros`/`survey_figure` aggregate all bands for the
  papers; satellite table extended (S-MSS, S-DARS, Ku FSS/DBS — X's only 4 survivors are
  Ku-band TV tones above the analysed passband); findings updated with the full funnel and
  measured cost (~95 h compute, ~3.7 TB transfer, one desktop); RNAAS updated to the
  full-band claim and a full-length `main.tex` added for arXiv.
- `plans/86-inner-rotation-curve.md` — new slice plan (`innerrc`): independent replication of
  Sofue & Kohno 2025 (PASJ 77, 1335; arXiv:2509.23581), the modern inner Milky Way rotation
  curve. GATE 0 done: full-text read (their 0.107 GeV/cm³ local DM density is halo-only and
  author-framed as a lower limit — the slice maps its sensitivity, not a refutation), novelty
  pass (cited, never audited), data verified (HI4PI CDS tiles ~15–20 GB for |b|≈0; their RC
  ASCII tables in the arXiv source as the anchor). Also retires the `hi` slice's documented
  2 K-threshold bias via a Gaussian-decomposition TVM head-to-head. First plan written under
  the committed-real-results integrity rule.
- `atlas3i`: S/C/X band support — the pinned node→band tables for all four receivers
  (`BAND_NODES`, with the six pairwise-duplicated bank-boundary recordings per C/X band
  deduplicated and documented in `DUPLICATE_NODES`), and `--sweep`/`sweep_summary`/
  `sweep_figure` generalised over `--band`. Extends the L-band reproduction toward the
  paper's full 1–12 GHz claim (S: 6 nodes; C: 23; X: 25 unique nodes).

### Fixed

- CI paper builds: `atlas3i` registered in the static-slice DAG (`workflow/Snakefile`, via its
  offline `--paper` entry point over the committed result JSONs), so `make figures`/CI
  regenerate `papers/atlas3i/generated/macros.tex` + the sweep figure before tectonic runs —
  the v1.1.0 `paper`/`release` workflows failed on the missing gitignored macros file.

## [1.1.0] — 2026-07-30

### Added

- `plans/85-bl3i-atlas-drift-reproduction.md` — new slice plan: independent reproduction of the
  Breakthrough Listen 3I/ATLAS GBT nondetection (RNAAS arXiv:2512.19763) from the public
  `bldata.berkeley.edu/ATLAS` archive, reusing the plan-11 `driftsearch` machinery. GATE 0 done
  in-plan (data verified public with sizes/layout; novelty pass: no independent reanalysis of
  the released data exists — ATA/FAST/MeerKAT follow-ups are all original observations).
- `src/jansky_research/atlas3i.py` + `tests/test_atlas3i.py` — the plan-85 tool: archive-index
  parser and ABACAD cadence selector (validated against a vendored real GB_ATLAS listing),
  bandpass normalisation with DC-spike excision (the artifact that defeated plan 11's teaching
  detector on real BL data), a physical-units (Hz/s) brute-force de-Doppler search with robust
  per-drift S/N and neighbour suppression, the ABACAD on/off sky-localisation filter, a
  narrowband EIRP limit that reproduces the paper's ~100 mW headline from GBT L-band
  parameters, and a fully offline synthetic-cadence round-trip (`run`) that recovers an
  injected drifting tone while rejecting always-on RFI and the DC spike. Real-data leg
  (`--node`, network + `voyager` extra) streams each 10 GB fine-resolution file in bounded
  frequency chunks and deletes scans after searching (peak disk ≈ one cadence). The node→band
  mapping (L: blc21–26 = 939–2064 MHz, etc.) was pinned by remote HDF5 header reads and is
  embedded as constants. Candidate vetting: per-survivor drift-coherence stamps
  (`vet_stamps`, inline during the search) plus a satellite-allocation exclusion axis
  (`classify_band`; Iridium/Inmarsat/GNSS downlinks defeat two-position filters by design).
- `papers/atlas3i/` — the RNAAS note (`rnaas.tex` + `refs.bib`, compiles in the tectonic
  container): the reproduction result + the filter-evasion taxonomy, with every number
  `\input` from `sweep_macros`/`sweep_figure` (new module functions reading the committed
  per-node result JSONs — nothing typed by hand). README slice tables updated. GATE-2
  review passed with fixes (threshold asymmetry implemented + honest caveats: raw-hit-count
  comparability, analysed-band scope, β=1, notch filters); distance Horizons-pinned to
  1.79801 au → final EIRP limit 99.2 mW.
- `survey/atlas3i-findings.md` — the L-band real-data result: all six nodes (939–2064 MHz)
  searched and vetted, 261 ABACAD survivors, 0 confirmed — **the paper's L-band null
  reproduces independently**, with a matching ~99 mW EIRP limit; includes the two instructive
  filter-evasion modes (sub-threshold OFF carriers, satellite downlinks) and honest caveats
  (L band only, low-drift blind spot, S/N-convention and SEFD/distance nominals).

## [1.0.1] — 2026-07-25

### Changed

- **Install: a single clone now works.** The `jansky` dependency defaults to its pinned git tag
  (`v0.2.0`) in `[tool.uv.sources]` instead of the `../jansky` path source, so `git clone && uv sync`
  succeeds with no sibling checkout — removing a JOSS review-checklist blocker and the most likely
  first-run failure for new users. Cross-repo development is now an explicit opt-in:
  `eval "$(make -s dev-env)"` puts a sibling `../jansky/src` on `PYTHONPATH` ahead of the pinned tag
  (`uv run` preserves `PYTHONPATH`; a `uv pip install -e` overlay does *not* survive, because
  `uv run` re-syncs the environment). CI/paper/release workflows drop their second
  `actions/checkout`, so CI installs exactly what a stranger gets.
  Docs updated: `README.md`, `CONTRIBUTING.md`, `REPRODUCING.md`, `docs/usage.md`,
  `joss/SUBMISSION.md`, `CLAUDE.md`.
- Documentation: added the Zenodo DOI badge (concept DOI `10.5281/zenodo.21482378`) to the README
  following the v1.0.0 release.
- Documentation: added community guidelines for JOSS review — `CONTRIBUTING.md` (how to contribute,
  report issues, and get support; dev setup; the slice pattern; the pre-PR checks),
  `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), and `docs/usage.md` (a short install-and-run
  guide). Linked from the README.
- JOSS paper (`joss/paper.md`): added a substantial-effort paragraph to the Statement of need
  (breadth across domains + the dual Snakemake/Airflow reproducibility layer + "the software is
  itself the research"). Refreshed the stale "twenty-six slices" tally → "forty slices plus a
  synthesis" in the README Results header and `CLAUDE.md` (the slice table has 41 rows).
- Papers now cite the toolkit: the `vgpra` and `spectra` RNAAS notes `\software{}`-cite
  **`jansky-research`** via its Zenodo concept DOI `10.5281/zenodo.21482378` (new `@misc{janskyresearch}`
  `refs.bib` entry). Added a related-work disclosure to `joss/paper.md` (the two notes are in
  preparation for RNAAS, distinct from the software paper), and recorded the "cite jansky-research
  going forward" convention in `CLAUDE.md` and the `research-publish` skill.
- `research-publish` skill: the readiness check now **auto-discovers every `papers/*/rnaas.tex`**
  note (was hardcoded to `frbstats`) with a per-note 1000-word sanity check, so `vgpra`/`spectra`
  and any future note are covered without editing the script.
- Added `joss/SUBMISSION.md` — a paste-ready helper for the manual JOSS submission (form fields,
  the comments-to-editor text with the substantial-effort + related-work + AI-use disclosures, a
  verified suggested-reviewers shortlist, and RNAAS Editorial-Manager cover text).
- Added `docs/faq.md` — how others use the toolkit (depend, don't fork), how the in-repo papers work
  (authorship, provenance/priority via Zenodo DOIs, what "unpublished" does and doesn't mean), and a
  licensing note. Linked from the README and `CONTRIBUTING.md`.
- **Dual-licensed the repository:** the papers in `papers/` are now **CC BY 4.0** (new
  `papers/LICENSE`); the code remains **MIT**. README, `docs/faq.md` updated to state it.

## [1.0.0] — 2026-07-21

**Initial public release.** With no prior tag this records the full toolkit as it stands rather
than a diff from an earlier version; every later section is a diff from its predecessor.

### Added

- The `jansky-research` toolkit — ~40 self-contained, reproducible research slices (a tested tool
  → real public data → adversarial science-review gate → honest AASTeX write-up), grouped by
  domain:
  - *FRB & time-domain* — `frbstats`, `frbperiod`, `frbwait`, `frblens`, `singlepulse`;
  - *Pulsars* — `pulsarspec`, `ppdot`, `pte2`, `glitchpop`, `wdpulsar`;
  - *HI & spectral line* — `hi`, `fashienv`;
  - *Solar & heliospheric* — `solarbursts`, `windwaves`, `swaves`, `triangulate`,
    `type3synthesis`, `typeii`, `ecallisto_census`, `ecallisto_catalog`, `rfitrend`;
  - *Planetary radio* — `junodam`, `skr`, `vgpra`;
  - *RM / Faraday & cosmology* — `rmsky`, `rmstructure`, `rmdipole`;
  - *Continuum & variability* — `vlass`, `vlbi`, `stacking`, `sourcecounts`, `peaked`,
    `southern`, `spectra`, `offsets`, `stokesv`, `stokesv_discovery`, `lpt`, `lptv`;
  - *SETI* — `driftsearch`;
  - *GPU / ML* — `fdmt` + `torchdsp` (device-portable pure-PyTorch DSP), `svsbi`
    (neural simulation-based inference).
- The shared slice scaffolding in `src/jansky_research/`: `data.py`, `pipeline.py`, `report.py`.
- Reproducibility paths: `make reproduce`, the Snakemake static-slice file-DAG
  (`workflow/Snakefile`, `make figures`), and the Airflow streaming e-Callisto ingest
  (`airflow/`).
- Publishing / data helpers under `.claude/skills/`: `arxiv-submit`, `casda-cutout-fetch`,
  `pull-station-data` (the `jansky-observe` rooftop-station bridge).
- **Release-versioning infrastructure** — `VERSIONING.md` (SemVer policy for this repo),
  this `CHANGELOG.md`, and `scripts/next_version.py` (recommends the next version from the
  `Unreleased` section, with its reasoning).
- `papers/vgpra/rnaas.tex` — an RNAAS short-form of the Voyager 2 PRA ice-giant rotation-period
  reanalysis (the recover-a-known → controlled-null showcase).
- `papers/spectra/rnaas.tex` — an RNAAS short-form showing raw TGSS×NVSS ultra-steep-spectrum
  selection is dominated by the TGSS flux-scale systematic (the "apparent signal is a
  systematic" cautionary note).

### Changed

- Optional GPU acceleration (`fdmt`, `sbi` extras, and the `torchdsp` slice) is pure-PyTorch and
  ROCm/CUDA-portable; the core install and CI remain CPU-only, GPU is opt-in.
- Refreshed the JOSS paper (`joss/paper.md` + `paper.bib`), `CITATION.cff`, and `.zenodo.json` to
  the current scope: "CPU-first with optional GPU" (was the false "CPU-only"), a domain-grouped
  capability list over the full >40-slice toolkit (was a stale six-module snapshot), and a
  Statement of need reframed around recover-a-known + honest-null at scale.

[Unreleased]: https://github.com/joebarbere/jansky-research/compare/v1.3.1...HEAD
[1.3.1]: https://github.com/joebarbere/jansky-research/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/joebarbere/jansky-research/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/joebarbere/jansky-research/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/joebarbere/jansky-research/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/joebarbere/jansky-research/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/joebarbere/jansky-research/releases/tag/v1.0.0
