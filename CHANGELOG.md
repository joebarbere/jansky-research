# Changelog

All notable changes to `jansky-research` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/) as codified in [`VERSIONING.md`](VERSIONING.md).

Every PR adds an entry to `## [Unreleased]`. `scripts/next_version.py` reads that section to
recommend the next version number.


## [Unreleased]

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
