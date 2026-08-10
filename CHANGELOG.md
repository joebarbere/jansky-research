# Changelog

All notable changes to `jansky-research` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/) as codified in [`VERSIONING.md`](VERSIONING.md).

Every PR adds an entry to `## [Unreleased]`. `scripts/next_version.py` reads that section to
recommend the next version number.

## [Unreleased]

### Fixed

- `arxiv-submit`: the assembler silently mangled abstracts containing textual citations or
  Greek letters, and the damage read as finished prose. `\citet{key}` was deleted outright,
  so an abstract opening "\citet{sofue2025} derived..." became "derived the definitive
  modern..." — a sentence with no subject. `\rho` and `\pi` were absent from the symbol
  table, so `$\rho_\mathrm{DM} = 0.107$` extracted as " = 0.107" and `$4\pi$` as "4".
  `\citet` now leaves a visible `[CITE:key]` marker, `validate()` treats that marker (and a
  lowercase first word) as **blocking errors**, and the Greek set is extended. Found while
  preparing the `innerrc` package; re-running the assembler across the other papers shows
  `spectra` and `vgpra` had the same latent defect.
- `arxiv-submit`: **the generated `metadata.yaml` required hand-editing and was overwritten
  by every `make arxiv`** — a footgun that silently destroyed the fix for the problem above.
  Hand-authored values now live in a tracked `papers/<slice>/arxiv.yaml` which the assembler
  merges over the auto-extracted ones, so regeneration is idempotent and the decision is
  reviewable in a PR.
- `arxiv-submit`: **`\citet{key}` is now resolved from `refs.bib`** to its natbib textual form
  ("de Gasperin et al. (2018)", "Zhu & Zheng (2025)"), handling comma-form and braced compound
  surnames and `and others`. This fixes the common case with no override at all — 9 papers.
- `arxiv-submit`: **macros whose value contained braces were never loaded**, because the
  `\newcommand` regex used `[^{}]*`. Any value with an exponent, subscript or `\mathrm{}` was
  dropped, and the abstract's `\rfRealHUMAINP` then fell through the generic macro strip and
  vanished — leaving prose reading `p=,`: a stated significance with the number silently gone.
  Now brace-matched. Affected `dr20radio`, `frbperiod`, `frbstats`, `glitchpop`, `offsets`,
  `rfitrend`.
- `arxiv-submit`: `_apply_macros` passed macro values to `re.sub` as **replacement templates**,
  so any value containing a backslash raised `bad escape` (or silently expanded a group
  reference). Only surfaced once brace-nested values started loading. Now a literal callable.
- `arxiv-submit`: the generator emitted **invalid YAML** for any multi-line override (only the
  abstract's first line was indented) and now parses its own output as a validation step.
- `arxiv-submit`: **new blocking check for macros with no committed value.** A metric that is
  `null` in the results JSON renders as `--`, which in a sentence reads as a hole — *"the
  detector separates type II from type III and RFI at purity --"* — while the prose around it
  still parses as a claim. Flags `typeii` (`\tiiPurity` + three completeness macros),
  `torchfdmt` (`\spBruteGpu`), `torchdsp` (two GPU benchmarks), `rmstructure` (an injected
  recovery and its error). **These are paper problems, not packaging problems**, and are left
  flagged deliberately: each needs the metric computed or the sentence rewritten.
- `make arxiv` **stopped at the first paper with a blocking error and silently skipped every
  paper after it**, so a partial run looked like a complete one. It now packages all of them
  and fails at the end with the list.
- Abstracts over arXiv's ~1920-char limit trimmed into tracked `arxiv.yaml` overrides for
  `frblens`, `glitchpop`, `peaked`, `pte2`, `rfitrend`, `rmdipole` and `vgpra` (and `innerrc`
  for its citation). Trims cut connective tissue and closing "…are the contribution" summaries
  only; a mechanical guard rejected any trim that dropped a numeric result.

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
