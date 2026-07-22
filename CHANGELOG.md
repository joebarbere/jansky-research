# Changelog

All notable changes to `jansky-research` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/) as codified in [`VERSIONING.md`](VERSIONING.md).

Every PR adds an entry to `## [Unreleased]`. `scripts/next_version.py` reads that section to
recommend the next version number.

## [Unreleased]

### Changed

- Documentation: added the Zenodo DOI badge (concept DOI `10.5281/zenodo.21482378`) to the README
  following the v1.0.0 release.
- Documentation: added community guidelines for JOSS review — `CONTRIBUTING.md` (how to contribute,
  report issues, and get support; dev setup; the slice pattern; the pre-PR checks),
  `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), and `docs/usage.md` (a short install-and-run
  guide, incl. the sibling-checkout install a fresh clone needs). Linked from the README.
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

[Unreleased]: https://github.com/joebarbere/jansky-research/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/joebarbere/jansky-research/releases/tag/v1.0.0
