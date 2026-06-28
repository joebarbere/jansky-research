# jansky-research — guide for Claude

**What this is.** Amateur radio-astronomy *research*, end to end. A public sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course (checked out at `../jansky`): where
jansky *teaches*, this repo *does original work* as reproducible "slices" — one gap → one tested tool
→ real public data → honest write-up. It **depends on `jansky` as a library** (`from jansky import …`)
and mirrors its `uv`/ruff/mypy/pytest conventions.

## The two repos (each should point at the other)

- **`../jansky`** — the course (library we depend on). Its `.claude/skills/` hold the general
  research helpers: `find-radio-papers`, `radio-source-lookup`, `dataset-watch`, `radio-mastodon`.
- **this repo** — the research. Its `.claude/skills/` hold the publishing/data helpers: `arxiv-submit`,
  `casda-cutout-fetch`. (Ports of `find-radio-papers` + `radio-source-lookup` live here too so they
  work without the course checked out.)

A new session in *either* repo: read this file (or jansky's `CLAUDE.md`) to learn the other exists.

## Working rules (non-negotiable)

- **Branch before committing — never commit on `main`.** Squash-merge PRs; delete the branch.
- **85% coverage floor**, `ruff` (line-length 100) + `mypy` clean, before every PR.
- Generated artifacts are **gitignored** (`papers/*/figures/*`, `generated/*`, `results/*`); only
  `main.tex`/`refs.bib`/`.gitkeep` are tracked per paper. Macros come from the pipeline.
- **Honest framing**: validations, limits, and negatives reported plainly; no overclaiming. The
  `science-reviewer` agent gates this.
- Commit footer: `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
  PR footer: `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

## The slice pattern (how every result is built)

tested helper (pure NumPy/SciPy/astropy + synthetic offline fixture) → real-data run (network,
`# pragma: no cover`) → GATE-2 science review → AASTeX paper (`papers/<slice>/`) → arXiv package
(`make arxiv`). Plan in `plans/NN-*.md`; findings in `survey/<slice>-findings.md`.

**Merged slices:** frbstats, spectra (USS), frbperiod, driftsearch, hi, vlass, peaked. Each has a
paper under `papers/<slice>/`. Publishing steps are in `TODO.md` (Zenodo → JOSS → RNAAS → arXiv).

## Active direction (2026-06)

- **`stokesv`** — RACS Stokes-V coherent-emitter slice (`plans/15`, `src/jansky_research/stokesv.py`).
  Tooling + the SRSC recover-a-known + forced-photometry core (`measure_circular_pol`) are done and
  tested. The live forced-photometry leg is **blocked on a CASDA outage** (VO TAP/SIA2 *and* the web
  cutout service return 0/errors for every position; OPAL login works). See `survey/stokesv-findings.md`
  and the `casda-cutout-fetch` skill. Retry CASDA; when it recovers, finish the leg.
- **Queued runner-up:** southern GLEAM-X×RACS spectral-curvature catalogue — reuses `peaked`/`spectra`,
  fixes peaked's TGSS-upper-limit limitation with real GLEAM-X in-band turnovers; image access
  de-risked (Data Central serves GLEAM-X Stokes-I, VizieR serves the catalogues). See
  `survey/new-findings-scan.md`.

## Layout

`src/jansky_research/` (slice modules + `data.py`/`pipeline.py`/`report.py`) · `tests/` · `plans/` ·
`survey/` (committed findings) · `papers/<slice>/` · `airflow/` (Podman DAG) · `.claude/` (agents +
skills) · `Makefile` (`setup`/`test`/`cov`/`lint`/`typecheck`/`figures`/`paper`/`arxiv`/`reproduce`).
