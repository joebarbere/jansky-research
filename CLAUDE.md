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

**Merged slices:** twenty-six plus the type III synthesis — the authoritative per-slice table
(tool + outcome) is in `README.md`. Each has a paper under `papers/<slice>/`. Publishing steps
are in `TODO.md` (Zenodo → JOSS → RNAAS → arXiv).

## Active direction (2026-07)

- **Pick the next slice from `fable-ideas.md`** (2026-07-05, a 12-agent deep re-scan; supersedes
  the shortlist in `survey/opportunity-scan-2026-07.md`, whose Tier-1 items are now merged:
  `stokesv_discovery`, `lpt`, `rmstructure`, `torchfdmt`, `junodam`). Suggested first moves
  there: F6 (`torch-dsp` GPU suite), F4 (WD-pulsar sweep), F8 (FASHI DR2). Executed so far:
  F1 → `rmdipole` (plan 38); F2+F5 → `frbwait`+`frblens` (plans 39+42; Cat 2 mirrored from
  CANFAR DOI 10.11570/25.0066 — chime-frb.ca itself is still 503). All fable-ideas have plans
  (`plans/38`–`84`).
- **Standing GATE-0 for anything from that file:** the scan session couldn't fetch primary
  sources (egress-blocked), so do a full-text novelty pass + a data-URL check before writing the
  plan. Its "Corrections & closed doors" section lists ideas already killed — check it first.
- Earlier blockers are resolved: CASDA recovered (`stokesv` complete, both legs; the discovery
  census is merged as `stokesv_discovery`); the southern GLEAM-X×RACS curvature catalogue is
  merged (`southern`).

## Layout

`src/jansky_research/` (slice modules + `data.py`/`pipeline.py`/`report.py`) · `tests/` · `plans/` ·
`fable-ideas.md` (current plan-ready idea list) · `survey/` (committed findings) ·
`papers/<slice>/` · `station/` (build guides for the physical rooftop
station — self-collected data, WIP; owner's working notes live in an Obsidian vault, not here) ·
`workflow/Snakefile` (static-slice file-DAG, drives `make figures`) · `airflow/` (streaming e-Callisto
ingest, Podman DAG) · `.claude/` (agents + skills) ·
`Makefile` (`setup`/`test`/`cov`/`lint`/`typecheck`/`figures`/`paper`/`arxiv`/`reproduce`).
