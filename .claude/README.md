# Claude Code automation for jansky-research

This directory holds the project's **agents** and **skills**. See the repo `CLAUDE.md` for the
working rules and the relationship to the sibling `jansky` course (`../jansky`).

## Agents (`agents/`)

Data-analysis agents added for this repo:

- `dataset-analyst` — drives `jansky_research.pipeline` over a fetched public dataset.
- `pipeline-runner` — operates the local Airflow-on-Podman stack.
- `results-interpreter` — writes the honest results interpretation, flags overclaiming.

Plus copies of the `jansky` course's reusable agents (kept here so they work without `../jansky`):
`archive-scout` (data discovery), `radio-research-assistant` (literature synthesis),
`science-reviewer` (the GATE-1 / GATE-2 / paper-correctness gates).

## Skills (`skills/`)

- `arxiv-submit` — assemble + validate an arXiv submission package per paper (`make arxiv` uses it).
- `casda-cutout-fetch` — best-effort Playwright download of RACS Stokes-V FITS cutouts from the CASDA
  web Cutout Service when the VO APIs misbehave (prefer `astroquery.casda`; see its `SKILL.md`).
- `find-radio-papers` — literature search (ADS/arXiv), ported from the course and adapted to cite
  this repo's `survey/`/`refs.bib` first.
- `radio-source-lookup` — source/position profile (SIMBAD/NED + VizieR), ported and adapted to reuse
  this repo's `spectra`/`stokesv` fetchers and `spectra.spectral_index`.

## Skill sync with the `jansky` course

Skills are discovered only from *this* repo's `.claude/skills/` (a sibling repo's skills are not
auto-available), so the two general helpers above are **copied** here, not symlinked — keep them
roughly in step with `../jansky/.claude/skills/` when either changes. The course additionally has
`dataset-watch` and `radio-mastodon`, which are course-flavoured (they read `docs/mastodon.md` etc.)
and are intentionally **not** ported. This repo's `arxiv-submit` / `casda-cutout-fetch` are
research-specific and intentionally **not** in the course.
