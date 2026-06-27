# Plan 07 — The amateur research paper 📋

> Context: the headline deliverable. Depends on 06 (GATE 2). Scope: medium.

## Context

Present the whole slice — survey, new code, data analysis, and findings — as an amateur research
paper in PDF form, authored by **Joe Barbere** and **Claude**. It is framed as a
**reproducibility / tooling contribution**, not a discovery claim, and every number is auditable.

## Deliverables

- `paper/main.tex` — **AASTeX 6.3** (`aastex631`), arXiv-acceptable.
- `paper/sections/*.tex` — intro, survey, methods (the new tool), data & analysis, results,
  discussion, reproducibility, conclusions.
- `paper/refs.bib` — seeded from jansky's `docs/references.md` + `docs/papers-timeline.md` (ADS
  bibcodes/DOIs → BibTeX), plus the survey's verified citations; every entry resolves live.
- `containers/paper.Dockerfile` — **tectonic** build (hermetic; pre-warmed package cache for
  offline); `make paper` → `paper/main.pdf`.
- `.github/workflows/paper.yml` — builds the PDF, uploads the artifact.

## Approach

`main.tex` `\input`s `paper/generated/macros.tex` and `\includegraphics` the generated figures —
**no number or figure is hand-produced**, so the paper regenerates from the pipeline.
**science-reviewer** gates correctness and overclaiming; **radio-research-assistant** verifies
every citation.

## Verification

- `make paper` builds `paper/main.pdf` in the tectonic container; authors render as **Joe Barbere**
  and **Claude**.
- Every headline number in the PDF traces to `paper/generated/macros.tex`; every `refs.bib` entry
  resolves; science-reviewer signs off on claims and citations.
