---
name: research-publish
description: Guide a research repo's papers + software out the door — Zenodo (DOI), JOSS (software paper), RNAAS (note), and arXiv (preprints), in dependency order. Runs a static readiness check on the repo's artifacts and walks through the manual, account-bound steps. Use when the user wants to publish or asks "what's left to submit".
---

# Research publishing helper

**Reality first (state this to the user):** none of these submissions can be done by the repo — each
needs the user signed in to an external service, and DOIs/endorsements are issued by people, not code.
What this skill *does*: run a **static readiness check** on the repo's artifacts (no network, no
accounts), report what's present/missing, and walk the venues in the right **dependency order**. It
complements the per-paper `arxiv-submit` skill (which assembles one arXiv package) by covering the
whole Zenodo → JOSS → RNAAS → arXiv flow. The canonical account-bound checklist lives in Joe's
personal notes, outside this repo (Obsidian vault: `efforts/radio_astronomy/research_paper_todo.md`).

## When invoked

1. **Run the readiness check:**
   `uv run python .claude/skills/research-publish/check_publish_readiness.py`
   It inspects `LICENSE`, `CITATION.cff`, `.zenodo.json`, `joss/paper.md`+`paper.bib`, every
   auto-discovered `papers/*/rnaas.tex` note, and each fresh-angle paper's `arxiv-submission/`
   package, printing `[x]`/`[ ]`/`[~]` with the next action per venue.
2. **Consult the checklist** (`efforts/radio_astronomy/research_paper_todo.md` in the Obsidian
   vault, if reachable; otherwise the order below) for the full manual walkthrough, then guide the
   user venue by venue.

## The order (and why)

1. **Zenodo first** — the only hard dependency: it mints the permanent DOI **JOSS requires** at
   acceptance, and it's the fastest (connect the repo, cut a GitHub release, copy the concept DOI).
   Then add the DOI badge to the README.
2. **JOSS** — reviews the *software*, not novelty; needs the step-1 DOI. Lead the submission notes
   with the validated analyses + the Airflow/Podman reproducibility layer (clears the "substantial
   scholarly effort" bar).
3. **RNAAS** — independent of 1–2; quickest (editorial screening, no peer review). Confirm the current
   AAS fee; submit via Editorial Manager. Not posted to arXiv.
4. **arXiv** — the fresh-angle papers only (`frbstats`, `vlass`, `peaked`, `southern`). Register and
   **line up an endorser early** (multi-day lead time on a first `astro-ph` submission). Run
   `make arxiv` to assemble each `papers/<slice>/arxiv-submission/`, fill the `metadata.yaml` TODOs
   (`primary_category`, `comments`), then upload the `arxiv-source.tar.gz`. Among them, submit `vlass`
   or `peaked` first. The pure reproductions/negatives stay repo + Zenodo only.

## Notes

- Author is **Joseph Barbere**, ORCID `0009-0008-3289-4447`; an AI/LLM is **not** an author (it's in
  each paper's AI-use disclosure + `\software{}`). Keep the name consistent on every account.
- **Cite the toolkit.** Since the v1.0.0 Zenodo archive, every new slice paper / RNAAS note should
  `\software{}`-cite **`jansky-research`** (the toolkit that produced the result) alongside `jansky`,
  via its Zenodo **concept DOI `10.5281/zenodo.21482378`** (all-versions). Add this `refs.bib` entry
  (already in `papers/vgpra/` and `papers/spectra/`):
  ```bibtex
  @misc{janskyresearch,
    author       = {Barbere, Joseph},
    title        = {jansky-research: A CPU-first, reproducible toolkit for amateur radio-astronomy analyses},
    year         = {2026}, version = {1.0.0}, doi = {10.5281/zenodo.21482378},
    howpublished = {\url{https://doi.org/10.5281/zenodo.21482378}}
  }
  ```
- **JOSS co-publication:** when submitting the JOSS software paper, **disclose related publications**
  (published / in review / nearing submission) — currently the `vgpra` + `spectra` RNAAS notes. The
  software clears JOSS's bar on its own (40+ slices + the dual reproducibility layer), so a parallel
  RNAAS submission is fine and is not a reason to delay JOSS.
- The check is artifact-only: it can't verify a DOI resolves or an endorser is lined up — those are
  the `[~]`/manual items. Re-run it after each milestone to see the chain advance.
