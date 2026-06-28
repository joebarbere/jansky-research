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
whole Zenodo → JOSS → RNAAS → arXiv flow. The repo's `TODO.md` is the canonical account-bound checklist.

## When invoked

1. **Run the readiness check:**
   `uv run python .claude/skills/research-publish/check_publish_readiness.py`
   It inspects `LICENSE`, `CITATION.cff`, `.zenodo.json`, `joss/paper.md`+`paper.bib`, the RNAAS note,
   and each fresh-angle paper's `arxiv-submission/` package, printing `[x]`/`[ ]`/`[~]` with the next
   action per venue.
2. **Read `TODO.md`** for the full manual walkthrough, then guide the user venue by venue.

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
- The check is artifact-only: it can't verify a DOI resolves or an endorser is lined up — those are
  the `[~]`/manual items. Re-run it after each milestone to see the chain advance.
