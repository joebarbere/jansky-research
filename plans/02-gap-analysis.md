# Plan 02 — Gap analysis & domain pick (GATE 1) 📋

> Context: the convergence gate. Pairs with 01 (survey) and 03 (tooling). Scope: small.

## Context

The survey produces several plausible gaps across five domains. To keep the work a single,
honest vertical slice, exactly **one** must be chosen — with one openly-downloadable dataset and
one `jansky` helper to build on. This plan ranks the gaps and records the decision. **No tooling
is written until this gate passes.**

## Deliverables

- `survey/gap-analysis.md` — every shortlisted gap scored on:
  - **Impact** — does it fill a real, confirmed hole in the open-source ecosystem?
  - **Tractability** — buildable as a small, tested, CPU-only tool in one slice?
  - **Data availability** — a lightweight, openly-downloadable, programmatically-accessible dataset?
  - **Offline reproducibility** — works with a synthetic fallback; no network/GPU/HPC needed?
  - **jansky reuse** — composes an existing tested helper rather than reimplementing it.
- A single committed decision: **ONE domain + ONE dataset + ONE helper**, plus the `jansky` extra
  to request (`pulsar` / `hi` / `ml` / `seti` / `formats`).

## Approach

Synthesize the survey; score each gap; **science-reviewer** sanity-checks that the leading gap is
genuinely unsolved (not a maintained package the survey missed). Present the ranked table to the
human.

## Verification

- **GATE 1 (hard stop):** the human approves the domain + dataset before Plan 03 begins.
- The chosen dataset is confirmed openly downloadable and lightweight (or has a synthetic
  fallback); the chosen gap is confirmed not already covered by a maintained tool.
