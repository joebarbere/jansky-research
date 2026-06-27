# Plan 01 — Deep-research survey 🚀

> Context: the evidence base for the gap pick. Pairs with 02 (gap analysis). Scope: medium.

## Context

Before building anything we survey the landscape so the gap we fill is *real* (not already
solved) and *tractable* (CPU-only, offline-reproducible, lightweight public data). Three strands:
the GitHub open-source tooling landscape, the current scientific literature, and the
data/automation ecosystem.

## Deliverables

- `survey/literature.md` — state-of-the-art + active questions per candidate domain
  (FRBs/transients, pulsar timing, HI/continuum, RFI/ML, SETI), each with verified citations.
- `survey/github-landscape.md` — the open-source tooling map (what exists, maturity, and what is
  MISSING/underserved), plus the state of pipeline/automation/orchestration tooling.
- A per-domain shortlist of concrete gap hypotheses with a candidate lightweight public dataset
  and the `jansky` helper each would build on.

## Approach

Fan out agents (single batch, parallel):

- **radio-research-assistant** — one invocation per domain: state-of-the-art, the GitHub tooling
  gap, lightweight datasets, and 2–4 gap hypotheses building on jansky helpers.
- **general-purpose** — the cross-cutting GitHub/automation ecosystem sweep (orgs, packages,
  coverage map, the Airflow/Prefect/Snakemake-in-radio-astronomy question, awesome-lists).
- Skills **find-radio-papers** / **dataset-watch** as needed.

Every package/dataset/paper is URL-verified — no invented citations.

## Verification

- The two survey docs exist, are concrete about exists-vs-missing, and cite live URLs.
- At least 4–6 distinct, plausible gap hypotheses are captured for 02 to rank.

**Status: in progress** — 6 agents fanned out; results being synthesised into `survey/`.
