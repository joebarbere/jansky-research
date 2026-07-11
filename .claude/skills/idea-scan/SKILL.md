---
name: idea-scan
description: Run a deep, multi-agent opportunity re-scan of the radio-astronomy research landscape and produce a plan-ready ideas file. Use when the user asks for new research opportunities, a fresh look at the field, "what should we work on next", or to refresh/replace an earlier ideas file (fable-ideas.md, survey/opportunity-scan-*.md).
---

# Idea scan — deep opportunity re-survey

Produce a **plan-ready ideas file** (each entry turnable into a `plans/NN-*.md` by a planning
model without re-doing the survey). This codifies the process behind
`survey/new-findings-scan.md`, `survey/opportunity-scan-2026-07.md`, and `fable-ideas.md` —
read whichever is newest first; the new scan must go *beyond* it, not repeat it.

## Procedure

1. **Ground in the repo's current state** (30 min, no network): `README.md` slice table,
   the newest `survey/*-findings.md` "caveats/future work" sections, open `plans/`,
   and the workstation profile (GPU/ROCm status, disk, station timeline) from the latest scan.
   List what is *merged* — the biggest failure mode is proposing already-executed ideas.
2. **Fan out parallel research agents** (Agent tool, background, one message). Proven split:
   one `radio-research-assistant` per domain (FRB/time-domain, LPTs/coherent transients,
   catalogue cosmology/statistics, solar/heliospheric, planetary/exoplanet, pulsars/PTA,
   HI/spectral-line, rooftop-station science), one `archive-scout` to live-verify data access,
   one deliberately speculative wild-card (general-purpose) told to **kill ideas that fail the
   "one person, public data, checkable, null-worth-something" bar and say why**, and optional
   GPU/ML sweeps (foundation models, SBI, DSP ports). Each prompt must be self-contained:
   repo track record, hardware constraints ("pure PyTorch runs; CUDA kernels don't"), and the
   required per-idea fields (below).
3. **Demand evidence-of-novelty discipline** from every agent: name the *closest* published
   work (arXiv IDs), say why it doesn't cover the idea, and flag every claim that rests on
   search snippets rather than full-text reads. In egress-restricted sessions WebSearch works
   but WebFetch to astronomy hosts may 403 — record that as a standing GATE-0 ("re-verify by
   full-text pass from the workstation"), don't hide it.
4. **Synthesize** into a single file with, per idea: gap · novelty evidence (IDs) · exact data
   products + access path + **GATE-0** · method sketch naming reusable `src/jansky_research`
   modules · recover-a-known validation · biggest risk · hardware fit. Include a
   **"Corrections & closed doors"** section (ideas killed and why — as valuable as the
   survivors) and an opinionated "first moves" ranking. Tier by novelty-ceiling × readiness.
5. **Commit on a branch** (never `main`), push, and cross-link: the new file supersedes the
   previous scan's shortlist — say so at the top of the new file.

## House rules that apply here

- Honest framing: a slice whose likely outcome is a null must be pitched as a limit/method
  paper from day one.
- Every idea needs a recover-a-known before any new claim — if none exists, that's a red flag.
- Time-sensitive items (collaboration-scoop risk) are marked and ranked accordingly.
