# 71 — PINT quick-look triage toolkit on NANOGrav 15yr + PPTA DR3

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

Pulsar-timing datasets ship as par/tim files, and every group re-invents the same first-look
triage: residual plots, a normalized power spectrum, a step-function/small-glitch matched
filter, and EFAC/EQUAD sanity checks. A small, tested, PINT-based quick-look toolkit run
uniformly over the two big *public* datasets — NANOGrav 15yr (Zenodo 7967584/8060824/8092346)
and PPTA DR3 (CSIRO DAP 10.25919/w0nw-jt05, 10.25919/23wj-1d69) — is methods/JOSS-adjacent
value (fable-ideas F34). IPTA DR3 is NOT public (fable-ideas corrections), so these two are the
right substrates. PINT's NGC 6440E tutorial dataset is vendored as the offline test fixture.
Any "candidate glitch" the matched filter flags must be gated against chromatic-noise nulls —
achromatic step vs DM-noise look-alikes is the known trap.

## Deliverables

- `src/jansky_research/pintlook.py`: `load_partim` (PINT wrapper, NGC 6440E fixture offline),
  `residual_psd` (normalized PSD with uneven-sampling handling), `step_matched_filter`
  (glitch-candidate scan + significance from scrambles), `noise_triage` (EFAC/EQUAD quick fit),
  `fetch_nanograv15` / `fetch_ppta_dr3` (network, `# pragma`), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/pintlook/`; `survey/pintlook-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text pass on the NANOGrav 15yr and PPTA DR3 data papers plus their 2026
   noise/profile-event papers (arXiv:2604.05453 closed the profile-event scan — scope around
   it); verify the Zenodo/DAP DOIs resolve; vendor NGC 6440E and pin the PINT version.
1. Tooling + synthetic recover-a-known: inject a step of known epoch/amplitude into NGC 6440E
   residuals; the matched filter recovers it at the stated significance; injected white-noise
   levels recovered by the EFAC/EQUAD triage.
2. Real leg: run the full quick-look uniformly over all 15yr + PPTA DR3 pulsars; publish the
   triage tables and per-pulsar summary figures.
3. Gate any candidate step against chromatic-noise nulls (refit with DM-noise terms; a candidate
   that dissolves is reported as dissolved).
4. GATE-2 (published-noise-model comparison; matched-filter trials factor) → paper
   (methods/JOSS-adjacent note; the toolkit is the deliverable, not any single candidate).

## Verification

Injected steps recovered in the NGC 6440E fixture; per-pulsar noise triage consistent with the
collaborations' published noise models where comparable; chromatic-noise nulls gate every
candidate; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Chromatic noise mimics steps** → the null gate is mandatory; no candidate is reported
  without surviving a DM-noise refit.
- **Collaborations' own pipelines are deeper** → framed as a quick-look/triage tool with tests
  and a uniform public run, not a competitor to enterprise noise analyses.
- **PINT API churn** → pin the version; the fixture-based test suite catches breakage.
