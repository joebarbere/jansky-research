# 48 — JBO glitch catalogue: post-2018 population statistics (≥120 unanalysed glitches)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — and an explicit check
for an in-press JBO population paper before committing

## Context

The Jodrell Bank glitch catalogue now stands at 664 glitches across 207 pulsars, but the last
population-statistics paper (Basu+2022, arXiv:2111.06835) stops at end-2018 (543 glitches / 178
pulsars) — ≥120 glitches have never entered a population analysis. Recent bimodality work is
single-pulsar only (arXiv:2502.20017). Data: scrape
`jb.man.ac.uk/pulsar/glitches/gTable.html` (access pattern already verified in earlier repo
scans). Method: for every pulsar with ≥5 glitches, BIC-selected exponential-vs-mixture
waiting-time fits plus size-distribution classification; the headline is which pulsars'
classifications *change* when the post-2018 glitches are added. This reuses the `ppdot`/`lpt`
compilation pattern (provenance-typed scrape, flag discipline, ATNF cross-match). CPU, days.

## Deliverables

- `src/jansky_research/glitchpop.py`: `scrape_glitch_table` (HTML → typed table with provenance
  flags, `# pragma`), `waiting_time_fit` (exponential vs mixture, BIC selection),
  `size_distribution_fit` (lognormal / power-law / bimodal classification),
  `classification_delta` (Basu+2022-epoch subset vs full catalogue, per pulsar),
  `synthetic_glitch_series` (injected exponential and mixture waiting-time processes →
  recovery), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/glitchpop/`; `survey/glitchpop-findings.md`; wiring.

## Approach

0. **GATE 0:** confirm `gTable.html` is live and its current row count (~664); full-text pass
   on Basu+2022 (arXiv:2111.06835) and arXiv:2502.20017; **search ADS/arXiv explicitly for an
   in-press JBO catalogue-update paper** — the team publishes every ~3–4 yr and one may be due;
   abort or rescope if found.
1. Tooling + synthetic recover-a-known: injected exponential and two-component-mixture waiting
   times at realistic per-pulsar glitch counts; BIC selection must classify both correctly.
2. Real leg: scrape, type, and flag the full table; per-pulsar (≥5 glitches) fits; the
   classification-delta table against the end-2018 subset (reproducing Basu+2022's own
   classifications on that subset first).
3. GATE-2 science review: detection-completeness caveats (small glitches missed unevenly across
   monitoring history), the ≥5-glitch cut's selection, honest framing of classification changes
   as data-driven not physical claims.
4. Paper: updated population statistics + the change table.

## Verification

Crab/Vela must come out with their known exponential waiting-time behaviour and J0537−6910 with
its known distinct behaviour before any new classification is reported; the end-2018 subset must
reproduce Basu+2022; checks green; GATE-2 sign-off.

## Risks & mitigations

- **JBO team publishes an update first** (their ~3–4 yr cadence) → the GATE-0 in-press check is
  mandatory; if something is in press, drop or narrow to an angle their paper does not cover.
- **Scrape fragility / transcription errors** → provenance flags and spot-checks against
  discovery papers, the `lpt` flag discipline.
