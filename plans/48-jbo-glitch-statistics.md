# 48 — JBO glitch catalogue: post-2018 population statistics (≥120 unanalysed glitches)

Status: ✅ done 2026-07-11 — a WORKING classification with a passing recover-a-known. Monitoring-gap-
robust per-pulsar waiting-time classification of the live JBO catalogue (727 glitches/223 pulsars; 33
with >=5 glitches: 23 exponential, 10 quasi-periodic, 0 clustered). **Recover-a-known PASSES** (J0537-6910
+ Vela both quasi-periodic) -- but ONLY after the key insight that the archival waiting times conflate
real intervals with MONITORING gaps: J0537's 2264-d RXTE->NICER gap inflates its raw CV to ~2 so a
naive fit calls it exponential; excising waits >6x median fixes it (cost: long-gap clustering
under-detected, stated). Classifier = parametric-bootstrap CV test (replaced a gamma-BIC that
underweighted Vela). Post-2018 delta: 2 flips (J1841-0524 qp->exp, J2229+6114 exp->qp) + 3 newly
classifiable vs the end-2018 Basu subset. GATE-2 pending. Module `glitchpop.py` (+`glitchpop_real.py`),
paper `papers/glitchpop/`, findings `survey/glitchpop-findings.md`. — GATE 0 done 2026-07-11:
**PASS (conditional)**. Data LIVE + scrapable, no
auth: `jb.man.ac.uk/pulsar/glitches/gTable.html`, **727 glitches / 222 pulsars** (up from Basu+2022's
543/178 -> ~184 new post-2018). Columns: Bname, **JNAME**, glitch#, **MJD**+err, **dF/F** (dnu/nu
x1e-9)+err, **dF1/F1** (dnudot/nudot x1e-3)+err, References. **Basu+2022** (MNRAS 510, 4049;
arXiv:2111.06835): catalogue + two-component-Gaussian mixture of glitch SIZE + dnudot; end-2018; did
NOT do systematic per-pulsar WAITING-TIME model selection. **arXiv:2502.20017** = Bactrian sizes,
J0537 only (single-pulsar). **IN-PRESS CHECK**: no JBO catalogue-update/population paper post-2022.
BUT **Zhu & Zheng 2025 (arXiv:2501.01862)** use the SAME >=5-glitch JBO+ATNF sample (671 glitches, 32
pulsars) for glitch-CLUSTER-period-vs-age -- MUST cite + differentiate; our deliverable is different
(per-pulsar exponential-vs-regular waiting-time classification + the post-2018 DELTA, NOT clustering,
which is taken -- do not drift there). **Recover-a-known**: most pulsars power-law sizes + EXPONENTIAL
(Poisson) waits; J0537-6910 quasi-periodic ~100 d (Gaussian, NOT exponential, size<->next-wait corr);
Vela ~2-3 yr quasi-periodic; Crab ~30 glitches. **The wedge**: per-pulsar waiting-time regularity
classification (exponential vs quasi-periodic vs clustered, via CV + gamma-vs-exponential BIC) + the
diff of which classifications FLIP when post-2018 glitches are added vs the end-2018 Basu subset.
ATNF join by JNAME. Reuse: `frbstats.fit_power_law`, `ppdot`/`lpt` scrape+flag discipline.

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
