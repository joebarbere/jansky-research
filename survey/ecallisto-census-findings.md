# Findings — type III occurrence census vs the solar cycle (method + recover-a-known)

`jansky_research.ecallisto_census` turns the coincidence-vetted type III event stream (from the
e-Callisto ingest, `survey/ecallisto-findings.md`) into an **occurrence census**: a rate that can be
compared honestly across epochs and correlated against a solar-activity index. This slice builds and
validates the census *statistic*; it does not report a measured solar-cycle correlation (the real
multi-cycle ingest is future work — see the honest limits below).

## Why a coverage correction is the whole point

The expectation is well established: the type III rate tracks solar activity over the ~11-year cycle
(Saint-Hilaire et al. 2013; Reid & Ratcliffe 2014). But a raw confirmed-event count from e-Callisto
conflates real activity with **how many stations were watching** — a burst is confirmed only when
enough stations observe it, and the active-station count grew over the network's history and varies
day to day. The census removes that confound with the minimal correction

```
rate = N_events / C          (C = active-station coverage; zero coverage → undefined, not zero)
```

then correlates the corrected rate with the SILSO monthly sunspot number (Pearson r, Spearman ρ, OLS
slope). If a burst's confirmation probability scales with coverage, `N_events ∝ R·C`, so dividing by
`C` returns an estimator of the underlying rate `R` that is comparable across epochs with different
network sizes.

## Recover-a-known (synthetic event stream, real activity driver)

A full multi-cycle real ingest is a large data-collection task (below), so the statistic is validated
on a synthetic observing history for which the truth is known: a realistic fast-rise/slow-decay
sunspot cycle, and for each month with sunspot `S` and randomly varying coverage `C`, an event count
drawn from `Poisson(k·S·C)` with `k = 0.03`. By construction `N/C` has expectation `k·S`.

| quantity | value |
|---|---|
| synthetic months | 180 |
| total synthetic events | 4651 |
| Pearson r (rate vs sunspot) | **0.968** |
| Spearman ρ | 0.963 |
| OLS slope (events·station⁻¹ per unit sunspot) | **0.0302** (recovers injected k=0.03) |

The correlation emerges **through** the coverage variation: the raw count is contaminated by the
fluctuating station count, and only after the `N/C` correction does the clean activity signal appear.
The census statistic and its implementation are validated; what remains is data.

## Honest assessment & caveats

- **The event stream here is synthetic.** The sunspot driver is a realistic cycle (and the *real*
  SILSO series is parsed by the same code for the real run), but the monthly event counts are drawn
  from a model, not measured. This is a method validation, **not** a measured solar-cycle detection.
- **The real census is coverage- and detection-limited.** As the coincidence QC found, a small
  station subset with synthetic-tuned thresholds yields a coverage-limited lower bound (a six-station
  snapshot of a known burst day gave zero confirmed events). A real census needs the full active-station
  set per window and event-tuned detection — a many-years × full-coverage ingest, which is exactly what
  the scheduled, backfilling Airflow pipeline was built to run. This is a data task, not a method gap.
- **Type III only.** The underlying detector targets fast negative-drift ridges; a type II census
  would need a second, slower-drift template.
- **The coverage correction is first-order.** Dividing by the active-station count assumes confirmation
  probability scales with coverage; a fuller treatment would model per-station sensitivity and sky
  coverage. `rate = N/C` is the minimal, transparent correction, not the last word.
- **Reproducible:** `python -m jansky_research.ecallisto_census --offline` reproduces the validation
  (figure + macros + `results/ecallisto_census_metrics.json`); the same statistic runs on a real event
  stream once the multi-year ingest is collected.

GATE-2 (science-reviewer): no blockers; all six primary citations verified; the synthetic-vs-real
boundary is explicit in the abstract, results, figure caption, and discussion.

## Referee round on the style conversion (2026-08-23): scope words are load-bearing

Verdict *minor revision*, two MAJORs, both created by the restyle and both fixed:

1. `all event counts in **this paper** are synthetic` had become `in this validation` --- a
   document-scoped guarantee narrowed to a section-scoped one, leaving the abstract's
   `\ecsNevents` outside the stated scope.
2. `This paper validates the method; it does not report a measured solar-cycle correlation` had
   become agentless passive (`a measured solar-cycle correlation is not reported`), which read as
   a claim about the literature --- false on the paper's own account, since it cites
   \citep{sainthilaire2013} for exactly that. Restored an explicit subject.

Also fixed: a pronoun left pointing at "a method gap" after a sentence split; "the recovery works
**through** the coverage variation" (the italic carried *in spite of*, and unemphasised it reads as
*by means of*, reversing the claim) reworded to "survives"; the cross-station veto's premise
restored as a premise rather than an assertion.

**Outstanding, pre-existing:** `results/ecallisto_census_real_metrics.json` holds a committed real
leg (168 sampled days, 5 events, 2 days with events, `pearson_r` 0.28) that the paper never
mentions. It is consistent with the paper's stated position --- the real census is a data-volume
task --- and 5 events cannot support a correlation, but committed real evidence that the paper is
silent on should either be cited or explained.
