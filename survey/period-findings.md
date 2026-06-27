# Findings — FRB repeater activity-periodicity search (CHIME/FRB Catalog 1)

Run of `jansky_research.frbperiod.run` over the 18 repeating sources in the First CHIME/FRB Catalog
(CHIME/FRB Collaboration 2021, ApJS 257, 59; one row per event, `sub_num == 0`), using a
phase-folding **Rayleigh $Z^2_1$ periodogram** (Buccheri et al. 1983) on each source's burst arrival
times (MJDs). Full per-source table: `survey/period_results.csv`.

## Result: the known 16.35-day period is recovered; nothing new

- **FRB 20180916B — period $= 16.33$ days.** The periodogram peak coincides with the famous
  16.35-day activity period first reported by CHIME/FRB Collaboration (2020, Nature 582, 351;
  $16.35\pm0.15$ d) — recovered here **from the public catalogue alone**. The recovered value
  (16.326 d, 19 bursts) sits within the grid resolution (0.008 d) and within the published
  $1\sigma$ uncertainty ($0.16\sigma$ away). **That coincidence with the independently-published
  value is the validation** — not the analytic peak strength ($Z^2_1=33.4$; exposure-blind
  $\mathrm{FAP}=7\times10^{-6}$, which ignores survey sampling and is *not* a rigorous significance,
  see below). See `paper/figures/periodogram.pdf`.
- **FRB 20180814A** (11 bursts, the only other searchable source) shows **no significant period**:
  the best peak is $Z^2_1=9.4$ at 2.8 d with $\mathrm{FAP}\approx0.76$ — consistent with noise. No
  multi-day activity period has been reported for FRB 20180814A in the literature, so this null is
  expected. (The $\sim$157-day period in the literature belongs to a *different* source,
  FRB 20121102A — Rajwade et al. 2020 — not to FRB 20180814A.)
- **The other 16 repeaters have $\le 3$ bursts each** in Catalog 1 and are not searchable
  (`min_bursts = 8`). FRB 20121102A appears with a single Catalog-1 event only because most of its
  many hundreds of known bursts pre-date the catalogue's 2018–2019 survey window — it is not a
  sparse repeater in the wider literature.

So: **one detection, and it is the already-known one.** No new periodicity — the public catalogue is
simply too sparse (only 2 of 18 repeaters have enough bursts).

## Honest limitations

- **Exposure-blind significance.** CHIME is a transit survey: each source is seen $\sim$once per
  sidereal day with strongly non-uniform exposure. The analytic FAP treats trial frequencies as
  independent, which the aliased transit window violates, so it is an *approximate* number, not a
  rigorous detection significance. The real confirmation is the match to the published period.
- **Topocentric times.** The Catalog-1 MJDs are topocentric (after DM correction to infinite
  frequency); no barycentric correction was applied. For a 16-day period the resulting Roemer-delay
  modulation ($\sim$0.03% of the period) is negligible and the correct recovery confirms it — but
  this would matter for any search at sub-day periods.
- **Daily aliasing.** A catalogue-only periodogram aliases near 1 day and its beats; 16.33 d sits
  well above that and coincides with the known value, so the recovery is robust, but blind peaks at
  short periods (e.g. the 2.8 d for FRB 20180814A) are not trustworthy.
- **Sparse sampling.** With $\le 3$ bursts, 16 of 18 sources yield no constraint; even the two
  searchable sources are marginal by professional standards (CHIME's own detection used a dedicated
  monitoring campaign plus the exposure model).

## Bottom line

A clean, honest **validation** result: the open, tested, CPU-only tool recovers FRB 20180916B's
16.35-day activity period from public Catalog-1 data, and correctly finds nothing significant for
the rest — no over-claimed new periodicity. The contribution is the reproducible periodicity-search
tool plus this validation and the per-repeater limits; finding *new* repeater periods would need the
denser, exposure-modelled data the professional teams use, not Catalog 1.
