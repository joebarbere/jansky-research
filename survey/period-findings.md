# Findings — FRB repeater activity-periodicity search (CHIME/FRB Catalog 1)

Run of `jansky_research.frbperiod.run` over the 18 repeating sources in the First CHIME/FRB Catalog
(one row per event, `sub_num == 0`), using a phase-folding **Rayleigh $Z^2_1$ periodogram** on each
source's burst arrival times (MJDs). Full per-source table: `survey/period_results.csv`.

## Result: the known 16.35-day period is recovered; nothing new

- **FRB 20180916B — period $= 16.33$ days** ($Z^2_1 = 33.4$, 19 bursts, exposure-blind
  $\mathrm{FAP} = 7\times10^{-6}$). This **independently recovers, from the public catalogue alone,**
  the famous 16.35-day activity period first reported by CHIME/FRB Collaboration (2020, Nature 582,
  351; $16.35\pm0.18$ d). The agreement is essentially exact — this is the **validation** that the
  method works. See `paper/figures/periodogram.pdf`.
- **FRB 20180814A** (11 bursts, the only other searchable source) shows **no significant period**:
  the best peak is $Z^2_1 = 9.4$ at 2.8 d with $\mathrm{FAP} \approx 0.76$ — fully consistent with
  noise. (A $\sim$157-day period has been claimed for this source elsewhere, far beyond the
  $\sim$350-day Catalog-1 span and its sampling, so it is not recoverable here.)
- **The other 16 repeaters have $\le 3$ bursts each** in Catalog 1 and are not searchable
  (`min_bursts = 8`).

So: **one detection, and it is the already-known one.** No new periodicity — the public catalogue
is simply too sparse (only 2 of 18 repeaters have enough bursts).

## Honest limitations

- **Exposure-blind significance.** CHIME is a transit survey: each source is seen $\sim$once per
  sidereal day with strongly non-uniform exposure. The analytic FAP ignores this, so the
  $7\times10^{-6}$ is an *upper bound on confidence*, not a rigorous detection significance — the
  real confirmation is the **match to the independently-published period**, not the FAP.
- **Daily aliasing.** A catalogue-only periodogram aliases near 1 day and its beats; 16.33 d sits
  well above that and coincides with the known value, so the recovery is robust, but blind peaks at
  short periods (e.g. the 2.8 d for FRB 20180814A) are not trustworthy.
- **Sparse sampling.** With $\le 3$ bursts, 16 of 18 sources yield no constraint at all; even the
  two searchable sources are marginal by professional standards (CHIME's own detection used a
  dedicated monitoring campaign plus the exposure model).

## Bottom line

A clean, honest **validation** result: the open, tested, CPU-only tool recovers FRB 20180916B's
16.35-day activity period from public Catalog-1 data, and correctly finds nothing significant for
the rest — no over-claimed new periodicity. The contribution is the reproducible periodicity-search
tool plus this validation and the per-repeater limits; finding *new* repeater periods would need the
denser, exposure-modelled data the professional teams use, not Catalog 1.
