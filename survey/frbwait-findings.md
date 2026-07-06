# Findings — uniform Cat-2 repeater wait-time/duty-cycle census (plan 39)

`jansky_research.frbwait` computes the first ONE-statistic timing census across every CHIME/FRB
Catalog 2 repeater: Weibull wait-time shape k + Rayleigh-Z² periodogram + activity-window duty
cycle, with significance from a transit-comb-preserving scramble.

## GATE 0 (2026-07-06)

- **chime-frb.ca is still 503** (whole domain) — but the official release lives on CANFAR:
  DOI 10.11570/25.0066, VOSpace `AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/` on
  `cadc-west-01.canfar.net` (the `ws-cadc` vault path 404s; resolve via `/reg/resource-caps`).
  Mirrored: `data/chimefrbcat2.csv` (4.1 MB, 5,045 rows incl. sub-bursts → 4,539 events, 83
  repeaters, 1,282 repeater rows) + `data/chimefrbcat2_exposure.h5` (216 MB). The catalog paper
  is ApJS doi:10.3847/1538-4365/ae3828 (arXiv:2601.09399).
- **Novelty confirmed by full-text pass**: the Cat 2 paper defers all repeater population work
  to Cook et al. (arXiv:2605.08410), which computes rates/DM-drift/repeater-fraction but **no
  Weibull k, no periodograms, no duty cycles** (quote: "Repeater burst arrival times are
  stochastic; however, they cluster in time" — unmodelled). No other uniform census exists.
- **Plan correction**: FRB 20240114A (and 20240209A) post-date the Cat-2 cutoff (2023-09-15)
  and are NOT in the table — the plan's second anchor is impossible from this dataset. Also its
  "~112.9 d" is a *chromatic* (central-frequency) period (arXiv:2605.12098), not an activity
  period. Single in-catalogue anchor: FRB 20180916B, 16.35 d.
- **Exposure product limitation**: `chimefrbcat2_exposure.h5` = two nside-4096 HEALPix maps of
  TIME-INTEGRATED exposure (2018-09-04→2023-09-15, upper/lower transit). No per-epoch exposure
  history is public → the plan's "exposure-corrected periodogram" is scoped to (a) per-source
  total-exposure rate normalisation (the catalogue's own `exp_up` hours) + (b) a
  **transit-comb-preserving scramble null**: keep each burst's sidereal phase (0.99727-d comb),
  redraw its day number uniformly over the span. Exposure-blind analytic FAPs (the Cat-1
  lesson) are never used.

## Recover-a-known

- Synthetic transit-sampled Weibull train (k=0.4, 16.35-d period, 30% duty): all three
  statistics recovered (k within CI, period to grid resolution, duty within transit-quantised
  tolerance); Poisson control stays null. Runs offline in CI.
- **Real anchor: FRB 20180916B → 16.33 d at the scramble floor (p=0.001, 999 scrambles), 107.1
  cycles, duty(90%)=0.214, k=0.42 [0.36, 0.50]** — against the published 16.35 d and ~5-d
  window (CHIME/FRB 2020, Nature 582, 351).

## Result (999 scrambles; full census in `results/frbwait_metrics.json`)

- 83 repeaters censused (rate rows); **15 above the completeness cut** (≥10 bursts, ≥30-d
  span) get k/periodogram/duty.
- **Median k = 0.83**; 3 sources clustered at 95% (k CI < 1): 20220912A (0.30), 20201124A
  (0.34), 20180916B (0.42) — the well-sampled sources are all strongly sub-Poissonian; most
  N~10–35 sources have k CIs straddling 1 (honestly wide).
- **3 sources beat the scramble at p≤0.01** — but only the anchor has many cycles:
  20180916B 16.33 d @ 107 cycles (real), 20201124A 60.1 d @ 5.0 cycles and 20220912A 122.7 d @
  3.0 cycles are **activity-epoch degeneracies, labelled as such, NOT new period claims** (the
  paper's stated N_cyc≲5 warning label).

## Interpretation limits (GATE-2 material)

- Duty cycles inherit the transit-sampling approximation (no public per-epoch exposure); the
  scramble null carries the comb but not instrument downtime.
- **Transit-visibility censoring of waits (GATE-2 required disclosure):** waits between ~15 min
  and ~1 sidereal day are unobservable (source set), so observed waits are bimodal
  (intra-window | ≥1 d) and the Weibull k is biased DOWNWARD (toward apparent clustering),
  rate-dependently — worst for the highest-rate source (20220912A, k=0.30, many 306-s
  intra-window waits). The census k values are conditioned-on-CHIME-sampling comparative
  statistics, not deconvolved intrinsic shapes; stated in the paper.
- Fluence incompleteness additionally clips faint bursts (waits merged, k biased upward) —
  smaller and roughly uniform across sources.
- min 10-burst cut: k for low-N sources is honestly unconstrained; population claims are made
  only above the cut.

## Reproduce

`uv run python -m jansky_research.frbwait --n-scramble 999 --out .` (~4 min CPU with the local
mirror). Offline CI leg: `--offline`.
