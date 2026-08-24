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

## Full referee round (2026-08-24): MAJOR REVISION, 15 findings, one BLOCKER

The anchor recovery is real (16.325 d vs published 16.35 +/- 0.15 at 107 cycles) and the
refusal to call 3-5-cycle peaks periods is right. The population claim is the problem.

**BLOCKER: the three "clustered" sources are exactly the three highest-rate sources** ---
probability 1/C(15,3) = 0.0022 under no rate dependence --- which is precisely the ordering the
paper's own disclosed, unquantified, rate-dependent censoring bias predicts. The three periodic
detections are the same three sources. "We disclose rather than correct it" does not survive the
disclosed effect reproducing the entire result set. The decisive experiment: recover k from
transit-sampled synthetics with known k=1 across the observed rate range (0.02-3.23 bursts/hr).

**MAJORs:** two sources are significantly SUPER-Poissonian (k_ci_low > 1) and never mentioned ---
the 95% exclusions split 3:2 in opposite directions, which "the population norm" cannot survive;
the median k=0.831 has no uncertainty and a sign test gives two-sided p=0.12; the offline
validation injects at ~25-30x below the worst-case source's rate, so the censoring bias is absent
from the validation and it cannot fail; the scramble redraws each burst's day independently,
destroying within-transit multiplicity, so the FAP tests "more than one burst per transit"
rather than periodicity --- largest for exactly the three detected sources; FRB 20220912A's
"122.67 d" peak is railed at the grid edge (span/3 exactly), so its N_cyc=3 is forced, not
measured.

Nine MINOR/NIT: duty-cycle consistency with the published 5-d window asserted but the
convention-matched containments never computed; n_scramble/seed/n_boot not in the JSON; the
dec-completeness claim traces to no committed number (rows carry no dec); n_boot=200 decides
boundary cases (two sources sit on the CI edge); the N_cyc rule is stated as <~5 but coded >=10
with 20201124A at 5.008 rendered as "5"; the Discussion calls non-deconvolved k values a
"population-synthesis input" two sections after saying they are not; cat2 bib entry lacks the
published volume/page (ApJS 283, 34); cat1 invoked but never cited; the macros header claims the
opposite of the guard it calls, and a test asserts the placeholder behaviour.

**Status: fixes pending** (data local: data/chimefrbcat2.csv + exposure h5).
