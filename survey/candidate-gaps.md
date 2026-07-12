# Candidate gaps — preserved backlog (DO NOT DELETE)

The deep-research survey (P1) surfaced ~24 concrete open-source gap hypotheses across five
domains. GATE 1 committed the vertical slice to **one** of them (FRB burst-statistics →
`jansky_research.frbstats`). **This file preserves all the *other* candidates** — each is a
plausible future vertical slice (its own tool + dataset + paper). Keep this file even after the
`plans/` files are deleted post-merge; it is the project backlog, not a plan.

Every gap below is CPU-only, offline-reproducible, lightweight-data, and composes an existing
tested `jansky` helper. Datasets and packages were URL-verified in the survey (see
`literature.md`, `github-landscape.md`).

## ✅ Completed slices (built, reviewed, merged)

1. **FRB burst-statistics** → `jansky_research.frbstats` — reproduced the CHIME width result.
2. **Multi-survey spectral-index / USS hunt** → `jansky_research.spectra` — candidates did not
   survive the de Gasperin cross-check (honest negative).
3. **FRB repeater periodicity** → `jansky_research.frbperiod` — recovered FRB 20180916B's 16.35-day
   period.
4. **SETI injection-recovery benchmark** → `jansky_research.driftsearch` — benchmark built; the
   Voyager real-data check is an honest negative (DC-spike, not the carrier).
5. **HI tangent-point rotation curve** → `jansky_research.hi` — recovered the flat Milky Way curve
   (dark-matter signature).

The remaining gaps below are unbuilt.

---

## FRBs & transients (chosen domain — non-chosen gaps)

- **[CHOSEN] burst-statistics** → `jansky_research.frbstats` (Weibull wait-times, power-law
  energy, repeater vs non-repeater). Data: CHIME/FRB Cat 1 CSV.
- **DM → z pipeline with Macquart scatter** — offline `(DM, l, b) → DM_MW (pygedm) → DM_excess →
  z ± σ` from the DM_host log-normal prior. Fills the stale `fruitbat`. Reuse:
  `jansky.transients.macquart_redshift`, `jansky.constants.MACQUART_SLOPE`. Data: TNS / FRBSTATS
  CSV (tens of KB). Refs: Macquart 2020 (Nature 581, 391); arXiv:2511.01195 (F≈0.32).
- **Burst morphology from filterbank** — batch drift-rate ("sad trombone"), component count,
  temporal/spectral width from a `.fil` + DM. Fills stale `burstfit`/`dfdt`. Reuse:
  `jansky.transients.dedisperse/boxcar_snr`, `jansky.formats.read_filterbank`. Data: FRB 140514
  Parkes `.fil` (~100–200 MB/beam); BL FRB 121102 example (jansky vendors a ~1.6 MB `.fil`).
- **Population-consistent injection → recovery efficiency** — inject FRBs drawn from a fluence
  power-law + Macquart DM(z), run a DM search, report recovery vs fluence. Reuse:
  `jansky.transients.disperse_pulse/dm_search`. Telescope-agnostic; no open tool does this.

## Pulsar timing & PTAs

- **Residual quick-look** — `par,tim → post-fit residuals + plot + RMS/χ²` in <10 lines
  (thin PINT wrapper; none exists). Reuse: `jansky.timing`, `jansky[pulsar]` (PINT). Data:
  NGC 6440E `.par`/`.tim` **already vendored** (~4 KB); any of 68 NANOGrav 15yr pulsars (KB each,
  Zenodo 8423265).
- **Residual red-noise spectrum** — correctly-normalised Lomb–Scargle PSD of residuals + fitted
  power-law slope (the documented astropy `psd`-mode gotchas, packaged). Reuse: `jansky.timing`,
  astropy `LombScargle`. Refs: bgoncharov.com/pulsar-timing-software.
- **Small-glitch scan** — step-function / matched-filter search for phase jumps in residuals
  (small glitches are missed by eye — IAR 2024, A&A). Reuse:
  `jansky.transients.epoch_folding_search/fold_profile`.
- **PTA noise budget** — laptop EFAC/EQUAD + red-noise-index triage ("is this MSP PTA-quality?")
  without enterprise/HPC. Reuse: `jansky.timing.simulate_pta_residuals`, `scipy.optimize`.

## HI 21 cm, Galactic structure & continuum source counts

- **Tangent-point rotation-curve extractor** — `(l,v) FITS slice → (R, v_rot, σ)` with
  uncertainty; the Ch 11 `# TODO` made into a tested function. Reuse:
  `jansky.data.synthetic_hi_cube`, spectral-cube. Data: **LAB (l,v) slice vendored** (366 KB);
  more longitudes at VizieR VIII/76 (KB each).
- **log N–log S fitter with uncertainty** — Poisson-weighted / Crawford-Murdoch ML slope with
  asymmetric CIs (current `jansky.sourcecounts.count_slope` is OLS only) + de Zotti 2010 overlay.
  Data: NVSS cone search via astroquery/VizieR (VIII/65), no auth.
- **Multi-survey spectral-index cross-matcher** — NVSS × FIRST cone search → α ± σ_α (peaked /
  flat / steep sources). Reuse: `jansky.sourcecounts`, `astroquery.vizier` (VIII/65, VIII/92).
- **Unified offline HI (l,v) pipeline** — one dict from synthetic *or* the vendored LAB slice,
  feeding the tangent-point extractor; fully offline. Doubles as the extractor's test harness
  (synthetic cube has a known ground-truth velocity field).

## RFI mitigation & ML

- **CPU-only FRB/RFI sklearn classifier** — feature-vector (`jansky.rfi` SK/MAD/flag-fraction +
  `jansky.transients` DM-SNR curve) → RandomForest joblib; fills the GPU/TF-heavy FETCH gap.
  Reuse: `jansky.rfi`, `jansky.transients`, `jansky[ml]`. Data: HTRU2 CSV (<1 MB, Figshare
  3080389); HTRU1 images (183 MB, Zenodo 3205409).
- **SK-vs-SumThreshold reproducible comparison harness** — both flaggers on the same data, same
  threshold convention → ROC / flag-map; addresses the field's reproducibility hole. Reuse:
  `jansky.rfi.spectral_kurtosis/sumthreshold2d`. Data: synthetic (`disperse_pulse`) or BL open
  filterbank.
- **Tiny canonical pixel-level RFI benchmark builder** — generate `.npz` dynamic spectra + masks
  with a fixed split (a sub-100 MB HTRU2-equivalent for flagging; none exists). Reuse:
  `jansky.transients.disperse_pulse`, `jansky.rfi`. No download needed.
- **Flag-quality scorer / threshold tuner** — score a mask by residual Gaussianity
  (Anderson–Darling) + flagged-fraction + residual SK; grid-search the Pareto front when no
  ground truth exists. Reuse: `jansky.rfi`, `scipy.stats.anderson`.

## SETI / technosignatures

- **Injection–recovery benchmarking harness** — `(waterfall, drift grid, SNR grid, detector) →
  P_detect(SNR, drift)` matrix; the first public, reproducible SETI recovery benchmark. Reuse:
  `jansky.seti.drifting_tone/drift_search`.
- **CPU-only drift-search CLI** — `drift-search file.h5 --drift-min … --snr 8` on a small
  filterbank in <60 s, pure NumPy (turboSETI is heavy; hyperseti/BLISS need GPU/C++). Reuse:
  `jansky.seti.drift_search`, `jansky.formats.read_filterbank`. Data: **Voyager-1 GBT `.h5`
  (~50 MB) — has a *verified real* drifting tone**.
- **Paired setigen-injection / `jansky.seti`-detection scored benchmark** — close the loop
  setigen leaves open (it has no detector); CI-runnable sensitivity regression test.
- **Format-agnostic cadence-vetting pipeline** — extend `jansky.seti.cadence_detection` from a
  bool to a scored `CadenceVerdict` (per-scan S/N, zero-drift RFI-zone flag, audit trail);
  foundation for reproducing the BLC1 logic. Refs: Sheikh 2021 (Nat. Astron. 6, 352).

## Cross-cutting ecosystem gaps (from the GitHub/automation sweep)

These are domain-independent and partly addressed by *this* project's automation layer:

- **Airflow-on-Podman reproducible workflow template** for public radio data on a laptop —
  rootless, daemonless, offline, CPU-only. **Zero prior art** (CHIME/FRB uses Airflow + Docker +
  Swarm; Podman unused in radio astronomy). *This project's P5 is the first instance.*
- **Offline-by-default public-data wrapper** around CHIME `cfod` / ATNF psrcat / NRAO-ALMA /
  NASA LAMBDA (current tools fetch catalogues at runtime and assume connectivity).
- **A maintained "awesome-radio-astronomy" index** flagging maturity / CI / GPU / offline status
  — the meta-resource itself is missing.
- **CPU-only FRB/pulsar candidate on-ramp** composing the healthy CPU tools (PRESTO, riptide,
  sigpyproc3, `your`, fitburst) into one tested pipeline, avoiding the GPU-only/abandoned corner.
- **Reproducible, offline radio-ML tooling** with CI + pretrained models (FETCH/sclassifier/ClaRAN
  are GPU-assumed, untested, or abandoned).

## Station hardware methods (deferred from the station track)

- **LNA-first vs filter-first: measured NF/overload A/B** (deferred 2026-07-11, full plan
  preserved in `plans/77-station-lna-filter-ab.md`) — the folk trade-off nobody has published
  measured. Deferred because the station's sealed Discovery feed fixes its own ordering, so the
  measurement no longer informs any build decision; needs a dedicated ~$85 discrete-parts buy
  (SPF5189Z, standalone 1420 MHz filter, SAWbird+ H1) + the owned bench. Its on-sky Discovery-feed
  Tsys leg moved to `plans/78`. Still the quickest citable result in the station track
  (RNAAS-able, bench-only); re-run the novelty check before picking up.
