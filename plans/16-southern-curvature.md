# 16 — Southern peaked-spectrum catalogue via GLEAM-X + RACS multi-band curvature

Status: 🚧 in progress (tooling + real GLEAM-X×RACS fetch + run done; recover-a-known + paper next)

## Context

The northern `peaked` slice (`plans/14`) selected GPS/CSS candidates from three frequencies (TGSS
150 MHz, NVSS 1.4 GHz, VLASS 3 GHz), but with a fundamental limitation: TGSS is shallow, so the low
-frequency point had to be used as an **upper limit**, giving only a *bound* on $\alpha_\mathrm{low}$
and no measured turnover frequency. This slice removes that limitation in the **southern** sky.

GLEAM-X DR2 (Ross et al. 2024) measures each source in **20 in-band sub-bands across 72–231 MHz**, so
the low-frequency spectral shape is *measured*, not bounded. Combined with the three RACS bands
(RACS-low 887.5 MHz, RACS-mid 1367.5 MHz, RACS-high 1655.5 MHz; McConnell/Hale/Duchesne 2020–2025),
a source has **up to ~23 flux points over a ×23 frequency span** — enough to fit a real log-parabola
SED and *measure* the turnover frequency. No published all-sky-south homogeneous GLEAM-X×all-three-RACS
peaked-spectrum catalogue exists (Callingham 2017 predates RACS; RadioSED II covers only 300 deg² of
Stripe 82), so a value-added scout catalogue is a genuine contribution (see `survey/new-findings-scan.md`).

## Reuse (the point of this slice)

- `spectra.spectral_index`, `spectra.crossmatch`, `spectra.fetch_survey` — indices + matching + fetch.
- `peaked.classify_sed` — peaked/steep/flat/inverted from two indices (the coarse fallback class).
- The parabolic-turnover idea from `peaked.peak_frequency`, generalised here to a **weighted N-point**
  log-parabola fit (`fit_log_parabola`) — the new capability the southern data enable.
- Image access (later, for compactness vetting): Data Central serves GLEAM-X / RACS Stokes-I cutouts
  (verified during the stokesv slice); VizieR serves the catalogues.

## Deliverables

- `src/jansky_research/southern.py` — tested logic:
  - `fit_log_parabola` — weighted least-squares $\log S = a\,x^2 + b\,x + c$ ($x=\log_{10}\nu$);
    returns turnover $\nu_\mathrm{pk}=10^{-b/2a}$, curvature $a$ (concave $a<0$ = a real peak),
    and a reduced-$\chi^2$ goodness; the measured-turnover upgrade over the upper-limit method.
  - `classify_curved` — peaked (concave, $\nu_\mathrm{pk}$ in band) / steep / flat / inverted / uss
    ($\alpha<-1.2$, candidate high-$z$ radio galaxy) from the fit + a two-point high index.
  - `find_peaked_south` — cross-match GLEAM-X (multi-band) × RACS (3-band), build each SED, fit,
    classify, return the candidate table.
  - `synthetic_field` — offline GLEAM-X(20-band)+RACS(3-band) fixture with injected peaked / steep /
    USS / flat SEDs.
  - `run(offline=...)`, `_figure`, `_write_macros`, `_main`.
  - network fetchers (`fetch_gleamx`, `fetch_racs_bands`) — VizieR, `# pragma: no cover`.
- `tests/test_southern.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** The weighted log-parabola fit + classification + 3-survey orchestration,
   validated on a synthetic field that recovers injected peaked sources (with a *measured* turnover)
   at low contamination, and separates USS and steep.
2. **Real data (next).** A GLEAM-X-footprint cone (Dec $-80°$ to $+30°$) cross-matching GLEAM-X DR2 +
   RACS-low/mid/high; resolution-mismatch guard (GLEAM-X ~45″ vs RACS ~10–25″ → single-RACS-component
   + compactness cut); ~10% cross-survey flux-scale caveat. Recover known southern GPS (Callingham
   subset) as a recover-a-known — this time *measuring* $\nu_\mathrm{pk}$, not bounding it.
3. **GATE-2** before write-up — candidates survive flux-scale, resolution, and turnover-fit caveats.
4. **Write-up** as `papers/southern/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- `fit_log_parabola` recovers an injected turnover frequency to within a few % on noisy synthetic SEDs.
- Offline `run` recovers ≳80% of injected peaked sources with low false-positive contamination, and
  flags injected USS sources correctly.
- (Real-data, later) recovers known southern GPS sources with measured $\nu_\mathrm{pk}$; GATE-2.
