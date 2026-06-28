# 20 — Multi-decade VLBI flux variability of compact AGN (Astrogeo)

Status: 🚧 in progress (tooling)

## Context

The **Astrogeo VLBI image database** (Petrov; astrogeo.org/vlbi_images) is the largest public VLBI
archive — ~139k images of ~21k compact radio sources from decades of geodetic and absolute-astrometry
VLBI, almost all dual-band **S/X (2.3 / 8.4 GHz)**. Because each source is observed across many
sessions spanning years to decades, the per-source total-flux-density histories form **multi-decade,
parsec-scale light curves** — a variability regime complementary to our arcsec-scale VLASS three-epoch
slice (#13). The data are fully public over plain HTTP with **no authentication** (unlike CASDA, which
the Stokes-V slice #15 is blocked on).

This slice reuses the transient-survey variability statistics we already built and tested for VLASS
(`vlass.eta_metric` / `v_metric` / `variability_metrics` / `select_candidates` / `injection_recovery`)
and applies them to Astrogeo flux histories, plus the dual-band data give a per-epoch S/X **spectral
index** (`spectra.spectral_index`) so a source can be characterised as variable *and* by spectral
behaviour. The deliverable is a reproducible tool + a variability-ranked **candidate catalogue** with a
recover-a-known validation on a famously variable blazar — not a discovery claim.

## Reuse

- `vlass.variability_metrics` (η, V, m_debiased, χ², p), `vlass.select_candidates` (2-D log η–log V
  outlier cut), `vlass.injection_recovery` (data-driven completeness) — applied verbatim to the VLBI
  light curves.
- `spectra.spectral_index` for the two-point S/X index per epoch.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts).

## Deliverables

- `src/jansky_research/vlbi.py`:
  - `synthetic_lightcurves` — offline fixture: a steady population (η≈1) + an injected variable subset
    with known truth labels, dual-band, uneven epoch sampling, so tests recover the injected variables
    with low contamination.
  - `band_variability` — per-source `variability_metrics` on one band's flux history (reuses `vlass`).
  - `sx_index` — per-epoch / mean S/X spectral index (reuses `spectra.spectral_index`).
  - `select_variable` — the 2-D log η–log V candidate cut (reuses `vlass.select_candidates`), with a
    minimum-epoch gate.
  - `fetch_astrogeo` (network, `# pragma: no cover`) — per-source flux histories from Astrogeo.
  - `run(offline=...)` writing `results/vlbi_metrics.json`, an η–V candidate figure, and macros.
- `tests/test_vlbi.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for the Astrogeo VLBI database.

## Approach

1. **Tooling (this step).** Pure-NumPy metrics composing the `vlass`/`spectra` helpers, validated on a
   synthetic dual-band population (steady sources at η≈1; injected variables at high η and high V).
   Offline `run` recovers the injected variables at low contamination.
2. **Real-data fetch (next).** `run(offline=False)` pulls a curated source list's Astrogeo flux
   histories (S and X separately), assembles per-band light curves, applies the variability cut, and
   computes the mean S/X index. Honest systematics: VLBI total flux density depends on `(u,v)` coverage
   and resolved-out flux per session, so an apparent change can be structural, not intrinsic — guard
   with a per-session-quality / minimum-epoch cut and report it as a caveat.
3. **Recover-a-known validation.** Reproduce the known strong variability of a well-studied blazar
   (e.g. a CGRaBS/3C blazar with published VLBI flares) from its Astrogeo history — the slice's
   ground-truth check, like FK Com for #13.
4. **GATE-2** before write-up — candidates survive the `(u,v)`-coverage / resolution caveats and a
   literature cross-check.
5. **Write-up** as `papers/vlbi/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers ≳80% of injected variables with low false-positive contamination.
- η and V match the `vlass` implementations on a constant and a known-variable light curve.
- (Real-data, later) the blazar recover-a-known reproduces its published variability; candidates
  survive the VLBI caveats and a literature cross-check; GATE-2 sign-off before write-up.
