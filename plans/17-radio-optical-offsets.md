# 17 — Radio–optical position offsets of AGN (Gaia DR3 × ICRF3)

Status: ✅ done (tooling + real ICRF3×Gaia run reproducing the 24× excess + paper)

## Context

The radio (VLBI) and optical (Gaia) positions of the same AGN do not perfectly coincide. The
**radio–optical offset** is real and physical: the VLBI position marks the synchrotron self-absorbed
jet base (core), while the Gaia optical position can be pulled toward optical jet/host structure or
the accretion disk. Across the AGN population the *normalised* offset (separation ÷ combined position
error) shows a heavy tail far beyond the Rayleigh expectation for pure Gaussian astrometric noise —
a well-established result (Mignard et al. 2016; Petrov & Kovalev 2017; Kovalev et al. 2017;
Lindegren et al. 2018; Petrov et al. 2019). This slice **reproduces** that result with a small,
tested tool, and adds a reproducible offset catalogue.

It is catalogue-only (no images, no blocked archives): **ICRF3** (the VLBI radio reference frame;
Charlot et al. 2020) cross-matched to **Gaia DR3** (Gaia Collaboration 2022). It reuses
`spectra.crossmatch` and the project's fetch/figure/macro conventions.

## Reuse

- `spectra.crossmatch` — ICRF3 ↔ Gaia positional matching.
- the CDS **XMatch** service (`astroquery.xmatch`) — efficient server-side ICRF3 → Gaia DR3 match.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, AASTeX paper).

## Deliverables

- `src/jansky_research/offsets.py`:
  - `radio_optical_offset` — separation (mas) + position angle from radio→optical, given two positions.
  - `normalised_offset` — $X=\sqrt{(\Delta\alpha^*/\sigma_{\alpha})^2+(\Delta\delta/\sigma_\delta)^2}$
    with $\sigma^2=\sigma_\mathrm{radio}^2+\sigma_\mathrm{Gaia}^2$ — the standard significance measure.
  - `offset_statistics` — n, median offset, the fraction with $X>3$, and the Rayleigh expectation
    (the excess tail is the reproduced result).
  - `synthetic_field` — ICRF3-like + Gaia-like positions where most offsets are pure Gaussian noise
    and an injected minority carry a real structural offset; tests recover the excess tail.
  - `fetch_icrf3_gaia` — ICRF3 (VizieR) XMatched to Gaia DR3 (positions + errors + G mag).
  - `run(offline=...)`, `_figure` (the $X$ distribution vs Rayleigh), `_write_macros`, `_main`.
- `tests/test_offsets.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** The offset + normalised-offset + statistics, validated on a synthetic
   field that recovers an injected excess tail and rejects a pure-noise field.
2. **Real data (next).** ICRF3 × Gaia DR3 all-sky; compute the $X$ distribution; **reproduce the
   excess tail** ($X>3$ fraction $\gg$ the ~1.1% Rayleigh expectation) — the recover-a-known. A light
   value-add: how the offset excess varies with Gaia $G$ magnitude / source brightness.
3. **GATE-2** before write-up — the excess survives the error-model and selection caveats.
4. **Write-up** as `papers/offsets/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- `radio_optical_offset` / `normalised_offset` match hand-computed values on textbook inputs.
- Offline `run` recovers the injected excess-tail fraction and reports the Rayleigh baseline.
- (Real-data, later) reproduces the literature radio–optical offset excess; GATE-2 before write-up.
