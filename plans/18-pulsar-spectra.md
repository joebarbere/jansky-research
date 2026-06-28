# 18 — Pulsar radio spectral indices (ATNF catalogue)

Status: ✅ done (tooling + real ATNF run reproducing alpha~-1.8 + paper)

## Context

Pulsars have steep radio spectra: their flux density falls off as $S_\nu\propto\nu^\alpha$ with a mean
two-frequency spectral index $\alpha\approx-1.6$ to $-1.8$ (Maron et al. 2000; Bates et al. 2013;
Jankowski et al. 2018) --- much steeper than the typical $-0.7$ of synchrotron sources. This slice
**reproduces** that result from the public ATNF Pulsar Catalogue (Manchester et al. 2005), using its
tabulated 400 and 1400 MHz flux densities, and adds the value-add comparison of **millisecond vs
normal** pulsars and the $\alpha$--period relation.

It is catalogue-only and maximal-reuse: it composes `spectra.spectral_index` directly (the 400$\to$1400
MHz two-point index) --- no blocked archives, no new fetch toolchain.

## Reuse

- `spectra.spectral_index` — the two-point index + error (S400 at 0.4 GHz, S1400 at 1.4 GHz).
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, AASTeX paper).

## Deliverables

- `src/jansky_research/pulsarspec.py`:
  - `pulsar_alpha` — $\alpha^{1400}_{400}$ from the two flux densities (reuses `spectra.spectral_index`).
  - `is_millisecond` — $P<30$ ms classifier (the standard MSP cut).
  - `spectral_distribution` — n, mean, median, std of the $\alpha$ distribution.
  - `find_spectra` — per-pulsar $\alpha$ + MSP/normal class for those with both fluxes.
  - `synthetic_field` — synthetic pulsars with steep spectra (and a flatter-MSP sub-population).
  - `fetch_atnf` — ATNF catalogue (VizieR `B/psr`) P0/S400/S1400.
  - `run(offline=...)`, `_figure` ($\alpha$ histogram + $\alpha$ vs $P$), `_write_macros`, `_main`.
- `tests/test_pulsarspec.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** The $\alpha$ + classification + distribution, validated on a synthetic
   pulsar population that recovers an injected mean $\alpha$ and the MSP/normal split.
2. **Real data (next).** ATNF `B/psr`; the $\sim$470 pulsars with both S400 and S1400. **Reproduce**
   the steep mean ($\alpha\approx-1.8$) --- the recover-a-known --- and report the MSP-vs-normal means.
3. **GATE-2** before write-up --- the distribution survives the flux-scale and selection caveats.
4. **Write-up** as `papers/pulsarspec/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- `pulsar_alpha` / `is_millisecond` match hand-computed values.
- Offline `run` recovers the injected mean $\alpha$ and the MSP/normal split.
- (Real-data, later) reproduces the literature mean pulsar spectral index $\approx-1.8$; GATE-2.
