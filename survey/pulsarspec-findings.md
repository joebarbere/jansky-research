# Findings — pulsar radio spectral indices (ATNF Pulsar Catalogue)

`jansky_research.pulsarspec` reproduces the well-known result that pulsars have **steep** radio
spectra, $S_\nu\propto\nu^\alpha$ with a mean two-frequency index near $-1.8$ (Maron et al. 2000;
Bates et al. 2013; Jankowski et al. 2018). It is catalogue-only and maximal-reuse: it computes the
$400\to1400$ MHz index with `spectra.spectral_index` from the ATNF Pulsar Catalogue's tabulated S400
and S1400 flux densities (VizieR `B/psr`).

## Result: the steep pulsar spectrum, reproduced

Of 2536 catalogued pulsars, **473** have both an S400 and an S1400 flux density. Their
$\alpha^{1400}_{400}$ distribution:

| quantity | value |
|---|---|
| mean $\alpha$ | **$-1.77$** |
| median $\alpha$ | $-1.87$ |
| scatter (std) | 0.75 |

This sits squarely in the literature range ($-1.4$ to $-1.8$ across Bates 2013 / Jankowski 2018 /
Maron 2000), reproducing the steep-spectrum result from public data with a $\sim$hundred-line tool.

## Value-add: millisecond vs normal pulsars

Splitting at the standard $P<30$ ms boundary, **43** millisecond pulsars have both fluxes:

| population | mean $\alpha$ |
|---|---|
| millisecond ($P<30$ ms) | $-1.75$ |
| normal | $-1.77$ |

The two are **indistinguishable** — millisecond pulsars are *not* significantly flatter than normal
pulsars, consistent with the literature (Kramer et al. 1998; Jankowski et al. 2018). The $\alpha$--
period relation likewise shows no strong trend.

## Honest assessment & caveats

- **Reproduction, not discovery** — a reproducibility/tooling contribution (a validated pulsar
  spectral-index statistic), on-brand for this repo.
- **Two-point index, selection-biased.** A single $\alpha$ from two frequencies misses spectral
  curvature and turnovers; many pulsars peak near $\sim$100--300 MHz and flatten/turn over below
  400 MHz, which a 400/1400 index cannot see (Jankowski 2018 fit multi-frequency models). The sample
  is also flux-limited (only the brighter, both-band-detected pulsars), which can bias the mean.
- **Catalogue flux scale.** S400/S1400 are heterogeneous literature values with $\sim$tens-of-percent
  uncertainties and intrinsic scintillation scatter; the broad std (0.75) reflects this as much as
  intrinsic spread.
