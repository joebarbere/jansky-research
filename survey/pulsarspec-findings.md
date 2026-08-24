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

## Referee round on the style conversion (2026-08-24)

**Restyle findings, fixed.** (i) The Introduction's rhetorical question had been recast as a purpose
infinitive, "a ... tool **to recover** the steep mean and the millisecond--normal comparison", which
is factive: it presupposes the recovery succeeded. That is defensible for the steep-mean leg and not
for the other, which is a **null** --- one does not "recover" a null. (ii) An appositive collapse put
"squarely in the literature range" inside the parenthetical series "(median ..., scatter ...)",
lending the *scatter* an endorsement the Discussion explicitly withholds. (iii) "essentially
identical" was left with "the normal pulsars" as its nearest antecedent rather than the pair of means.

## The null has no sensitivity (pre-existing, NOT fixed here)

The abstract says the two subsamples are "indistinguishable, i.e. millisecond pulsars are **not
significantly flatter**". `results/pulsarspec_metrics.json` has eight fields and contains **no
dispersion for either subsample, no N for the normal sample, no standard error and no test
statistic**. The code computes the dispersions and throws them away: `run()` calls
`spectral_distribution` on each subsample, which returns `std`, and stores only `mean`. The null
therefore rests on a bare comparison of two rounded means, -1.75 vs -1.77.

Taking the sample scatter 0.75 as a stand-in for both subsamples, with N_msp = 43 and
N_normal = 473 - 43 = 430: SE on the difference ~ 0.12, so the observed 0.02 is **0.17 sigma**, and
the smallest offset this sample could resolve at 2 sigma is **~0.24 in alpha**. That number is the
paper's actual result and it is nowhere stated.

Two consequences, both instances of "a null result divides by sensitivity, not by sample size":

- The offline fixture injects a millisecond--normal offset of **+0.2** --- *smaller than the real
  sample's 2-sigma threshold*. The real analysis could not have detected the very offset the
  validation builds in.
- No test asserts the injected offset is recovered. `test_run_offline` checks only
  `-2.2 < mean_alpha < -1.4` on a mean whose SE is ~0.014: a window ~28 sigma wide on an exact
  algebraic transform, which cannot fail. Meanwhile Methods tells the reader the tool is "tested
  against a synthetic pulsar population with injected steep spectra **and a flatter-millisecond
  sub-population**", which reads as though the split leg is validated. It is not.

**Fix (needs a real re-run, VizieR confirmed reachable):** write `std_alpha_msp`,
`std_alpha_normal`, `n_normal` and a two-sample statistic into the metrics, and restate the result as
a limit --- "the means differ by 0.02 +/- 0.12; offsets larger than ~0.24 in alpha are excluded". That
is a stronger paper than the unquantified "indistinguishable", and it makes the agreement with
\citet{kramer1998} checkable. Add a test asserting the fixture's +0.2 offset is recovered.

**Also:** "a validated pulsar spectral-index statistic **over the whole catalogue**" describes 473 of
2536 pulsars, i.e. **18.7%**. The parent count is in no macro and no results file, and the
double-detection requirement is exactly the selection that biases a two-frequency index **flat** (a
steep source at the S400 limit has already dropped below the S1400 limit) --- the same mechanism
recorded for `dr20radio`. The direction of that bias is not mentioned, and it acts on the headline
mean.

**Citation note, do not "repair":** Crossref returns a single author for
`10.1016/0370-1573(91)90064-S` (Bhattacharya alone). That is an incomplete legacy Elsevier record;
the paper is Bhattacharya & van den Heuvel and the `.bib` entry is correct as it stands.
