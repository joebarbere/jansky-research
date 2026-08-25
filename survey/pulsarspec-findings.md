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

### Resolved 2026-08-24 (real re-run)

`compare_subsamples` now computes the difference, its Welch standard error, the observed
significance and the smallest offset the pair of subsamples could resolve at 2 sigma; all four are
committed and macro-backed. The real leg reproduces the referee's reconstruction almost exactly:

| quantity | value |
|---|---|
| mean alpha, MSP / normal | -1.75 (scatter 0.76) / -1.77 (scatter 0.75) |
| N, MSP / normal | 43 / 430 |
| difference | **0.018 +/- 0.122** |
| observed | **0.15 sigma** |
| resolvable at 2 sigma | **0.24 in alpha** |
| parent catalogue | **2536** (sample is 473, i.e. 18.7%) |

The abstract and Results now state the bound rather than the agreement: an offset larger than 0.24
is excluded, one smaller is not constrained. The Data section gives the selection a denominator, and
the Discussion states the *direction* of the joint-detection bias (a steep source near the S400
limit has already dropped below the S1400 limit, so requiring both bands biases the index flat).
"over the whole catalogue" is now "over the both-band-detected subset".

Three tests added, and the one that matters is `test_offline_recovers_the_injected_msp_offset`:
the fixture's +0.2 offset was previously injected and never checked. At the fixture's size it is
recovered at >2 sigma with `resolvable` < 0.2, which is precisely the contrast with the real
catalogue, where 0.2 sits *below* the 0.24 resolution. The old `test_run_offline` window was ~28
sigma wide on an exact algebraic transform and could not fail.

## Full referee round (2026-08-25): MAJOR REVISION (light), 14 findings

Nothing is wrong arithmetically: the referee re-derived the entire chain independently from
VizieR `B/psr` TAP (2536 / 473 / −1.771 / 0.0184 ± 0.1220 / 0.2441) and every committed
number reproduces exactly. All five DOIs check clean. The revision is for two abstract-level
claims stronger than the numbers support and one missing check that moves the headline; every
fix is computable from the columns the pipeline already fetches.

**MAJORs:**
1. The selection is stated as a flux limit and is mostly observational coverage: 764 of 2536
   have any S400 vs 1676 with any S1400 — the binding constraint is that modern surveys don't
   observe at 400 MHz, and an absent S400 is *no published measurement*, not a non-detection.
   A coverage-driven selection plausibly biases steep — the opposite of the direction the
   paper asserts. Restate as coverage-dominated with a flux-limit component; give both
   single-band parent counts.
2. The named bias is never quantified, and the `dr20radio`-style completeness cut moves the
   mean −1.771 → −1.88 at α_min = −3 (2–5× the SE on the mean, comparable to the whole
   literature spread the paper "reproduces"). The null survives every variant (|diff| ≤ 0.066
   vs 2SE ≥ 0.22). Commit one completeness-limited subsample and quote "−1.77 raw, −1.88
   above the completeness limit". (The S1400-floor mirror is NOT a clean bias estimate — it
   conditions on the numerator of α — and must not be presented as one.)
3. "Offsets larger than 0.24 are excluded" conflates 2×SE with exclusion: an offset of
   exactly 2SE is detected with 50% probability; 90% power needs 0.40; what the data exclude
   is the complement of the CI [−0.221, +0.257]. Quote the CI, or label 0.24 the 50%-power
   threshold.
4. No flux-scale systematic in the error budget: a 10/20/30% inter-survey scale offset
   between the arms maps to Δα 0.076/0.146/0.209 — at 30% it is the size of the entire quoted
   2σ resolution. State a systematic floor; note the Welch SE assumes per-source independence
   ATNF's publication-clustered fluxes don't have.
5. The headline mean has no uncertainty (SE = 0.0345) and "reproducing"/"sits squarely in the
   literature range" is wrong at that precision — −1.77 is at the steep *edge* of the quoted
   −1.4 to −1.8 range. The recorded "recovers vs responds to" failure mode; downgrade the verb.
6. Sixteen rounded scalars are the whole evidence; `fetch_atnf` is untested and its choices
   invisible — including that 474 rows have both fluxes and J0540−6919 (S400 = 0.0) is
   silently dropped by `s400 > 0`. Commit `results/pulsarspec_sources.csv` (473 rows) +
   `n_rejected_nonpositive_flux`.

**MINOR/NIT:** per-arm joint-detection rates differ (43/337 = 12.8% MSP vs 430/~2199 = 19.6%
normal) and are uncommitted — the exclusion applies to the observed subsets; the 30 ms cut is
never varied (referee's sweep 10–100 ms: sign flips, null never close to breaking — commit
it, plus permutation p = 0.877 or the median difference 0.09); `test_run_offline`'s ±28σ
window cannot fail — assert within a few SE of the injected mean; "every function is tested"
is false for `fetch_atnf` ("every analysis function"); ddof=0/1 mixed in one sentence; the
sign of 0.018 (MSPs marginally *flatter*) never stated; README line 66 still carries the
pre-revision "MSPs ≈ normal pulsars"; no catalogue version/access date (B/psr is a frozen
2017-vintage snapshot well behind live psrcat).

Declined from the attack map: the σ-uncertainty on `resolvable` (0.24 ± 0.03 at n=43) — not
material, do not spend the revision on it.

**Status: fixes pending.** The single change: commit the completeness-limited subsample —
for a paper whose stated contribution is reproducibility, quantifying its own selection
function is the contribution.
