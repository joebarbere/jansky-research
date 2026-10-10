# `dmzcal` findings (plan 99)

## Step 0 — provenance (2026-10-10, before any posterior)

The decisions, counts and caveats are recorded in `plans/99-dmzcal-coverage.md`, in the section
"Step 0 outcome and pre-posterior refinements". The evidence is
`results/dmzcal_provenance.json`, produced by `scripts/dmzcal_provenance.py` with FRBs/FRB
pinned at 996fcda and zdm at e0f985b.

Headline: the certification sample is 36 bursts (CHIME 19, DSA 11, ASKAP 4, MeerKAT 2). Four
facts found along the way are worth keeping:

- **zdm's out-of-the-box redshift state cannot be audited for provenance.** `load_state()`
  defaults to `HoffmannHalo25`, which is unpublished. Its documented p(z|DM) script imports
  `imp`, so it does not run on Python ≥ 3.12. Anyone who calls zdm today without choosing a
  state gets parameters whose fit sample they cannot know.
- **The published fit sample is wider than "CRAFT".** Hoffmann et al. 2025 also used the DMs of
  25 DSA-110 commissioning bursts (redshifts for only 3). The DSA bursts that remain available
  for certification come from the later Sharma+2024 sample.
- **The FRBs/FRB host table under-flags spectroscopic redshifts.** The `Spectrum` column is blank
  for all 19 CHIME/KKO gold hosts of Leung et al. 2025, which do have spectra. A naive
  "flagged-spectrum only" rule would have left 20 certification bursts and an UNDERPOWERED
  verdict caused by a metadata gap.
- **One redshift is transcribed wrongly.** FRB 20231201A is 0.119 in FRBs/FRB and 0.1119 in the
  primary table. This could be reported upstream.

## Steps 1–2 — `dmzcal.py`, C0 and C1 (2026-10-10, before any real posterior)

The module is `src/jansky_research/dmzcal.py`. It contains the E2 posterior, PIT, coverage,
the calibration rule, the controls, and the step-0 sample-assembly helpers. The tests in
`tests/test_dmzcal.py` cover 99% of the module and are offline. The evidence is
`results/dmzcal_controls.json`, from `scripts/dmzcal_controls.py` with N = 36 and M = 2000.

**E2 checks.**
- The tabulated mean of DM_EG matches ⟨DM_cosmic⟩ + ⟨DM_host⟩/(1+z) to within 1% at
  z = 0.05–3.5.
- The remaining ~1% at z ≥ 1 is the Macquart Δ⁻³ tail beyond DM = 20000. That tail is kept out
  of the likelihood rather than renormalised into it. A first version renormalised it and lost
  2% of the mean at z = 2.
- A synthetic null drawn from E2 itself gives coverage 0.674 / 0.950 at 68% / 95% (N = 6000,
  KS p = 0.08), so the discretisation does not miscalibrate the null.

**Correction made during development.** I first wrote that the Macquart shape *cannot* reach
⟨Δ⟩ = 1 for σ ≳ 4. That was a fixed bracket at C₀ = −50, not a property of the model. Solved in
log space with an extending bracket, it reaches ⟨Δ⟩ = 1 at any σ (C₀ ≈ −1900 at σ = 40). The z
grid still starts at 0.01, for the reasons now stated in the docstring: the prior holds
2×10⁻⁷ of its mass below z = 0.01, and no certification burst lies below z = 0.019.

**C0 passes.**
- The bands are the central 95% of the null coverage distribution: [0.528, 0.833] at 68% and
  [0.861, 1.000] at 95%.
- On an independent null set, the false-fail rate is **0.045** against a gate of 0.07, so no
  widening was needed.
- The KS leg rejects only 0.0075 of null samples. Almost all of the false-fail rate comes from
  the coverage bands.

**C1 fails for the host case.**

| planted misspecification | power at N = 36 | required |
|---|---|---|
| host mean ×2 (μ_host + log10 2) | **0.12** | ≥ 0.8 |
| F ×2 | 0.95 | (none stated) |

- Under the frozen rule, the real-run outcome is therefore fixed in advance. Coverage numbers
  are quoted, but no calibration verdict is given ("the sample cannot distinguish calibrated
  from mis-calibrated at this N").
- The rule does have power against excess scatter: an F ×2 world fails it 95% of the time.

**A post-hoc explanation that did not hold.** After a preliminary M = 200 run, I attributed the
low host power to E2's comoving-volume prior, whose median is z ≈ 2.4. At that redshift the
host term is diluted by 1/(1+z). A supplementary run with the prior capped at z ≤ 1 was
labelled post hoc and non-gating before it ran. It gave host power **0.15**, so the prior is
not the main cause.

A doubled host mean shifts log DM_host by 0.30 dex. That is less than the 0.42 dex σ_host
the model already allows, so at N = 36 the central intervals barely move.

With the prior capped at z ≤ 1, F ×2 power drops to 0.40. At low z the DM is
host-dominated, so cosmic scatter matters less.

**Implication to settle before step 3.**
- The rule is asymmetric. Its false-fail rate is controlled at 4.5%, so a *failure* on real
  data is informative. A *pass* is not, because it cannot exclude a host mean off by 2×.
- The frozen text withholds every verdict when C1 fails. A narrower reading would report a
  failure, but never "CALIBRATED". That reading is a change to a frozen control. The real data
  have not been opened, but it is still the author's choice, and the user should decide it
  explicitly.

## Step 3 — real run (first pass: commit 6226ce5)

The first pass gave E1 / E2 / E3 coverage of 0.528 / 0.778, 0.222 / 0.500 and 0.167 / 0.444,
all OVERCONFIDENT. The committed write-up of that pass contained several claims that the GATE-2
round-1 review refuted or found unsupported. They are **retracted** here, not silently
replaced. The current numbers are in the next section.

| retracted claim (6226ce5 findings / CHANGELOG) | why |
|---|---|
| "E2's comoving-volume prior … alone predicts much of" E2's failure | Re-scoring with the prior capped at z ≤ 1 or z ≤ 0.7 gives the same result: median PIT 0.027 / 0.028. The posterior mass above z = 1 is negligible at these DMs. This is the second mechanism I proposed for E2 that did not survive a check. |
| "Whatever drives it is far larger than … a doubled host mean" | The C1 power was measured on E2's own synthetic population (z to 4) and does not transfer to a real z < 0.5 sample. |
| "The E1–E2 gap measures what survey-selection modelling buys" | E1 also differs from E2 in its SFR evolution, luminosity function and z-prior shape. |
| "Consistent with in-sample optimism" (E1 0.74 / 0.95 on its fit sample) | The production arm is itself biased low (median PIT 0.28, KS p = 0.018). Its 14 DSA bursts entered the fit by DM only, so their z is effectively held out too. |
| "Six bursts carry most of E1's 95% shortfall" | 4 of the 8 misses came from those six, which is half. |
| "Cause not separable … a host DM larger than modelled, MW DM … predict the sign" | A selection-free check separates the DM model from the z distribution; see below. |
| "Held-out" for E2 and E3 | Nothing was fitted for them. |
| "Recover-a-known" for FRB 20240304B | In this repo the term means a planted truth. This burst is a post-hoc single case study. |
| "E1 reproduces their zdm-based z_Macquart ≈ 2.8" | Caleb et al. attribute 2.8 to the Macquart relation, not to zdm (their "Clues" section). |

## GATE-2 round 1 — corrections, post-hoc checks, current result (2026-10-10)

The review verdict was PASS-WITH-FIXES with no blockers. The reviewer edited nothing: the
result hashes and `git status` were unchanged.

**Data corrections, applied after the first run.** All are in `results/dmzcal_provenance.json`
under `corrections`, each with a locator. Every burst was verified against primary sources by
an independent pass.
- **FRB 20201124A:** discovered by CHIME/FRB, not MeerKAT (Lanman+2022, arXiv:2109.09254,
  abstract). Its z becomes 0.0979 (Fong+2021, arXiv:2106.11993, abstract); FRBs/FRB's 0.0982
  has no located source. Scored under the CHIME model, its E1 PIT moves from 0.015 to 0.165.
- **FRB 20210410D:** DM 578.78 rather than 575.0 (Caleb+2023, arXiv:2302.09754, burst table).
  It was found in the MeerTRAP **incoherent** beam, so its survey model is now
  `MeerTRAPincoherent`. Its PIT is unchanged at 0.011.
- **Confirmed, not changed:**
  - FRB 20240201A's z = 0.042729. The 0.047279 lead is a digit swap in Shannon+2024's own host
    paragraph; SDSS DR18 gives 0.04273.
  - All 11 DSA redshifts (Sharma+2024, Extended Data Table 1).
  - All 4 ASKAP ICS detections and their bands.
  - FRB 20180916B as a CHIME discovery.
- **Unresolved, recorded but not applied:**
  - Three DSA DMs differ from Connor+2024 by 0.6–1.25 pc cm⁻³.
  - Connor+2024 gives z = 0.0700 for FRB 20231120A, against Sharma's 0.0368.

**Current verdicts.** N = 36 (CHIME 20, DSA 11, ASKAP 4, MeerKAT 1). Frozen bands; asymmetric
rule.

| estimator | cov. 68% | cov. 95% (floor 0.861) | median PIT | PIT < 0.16 / > 0.84 | KS p | verdict |
|---|---|---|---|---|---|---|
| E1 zdm HoffmannEmin25 | 0.556 | **0.806** | 0.185 | 0.42 / 0.03 | 1e-4 | OVERCONFIDENT (one-sided) |
| E2 minimal, no selection | 0.222 | 0.500 | 0.027 | 0.75 / 0.03 | 6e-13 | OVERCONFIDENT (one-sided) |
| E3 fruitbat Batten2021 | 0.167 | 0.444 | 0.016 | 0.83 / 0.00 | 9e-20 | OVERCONFIDENT (one-sided) |

**Reading the verdicts.**
- **"OVERCONFIDENT" is the frozen label, but the failure is one-sided.** PITs pile up low
  (28 of 36 below 0.5 for E1), and almost none fall in the upper tail. The posteriors are
  displaced towards high z. They are not symmetrically too narrow. Widening the intervals would
  be the wrong fix.
- E1's 68% coverage is in band, so its failure rests on the 95% leg and the KS test.
- **E1 is the result.**
  - E2 is the project's own no-selection comparison arm and was expected to fail.
  - E3 is fruitbat as documented. It has no host term: its PDF path never reads
    `subtract_host` (fruitbat `_frb.py:334-385`), and it subtracts no halo. Adding 50 pc cm⁻³ of
    halo still fails (reviewer check, KS p = 1.7e-10).
  - The three are scored on the same bursts, so they are not independent confirmations.
- **Scope.** The verdict concerns this host-identified, localised population. It is not a
  verdict on these estimators applied to unlocalised FRBs, which is the use they exist for.

**Post-hoc checks.** All of these were chosen after the verdicts were seen, and none is a
verdict. They are in `results/dmzcal_metrics.json` → `post_hoc`.

- **The DM model is fine at the known redshift.** The PIT of DM_EG given z_true is
  P(DM' < DM | z_true), from each estimator's own p(DM|z). It is independent of the z prior
  and of any redshift selection.

  | estimator | cov. 68% / 95% | median | KS p |
  |---|---|---|---|
  | E1 | 0.69 / 0.94 | 0.52 | 0.78 |
  | E2 | 0.75 / 0.94 | 0.52 | 0.63 |

  So the miscalibration of p(z|DM) is not in p(DM|z), and therefore not in the host-DM,
  Milky Way DM or cosmic-scatter terms at the level this sample can see. It lies in the
  **redshift distribution**: either the population model's p(z) for detected bursts, or which
  bursts receive a spectroscopic host. This sample cannot separate those two.
- **The bias appears as the DM limit rises.** This is post hoc, with nested, small subsets.
  Hoffmann et al. (Sec. 2.3, main.tex l.220-221) name the host-identification bias: fainter
  high-z hosts go missing, which biases a localised sample low in z at fixed DM. They use z only
  below DM_obs − DM_ISM = 183 pc cm⁻³ for DSA. E1 by that same variable:

  | DM_obs − DM_ISM | N | median PIT | KS p |
  |---|---|---|---|
  | < 183 | 9 | 0.59 | 0.80 |
  | < 250 | 12 | 0.44 | 0.31 |
  | < 300 | 16 | 0.29 | 0.15 |
  | < 400 | 24 | 0.25 | 0.027 |
  | < 500 | 27 | 0.28 | 0.015 |
  | all | 36 | 0.19 | 1e-4 |

  That pattern is what host-identification incompleteness predicts. It does not show it,
  because N is small and Hoffmann's 183 was derived for the DSA commissioning sample, not for
  this mix. E2 and E3 fail even below 183, as expected for estimators without a selection-aware
  prior.
- **Halo DM.** E1 at DM_halo = 25 / 75 has coverage 0.50 / 0.83 and 0.51 / 0.83, with KS p of
  2e-6 and 5e-4. The failure does not depend on the frozen 50.
- **Planted truth through the E1 and E3 code paths** (500 replications each, with z drawn from
  each burst's own posterior). The rule passes 0.948 (E1) and 0.958 (E3) of the time, against
  0.955 expected from C0. E2's C2 gives 0.962. The PIT implementations in the two out-of-repo
  drivers are not biasing the result.
- **Host-magnitude stratification**, which the plan pre-stated as conditional, was **skipped**.
  Host magnitudes exist only for the 20 CHIME/Leung hosts, so no sample-wide split is possible.
- **Point accuracy** on the certification sample. Median Δz/(1+z) is +0.057 (E1), +0.122 (E2),
  +0.154 (E3) and +0.056 (E4, N = 35). NMAD is 0.10, 0.15, 0.12 and 0.13.
- **E1 per telescope** (descriptive):

  | telescope | N | cov. 68% / 95% | median PIT |
  |---|---|---|---|
  | CHIME | 20 | 0.60 / 0.90 | 0.25 |
  | DSA | 11 | 0.64 / 0.82 | 0.34 |
  | ASKAP | 4 | 0.25 / 0.50 | 0.03 |
  | MeerKAT | 1 | 0.00 / 0.00 | 0.01 |

  The reviewer found that E1 still fails KS without ASKAP and MeerKAT (N = 30, p = 0.0035), and
  still fails with CHIME dropped (N = 15). This was measured on the pre-correction labels.

**Known provenance path, recorded rather than removed.** zdm's survey files at e0f985b list
some certification bursts (e.g. DSA 20221113A and 20221116A). Survey-file bursts enter the
survey efficiency through a median of DM_halo + DM_G (zdm `survey.py:1345`, `:1388`), so a
held-out burst can nudge E1's selection function. The effect is expected to be small, but it is
not measured. The lists are in `dmzcal_e1_zdm.json` → `survey_file_frbs`.
E1 also uses the state's fixed DM_G (`sigmaHalo = sigmaDMG = 0`). Hoffmann's fit modelled
σ_halo = 15 and σ_ISM = DM_NE/2.

**Not run.**
- NE2001 versus YMW16, because `pygedm` does not build here.
- An explicit low-z selection model. The DM-limit table is the nearest proxy.

**FRB 20240304B** (post-hoc case study, excluded from every statistic).

| | median z | 68% interval | 95% interval | PIT |
|---|---|---|---|---|
| E1 | 2.75 | [2.04, 3.18] ✓ | [1.06, 3.50] ✓ | 0.19 |
| E2 | 3.16 | [2.60, 3.57] ✗ | [1.80, 3.86] ✓ | 0.06 |
| E3 | 2.42 | [2.16, 2.64] ✗ | [1.77, 2.83] ✓ | 0.15 |
| E4 | 2.52 (point) | — | — | — |

Caleb et al.'s Macquart-relation range is [1.628, 3.397] (95%). All three estimators
overpredict this burst too. About 235 pc cm⁻³ of its DM is Virgo (Caleb+, Supplementary S1.4),
which a DM-only estimate cannot know about.

**Status.** These are *applied* fixes. Round 2 of GATE-2 decides whether they hold.
