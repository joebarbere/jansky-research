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
  `MeerTRAPincoherent`. Its E1 PIT moves from 0.0034 to 0.0113 (corrected at round 2; round 1
  misreported this as "unchanged"). Both corrections move towards less failure.
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

- **The DM model at the known redshift: no failure detected, with limited power.** The PIT
  of DM_EG given z_true is P(DM' < DM | z_true), from each estimator's own p(DM|z). It does
  not depend on the z prior or on selection on redshift. It does depend on any selection that
  acts on DM or on host properties at fixed z: for example, a host-identifiable sample that
  favours massive hosts with larger DM_host. The CHIME gold sample is that kind of selection.

  | estimator | cov. 68% / 95% | median | KS p |
  |---|---|---|---|
  | E1 | 0.69 / 0.94 | 0.52 | 0.78 |
  | E2 | 0.75 / 0.94 | 0.52 | 0.63 |

  Its power at these 36 redshifts, from `results/dmzcal_dmz_power.json` (synthetic, post hoc,
  KS p < 0.01, 1000 replications):

  | true DM model | rejection rate |
  |---|---|
  | unmodified (null) | 0.016 |
  | host mean ×1.5 | 0.24 |
  | host mean ×2 | **0.81** |
  | host mean ÷1.5 | 0.26 |
  | F (cosmic scatter) ×2 | 0.014 |
  | unmodelled +30 pc cm⁻³ (MW ISM/halo) | 0.08 |

  So **a host-DM mean error of about ×2 or more is disfavoured.** Smaller host errors, Milky
  Way errors of tens of pc cm⁻³, and *any* cosmic-scatter error are not constrained: at these
  low redshifts the DM is host-dominated.

  Within those limits, the miscalibration of p(z|DM) is more plausibly in the **effective
  redshift distribution** than in p(DM|z). That distribution is either the population model's
  p(z) for detected bursts, or the selection of which bursts receive a spectroscopic host.
  This sample cannot separate those two.
- **E1's PIT against DM.** These are post hoc, small, independent bins. Hoffmann et al.
  (Sec. 2.3, main.tex l.220-221) name the host-identification bias: fainter high-z hosts go
  missing, which biases a localised sample low in z at fixed DM. They use z only below
  DM_obs − DM_ISM = 183 pc cm⁻³ for DSA. E1 by that same variable, un-nested:

  | DM_obs − DM_ISM | N | median PIT | PIT < 0.5 | KS p |
  |---|---|---|---|---|
  | < 183 | 9 | 0.59 | 44% | 0.80 |
  | 183–300 | 7 | 0.17 | 71% | 0.10 |
  | 300–500 | 11 | 0.21 | 91% | 0.005 |
  | ≥ 500 | 9 | 0.14 | 100% | 0.002 |

  **No failure is detected in the 9 bursts below 183, but 9 bursts have little power.** The
  monotonic trend of PIT with DM across all 36 is **not significant** (Spearman ρ = −0.31,
  p = 0.07).

  The pattern is what host-identification incompleteness predicts, but it does not show it.
  A p(z) that is wrong in its high-DM tail predicts the same thing. Hoffmann's 183 was also
  derived for the DSA commissioning sample, not for this mix. E2 and E3 fail in every bin, as
  expected for estimators without a selection-aware prior.

- **Halo DM.** E1 at DM_halo = 25 (N = 36) and 75 (N = 35; one burst has DM_EG ≤ 0 at 75) has
  coverage 0.50 / 0.83 and 0.51 / 0.83, with KS p of 2e-6 and 5e-4. The failure does not depend on the frozen 50.
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

  Descriptive subsets on the corrected labels (refreshed at round 2):

  | subset | N | cov. 68% / 95% | median PIT | KS p |
  |---|---|---|---|---|
  | CHIME + DSA | 31 | 0.61 / 0.87 | 0.28 | 0.002 |
  | DSA + ASKAP | 15 | 0.53 / 0.73 | 0.18 | 2e-5 |
  | CHIME | 20 | 0.60 / 0.90 | 0.25 | 0.035 |
  | non-CHIME | 16 | 0.50 / 0.69 | 0.16 | 1e-5 |

  Without ASKAP and MeerKAT, E1's coverage is in band but its location still fails KS.

- **Unresolved source conflicts as a sensitivity row.** These are not applied to the headline.
  - Connor+2024's z = 0.0700 for FRB 20231120A against Sharma's 0.0368. It moves that burst's
    E1 PIT from 0.002 to 0.010.
  - Three DSA DMs differ by 0.6–1.25 pc cm⁻³. They move PITs by ≤ 0.01.
  - With all of them applied, E1's coverage, median PIT (0.185) and verdict are unchanged
    (KS p = 1e-4).

**Survey-file provenance path: measured.** zdm's survey files at e0f985b list **all 11**
DSA certification bursts in `DSA.ecsv`, and CHIME's FRB 20180916B in decbin 3. Listed bursts
enter the survey efficiency through a median of DM_halo + DM_G (zdm `survey.py:1345`, `:1388`).

Re-running E1 with every certification burst pruned from every survey file used
(`dmzcal_e1_zdm.py --prune` → `results/dmzcal_e1_zdm_pruned.json`) moves every certification
PIT by ≤ 6×10⁻⁵ and leaves coverage, median and verdict unchanged. The path exists, and its
effect is negligible.

E1 also uses the state's fixed DM_G (`sigmaHalo = sigmaDMG = 0`). Hoffmann's fit modelled
σ_halo = 15 and σ_ISM = DM_NE/2.
zdm reports its own version as `0.0.0`; the commit (e0f985b) is what identifies it.

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

**GATE-2 round 2** (verdict PASS-WITH-FIXES, no blockers): 10 of 13 round-1 findings were
closed. Six new items (N1–N6) are addressed above:
- N1: the DM|z claim is narrowed, and its power is quoted.
- N2: the 20210410D before/after number is corrected.
- N3: "all 11" is stated, and the survey-file path is measured.
- N4: the alternative-source sensitivity row is added.
- N5: the DM table is un-nested, a trend test is added, and the CHANGELOG is reworded.
- N6: the subsets are refreshed.

**Status.** These are *applied* fixes. Round 3 decides whether they hold.
