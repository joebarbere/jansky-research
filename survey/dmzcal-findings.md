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
