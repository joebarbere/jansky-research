# Plan 97 — `hiblend`: does FAST-beam blending inflate HI fluxes? A two-beam test with ALFALFA

Status: **plan, controls frozen 2026-10-05 before any real cross-match was run.** Nothing below
may be changed after the first real run without recording the change and the reason in
`survey/hiblend-findings.md` (CLAUDE.md, slice pattern: *freeze the controls before the first
real run*).

## Context

`fashienv` (merged #318) left one question it could not answer: the group−field knee offset is
carried by sources whose flux aperture contains another galaxy, and those sources carry
+0.160 ± 0.007 dex more HI at fixed optical luminosity than isolated galaxies, against
+0.026 ± 0.008 when the neighbour lies just outside the aperture (`results/fashienv_metrics.json`,
`robustness.blending.group.hi_at_fixed_optical`). That is the signature of blending, but the note
could only say "interferometric HI for flagged members would decide it", because a real HI excess
in the closest physical pairs predicts the same aperture dependence.

Two single-dish surveys observed the same galaxies with **different beams**: FASHI with FAST
(FWHM ≈ 2.9′ at 1.4 GHz) and ALFALFA with Arecibo/ALFA (FWHM ≈ 3.5′; the beam is elliptical,
≈3.3′ × 3.8′). Blending depends on the beam; the galaxies' physics does not. So for the same
target, the ratio of the two surveys' fluxes must track the *difference* of the two beams'
response to its neighbours if blending is real, and must not if the HI excess is physical. No
interferometer is needed.

Both beam sizes are working values to be **verified with locators** before the first run
(CLAUDE.md: facts read from a source carry a locator): FAST — the FASHI DR1/DR2 papers; ALFA —
Giovanelli et al. 2005 (AJ 130, 2598) and Haynes et al. 2018 (ApJ 861, 49). Changing a beam
value after the first real run is a change to the frozen controls and must be logged.

## GATE-0 (run 2026-10-05)

- **Data, both public, no account.** FASHI DR2 Table 2 (cached at `data/fashi_dr2/`, already used
  by `fashienv`). ALFALFA α.100 (Haynes et al. 2018), VizieR `J/ApJ/861/49/table2`: 31,502
  sources with RA/Dec (HI and optical), `Vhel`, `W50`, `W20`, `HIflux`, `e_HIflux`, `SNR`, `rms`,
  `Dist`, `logMHI`, detection code `HI` (1 = high S/N, 2 = lower-S/N "priors"). Verified to load
  through the CfA VizieR mirror.
- **Prior work.** Jones et al. 2015 (MNRAS 449, 1856, "Spectroscopic confusion: its impact on
  current and future extragalactic HI surveys") quantified confusion **by simulation** for future
  surveys. ADS (2026-10-05) finds no empirical two-beam test on matched galaxies; the nearest use
  of both catalogues, arXiv:2606.25367, studies AGN and satellite HI statistics, not blending.
  The full-text check that the FASHI DR2 paper's own FASHI–ALFALFA flux comparison does not
  already split by neighbour separation is the first task of the slice (locator required).

## The prediction (frozen)

For a target with catalogued HI flux $S_c$ and neighbours $n$ at separation $s_n$ and velocity
offset $\Delta v_n$, a survey with Gaussian main beam of FWHM $\theta$ measures, to first order,

$$S_{\rm obs}(\theta) = S_c + \sum_n S_n\,B_\theta(s_n)\,V_n,\qquad
B_\theta(s)=\exp\!\left(-4\ln 2\,s^2/\theta^2\right),$$

where $V_n = 1$ if the neighbour's profile overlaps the target's integration window
($|\Delta v_n| < (W_{50,c} + W_{50,n})/2$) and 0 otherwise. The predicted two-survey log ratio is

$$R_{\rm pred} = \log_{10}\frac{S_c + \sum_n S_n B_{3.5'}(s_n) V_n}{S_c + \sum_n S_n B_{2.9'}(s_n) V_n}.$$

The differential response $B_{3.5'} - B_{2.9'}$ is **largest near $s \approx 2'$** (0.137), is
0.087 at 2.9′, and is negligible beyond 5′. For an equal-flux neighbour the excess is at most
≈0.045 dex. So the test is about separations of 1–4′, not a 2.9–3.5′ ring.

**The primary statistic** is the slope $\beta$ of the observed $\log_{10}(S_{\rm ALF}/S_{\rm
FASHI})$ regressed on $R_{\rm pred}$, after removing the survey-to-survey flux-scale calibration
(C3). Blending predicts $\beta > 0$, with $\beta \approx 1$ if the first-order model is right. A
physical HI excess predicts $\beta \approx 0$, because both surveys measure the same gas.

**Primary sample.** Targets matched between FASHI DR2 and ALFALFA α.100 (code 1 only) whose
HI-detected neighbour lies at $1' \le s \le 6'$ in either catalogue. $S_n$ is taken from FASHI
(the narrower beam) when the neighbour is in FASHI, otherwise from ALFALFA; the choice is frozen.
Neighbours with no HI detection (optical-only) are a secondary analysis, not part of the primary
test.

## Controls and pass criteria (frozen before the first real run)

- **C0. Power.** On synthetic matched catalogues with the real flux-error distributions and the
  real neighbour-separation distribution, the analysis must detect a true $\beta = 1$ at
  $\ge 3\sigma$ in at least 80% of realizations. If it cannot, the slice stops and reports the
  sample as underpowered. No real result is quoted.
- **C1. Planted truth through the full path.** Inject blends per the model into *isolated* real
  matched targets, through the production entry point, and recover $\beta = 1.0 \pm 0.3$. Inject
  nothing and recover $\beta$ consistent with 0.
- **C2. Spectral negative control.** Same separations, but neighbours with $|\Delta v| > 600$
  km s$^{-1}$, which cannot blend spectrally. $\beta_{\rm null}$ must be consistent with 0
  ($|\beta_{\rm null}| < 2\sigma$). A non-zero value means an environment-linked flux-scale
  systematic, and the primary result is then reported as ambiguous, not as blending.
- **C3. Calibration on a disjoint sample.** The survey-to-survey flux scale, and its dependence
  on S/N, $W_{50}$, flux and declination, is fitted **only** on isolated targets (no neighbour of
  either catalogue within 6′ in any velocity). It is applied unchanged to the neighbour sample.
  Isolated targets never enter $\beta$; neighbour targets never enter the calibration (CLAUDE.md:
  *an oracle that fed a fit may never certify the result*).
- **C4. Match reliability.** The cross-match radius and velocity tolerance are fixed in advance:
  1.5′ and 100 km s$^{-1}$, both to be confirmed against the catalogues' stated positional
  accuracies, with locators. The chance-match rate is measured by shifting one catalogue by 10′
  and must be below 1%.
- **Uncertainty.** Bootstrap over targets, resampling Tempel groups as units where targets
  share a group, so correlated neighbours are not counted as independent.

**Outcomes, stated in advance:**
- **Blending supported:** $\beta > 0$ at $\ge 3\sigma$, C1 passes, and C2 is consistent with 0.
  Report $\beta$ and what it implies for the `fashienv` group offset.
- **Blending not supported:** $\beta$ consistent with 0. Report the 95% upper limit on $\beta$
  and the blending fraction it excludes. `fashienv`'s group excess then becomes more likely to be
  physical, which is itself publishable.
- **Ambiguous:** C2 non-zero, or C0 or C1 fails. Report which control failed and make no claim.

## Deliverables

- `src/jansky_research/hiblend.py`, tested with a synthetic offline fixture (pure NumPy):
  - beam response and $R_{\rm pred}$;
  - the cross-match and its shifted-match chance rate;
  - the isolated calibration;
  - the $\beta$ fit with grouped bootstrap;
  - the injection and power harness.
- `scripts/hiblend_real.py` (network; `--out`).
- `results/hiblend_metrics.json` with `source` set, including C0–C4.
- `survey/hiblend-findings.md`, starting with this plan's frozen criteria.
- A short note, either an RNAAS note or a section added to `fashienv` (decided after the result,
  not before).

## Order of work

1. Verify both beams, the match tolerances and the FASHI DR2 paper's own ALFALFA comparison, each
   with a locator. Log any change to the frozen values.
2. Build the module and tests; pass C0 (power) on synthetic data. Stop here if underpowered.
3. Run C1 on real isolated targets, then C3, C4 and C2, **then** the primary $\beta$.
4. Write up whichever outcome occurred. A referee round follows, as for every slice.

## Risks

- **Small effect.** At most about 0.045 dex per equal-flux neighbour, against roughly 0.1 dex of
  survey-to-survey flux scatter. C0 decides whether the sample can see it at all.
- **Catalogued fluxes are not beam-weighted.** Both surveys integrate over a source box in the
  cube, so the Gaussian model is first-order. The sign and the radial shape of the prediction are
  robust to this; the expected $\beta \approx 1$ is not. Report $\beta$ as measured.
- **Unresolved pairs.** At $s \lesssim 1'$ both surveys merge a pair into one source, which is why
  the primary sample starts at 1′.
- **Flux-scale nonlinearity** between the surveys, if it correlates with environment through S/N.
  C2 and C3 exist for this.
