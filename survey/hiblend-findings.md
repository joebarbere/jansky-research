# Findings: `hiblend` (plan 97). Does FAST-beam blending inflate HI fluxes?

Plan: `plans/97-hiblend-fashi-alfalfa.md`. The prediction, the five controls (C0–C4) and the
three outcomes were frozen in the plan on 2026-10-05, before any real cross-match. Every change
to a frozen value is logged here with its reason, in order.

## Step 1: frozen values verified against the sources (2026-10-05, before any real run)

Each value was read from the paper's own text: arXiv full text, page numbers from `pdftotext`.

| Frozen value | Source and locator | Verdict |
|---|---|---|
| FAST FWHM 2.9′ | FASHI DR2 (Zhang et al. 2026, arXiv:2606.31539), p. 4: "half-power beamwidth of 2.′9 at [1420 MHz]" and Table 1 "Beam size (FWHM) 2.′9 at 1420 MHz"; FASHI DR1 (arXiv:2312.06097) Table 1, same | **confirmed** |
| ALFA FWHM ≈ 3.5′ (3.3′ × 3.8′) | Giovanelli et al. 2005 (AJ 130, 2598; astro-ph/0508301), p. 8: "the beam sizes are 3.3′ along the azimuth direction and 3.8′ along the zenith"; p. 22: "an ALFA beam averaging 3.5′ width at half power" | **confirmed**: 3.5′ is the paper's own average; the ellipse is noted as a first-order simplification |
| Match radius 1.5′, \|Δv\| ≤ 100 km s⁻¹ | FASHI positional uncertainty ≈ 2.9′ / SNR (FASHI DR2 p. 5, after Koribalski et al. 2004): ≤ 0.6′ at SNR ≥ 5. ALFALFA HI centroids "on average good to only ∼20″ and their accuracy depends on SNR … For low SNR sources, offsets can exceed 1′" (Haynes et al. 2018, arXiv:1805.11499, Col. 3 notes and the matching advice). 1.5′ is ≥ 2.5σ of the combined error at SNR ≥ 5. FASHI DR2's own match used a 3′ × 3′ box with 100 km s⁻¹ (p. 8, §5.2) | **kept unchanged**. The tighter radius is deliberate: targets of interest have neighbours at 1–2′, which a 3′ box would confuse. The chance-match test (C4) checks it |

**The FASHI DR2 paper's own FASHI–ALFALFA comparison** (§5.2, pp. 8–9) cross-matches about
28,000 common sources. It reports that "Integrated flux densities agree well … for high-SNR
sources (SNR_FASHI ≳ 30). For faint sources with SNR_FASHI ≲ 20, however, the FASHI flux
measurements are systematically lower than their ALFALFA counterparts", and attributes this to
Eddington bias in ALFALFA. It does **not** split the comparison by neighbour separation, so the
plan's test is not anticipated. **Consequence for the frozen controls (no change needed):** the
S/N-dependent flux offset is exactly what C3 (the calibration on a disjoint isolated sample, with
S/N among its covariates) is there to remove. If S/N correlates with environment, C2 (the
spectral negative control) will show it.

**Recorded detail for the analysis.** Haynes et al. advise matching on optical counterpart (OC)
positions. This test matches HI to HI, so it uses the ALFALFA HI centroids (`RAJ2000`/`DEJ2000`,
not `RAO`/`DEO`). Neighbour separations are likewise HI-to-HI.

No frozen value changed in step 1.

## Step 2: three estimator biases found on synthetic data, fixed before any real run (2026-10-05)

The plan states the *model* (R_pred, beta, the samples); this step concerns how its inputs are
computed. All three changes were found by the planted-truth check on the offline fixture
(`synthetic_field`, blending planted at the model's own strength, so beta must come back as 1).
No real data had been opened.

1. **S_c taken from one survey's flux biases beta upward (3.4 for a planted 1).** The target's
   FAST flux sits in the denominator of the observed ratio *and* sets S_c in R_pred. When FAST
   measures low, both the ratio and R_pred rise, so the noise itself makes a slope. A
   noise-free check isolated a second, smaller bias of 1.24 from using an already-blended flux.
   **Change:** S_c is now the weighted geometric mean of the two surveys' fluxes, with weight
   w = σ_F² / (σ_A² + σ_F²), chosen so its noise is uncorrelated with the response noise to first
   order, minus the beam model's own blend contribution (`target_flux_estimate`).
   Noise-free result: 0.99.
2. **C3 covariates computed from blended quantities absorb the signal.** Calibrating on the
   catalogued flux and S/N of the target leaks blending into the correction applied to the
   primary sample. A first attempt to fix this computed S/N as snr × S_c / S_F, which divided by a
   measured flux again: beta came out 3.8 even with no blending planted. **Change:** C3 now uses the
   blend-corrected flux, and an S/N built from the catalogued flux *error* (a noise level, not a
   noise realization): snr_c = S_c / σ_F.

3. **Subtracting a blend that may not exist biases the null.** With fixes 1–2, the new unit test
   found beta = −0.32 ± 0.08 (3.8σ) on a field with *no* blending planted. The earlier
   8-realization mean of −0.19 had hinted at it, and I had read it as scatter. Cause: the C3
   covariates used the blend-*subtracted* flux, so under the null they shift with R_pred by
   construction. **Change:** the C3 covariates use the noise-decorrelated geometric-mean flux
   *without* blend subtraction (`subtract_blend=False`). R_pred keeps the subtraction, which only
   sets beta's scale if blending is real. If blending is present, the covariates carry some of it
   and the calibration absorbs a little, which biases beta toward 0 (the conservative direction).
   C1 measures how much.

**After all three changes** (16 synthetic realizations each, 3,000 targets):

| Scatter (dex) | Planted strength | beta (mean ± sd) | Mean quoted SE | beta_null (C2) |
|---|---|---|---|---|
| 0.08 | 1 | 0.992 ± 0.135 | 0.137 | −0.09 ± 0.21 |
| 0.08 | 0 | −0.014 ± 0.082 (none beyond 3σ) | 0.088 | −0.09 ± 0.21 |
| 0.15 | 1 | 0.891 ± 0.231 | 0.242 | −0.15 ± 0.40 |
| 0.15 | 0 | −0.014 ± 0.157 (none beyond 3σ) | 0.161 | −0.15 ± 0.40 |

The estimator is unbiased under the null, and the cluster-bootstrap SE matches the scatter across
realizations, so the quoted errors are honest. A planted beta of 1 comes back at 0.89–0.99,
slightly low at higher noise.

**Two implementation clarifications** of the frozen plan, also made before any real run:
- **Bootstrap units** are connected components of targets linked within 6′, not Tempel groups.
  Tempel covers only the SDSS footprint, and correlated targets are exactly those within 6′ of
  each other. This keeps the plan's intent (correlated targets are resampled together).
- **The beta fit is unweighted OLS.** The plan left the weighting open; unweighted adds no
  modelling assumption about the two surveys' error bars.

## Step 3: the real run (2026-10-06): AMBIGUOUS, the spectral null control (C2) failed

Run `scripts/hiblend_real.py --out .`, seed 97, 42 s; evidence in `results/hiblend_metrics.json`.

- **Inputs.** 156,269 FASHI DR2 sources and 31,500 ALFALFA α.100 detections (25,432 code 1),
  giving **22,889 matched targets**. The neighbour catalogue holds every FASHI source plus 3,801
  ALFALFA-only detections.
- **Samples.** 14,375 isolated (C3 calibration), 2,907 primary, 4,560 null (C2). The bootstrap
  resamples 2,330 and 3,684 clusters respectively.

| Gate (frozen order) | Result | Pass |
|---|---|---|
| C0 power, real geometry, synthetic fluxes | planted 1: β = 1.004 ± 0.038, detected in 20/20; planted 0: −0.007 ± 0.014, 0/20 | yes |
| C1 planted truth on real isolated targets | planted 1: β = **1.249** ± 0.038; planted 0: −0.012 ± 0.097 | yes (within the frozen ±0.3) |
| C4 match reliability | 90 of 22,889 matches survive a 10′ shift: 0.39% chance rate | yes |
| **C2 spectral null control** | **β_null = 0.256 ± 0.087 (2.95σ)**, against the frozen \|β_null\| < 2σ | **no** |
| Primary β | 0.385 ± 0.136 (2.83σ) | — |

**Outcome, by the frozen rule: ambiguous.** Neighbours that cannot blend (offset by more than
600 km s⁻¹) produce a positive slope against the R_pred blending *would* give them, of about the
same size as the primary sample's. So something tied to having a neighbour, but not to spectral
blending, moves the ALFALFA/FASHI flux ratio. The primary β (0.38 ± 0.14) is not distinguishable
from it: the difference is 0.13 ± 0.16. **No blending claim is made, and no limit either**,
because the null result the "not supported" outcome needs is contaminated by the same systematic.

**Recorded alongside, not changing the outcome:**
- **C1 recovered 1.25, not 1.0, on real isolated targets.** That is inside the frozen tolerance,
  and the synthetic fixture and C0 both give 1.00. The excess is specific to real fluxes, so the
  real calibration and flux distributions differ from the synthetic ones in a way that inflates
  β. That is a second sign the C3 calibration does not fully describe the real survey-to-survey
  differences.
- **The C3 coefficients are large and partly cancelling** (log S/N: −1.68, (log S/N)²: +0.37,
  log S: −0.38). A strongly S/N-dependent flux scale, as FASHI DR2 §5.2 reports (Eddington bias
  in ALFALFA at SNR ≲ 20), is being fitted on isolated targets and extrapolated to targets with
  neighbours.

**What could produce C2.** These are hypotheses for a follow-up, not tested here; per the plan
they would be post-hoc controls and must be reported as such:
1. **A residual S/N systematic.** R_pred ∝ S_n/S_c, so a large R_pred picks faint targets beside
   bright neighbours. If the C3 model leaves any S/N dependence in the ratio, it reappears as a
   slope against R_pred for *any* neighbour, blendable or not.
2. **Spatial confusion in the source-finding.** ALFALFA's flux extraction (a box in the cube) or
   its baselines may pick up a bright neighbour's emission at any velocity, for example through
   baseline ripple or sidelobes. That would be a non-spectral blending the C2 design assumes away.
3. **Real environment structure.** Targets with neighbours sit in denser regions, where one survey
   may have systematically different noise or RFI flagging.

Each predicts something different: (1) C2 depends on S_c and vanishes when the primary and null
samples are matched in target S/N; (2) it depends on the neighbour's brightness, not the
target's; (3) it depends on local density, not on any single neighbour. Running them would be a
new, post-hoc test set, and the plan requires it to be labelled that way.

## Step 4: post-hoc diagnostics of the C2 failure: predictions stated before running (2026-10-06)

**These are post-hoc.** They were designed after C2 failed, to locate its cause, and cannot
change step 3's outcome. Per CLAUDE.md (*freeze the controls before the first real run*), each
diagnostic's distinguishing prediction is written here **before** it is run.

All diagnostics use the step-3 samples and the same C3 calibration. "Residual" means the
calibrated log ratio, as in β.

- **D0. Out-of-sample calibration check.** Fit C3 on a random half of the isolated sample and
  look at the other half's residual against target S/N (5 bins).
  *If C3 is adequate:* every bin's median residual is within 2σ of 0.
  *If a residual S/N trend remains:* the faint bins deviate. Hypothesis (1) needs this.
- **D1. S/N-matched null.** Reweight the null sample to the primary sample's target-S/N
  distribution (10 bins) and refit β_null. Also fit β_null within S/N terciles.
  *(1) predicts* β_null shrinks toward 0 in matched or high-S/N subsets and is largest at low S/N.
  *(2)/(3) predict* it persists at every S/N.
- **D2. Which flux drives it.** On the null sample, regress the residual on R_pred plus
  log S_c (target flux) and log S_n,max (the brightest neighbour's flux).
  *(1) predicts* log S_c carries the effect, with a negative coefficient and the R_pred slope
  falling.
  *(2) predicts* log S_n,max carries it (positive), independent of S_c.
- **D3. Local density.** Count every catalogued HI source within 15′ (any velocity) of each
  target. Regress the null residual on R_pred plus log(1 + N_15).
  *(3) predicts* density carries the effect and the R_pred slope vanishes. Also checked in the
  isolated sample, which contains no neighbours within 6′ but varies in 15′ density: *(3)
  predicts* a density trend there too, while *(1)/(2)* predict none.

A diagnostic that matches none of its predictions is reported as such. None of these can turn
step 3 into a blending claim.

### Step 4 results (run 2026-10-06; `results/hiblend_diagnostics.json`, seed 97)

Judged against the predictions above, which were committed in `d35a374` before this ran.

**D0. Held-out calibration: inadequate, but not where (1) said.** Median residual of the
held-out isolated half by target S/N quintile (log S/N bin edges 1.48 / 2.04 / 2.17 / 2.30 /
2.47 / 4.41):

| quintile | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| median (dex) | −0.001 | −0.025 | −0.020 | −0.001 | **+0.039** |
| σ | −0.2 | −6.4 | −5.0 | −0.4 | **+8.9** |

C3 does not describe isolated targets it was not fitted on. The misfit reaches 0.039 dex, about
the size of the largest effect the test is looking for (≈0.045 dex per equal-flux neighbour). The
prediction said the *faint* bins would deviate; the faintest bin is in fact fine, and the worst
bin is the brightest, so the shape is not the one hypothesis (1) predicted.

**D1. S/N-matched null: the null slope goes away, as (1) predicted.** Reweighted to the primary
sample's S/N distribution, β_null = 0.117 ± 0.085 (1.4σ), against 0.256 ± 0.087 unmatched. By
null-sample S/N tercile: 0.380 ± 0.116 (3.3σ), 0.246 ± 0.136 (1.8σ), −0.030 ± 0.203 (−0.2σ),
low to high. Two limits: the lowest and highest terciles differ by only 1.8σ, and the
high-S/N tercile has the least lever arm (median R_null 0.0013 against 0.0080), so part of its
"vanishing" is lost power.

**D2. Which flux drives it: inconclusive.** With log S_c and log S_n,max added, the R_null slope
falls to 0.174 ± 0.095 (1.8σ). Both added coefficients have the signs their hypotheses predicted
(S_c −0.008 ± 0.008, S_n,max +0.007 ± 0.005), but neither is significant (−1.1σ, +1.3σ). The
three terms are collinear (R_null is close to a function of S_n/S_c), and this fit cannot
separate them.

**D3. Local density: rejected as the carrier.** In the null sample the density term is
−0.001 ± 0.009 and the R_null slope is unchanged (0.257 ± 0.081, 3.2σ). Hypothesis (3) predicted
the opposite on both counts. The isolated sample does show a density trend, −0.013 ± 0.004 per
unit log(1+N₁₅) (−3.3σ), which (3) predicted and (1)/(2) did not. But it has the wrong sign to
produce a positive β_null (a denser field gives a *lower* ALFALFA/FASHI ratio). It is a second,
small, environment-linked systematic, not the one that failed C2. (Median N₁₅: 2 isolated, 4
null, 4 primary.)

**What this supports.** The best-supported reading is hypothesis (1): the survey-to-survey flux
scale has S/N structure that C3's quadratic calibration does not capture (D0), and the null slope
disappears when the null sample is matched in S/N to the primary sample (D1). This is weaker
than it sounds. D0's misfit is not in the predicted place, D1's tercile contrast is 1.8σ, and D2
cannot separate the terms. Hypothesis (3) is rejected for the null slope. Hypothesis (2) is
neither supported nor excluded.

**What this does not do.** It does not rescue step 3. The calibration misfit (D0) is as large as
the signal, so a primary β measured with the frozen calibration is not interpretable at the
precision the test needs, whatever C2 had shown. Re-measuring the primary β with an S/N-matched
or re-calibrated sample would be a new analysis chosen after seeing these data. The honest
route is a new frozen protocol: a more flexible calibration (for example a spline in S/N) fixed
in advance and required to pass D0 on held-out isolated targets *before* C2 and β are rerun.
That is a decision for a follow-up plan, not something to do here.

## Step 5: plan 98, a recalibrated second attempt — stopped at C3′ (2026-10-06)

Plan 98 (`plans/98-hiblend-recalibrated.md`, frozen at `c8735cc`; code committed at `8657711`
before the run) replaced the quadratic S/N term with a linear spline (knots at the calibration
sample's log-S/N deciles), added log(1 + N₁₅), and put a new gate in front of everything. C3′
fits on isolated targets in even 2° RA strips, tests on odd strips, then the reverse. **Both
directions** must have every held-out S/N-decile median |Δ| < 0.010 dex and χ² p > 0.01. Seed 98;
`results/hiblend_v2_metrics.json`.

**Result: C3′ fails in both directions; the run stopped and no β is quoted.**

| | plan-97 form | plan-98 form |
|---|---|---|
| max \|held-out decile median\| (even→odd, odd→even), dex | 0.054, 0.052 | **0.015, 0.027** |
| χ² p (10 dof) | 0.0, 0.0 | **0.006, 0.001** |

The spline removes most of the misfit, and the middle eight deciles now sit within about 2σ of 0.
What remains is at the two ends of the S/N range, with the same signs in both directions:
- lowest decile: −0.015 (−2.7σ) and −0.027 (−3.5σ);
- highest decile: +0.009 (+2.5σ) and +0.015 (+3.5σ).

Those residuals are a third to a half of the ≈0.045 dex equal-flux blending signal, and at low S/N
they sit exactly where plan 97's null slope was concentrated (D1, lowest tercile).

**The frozen outcome applies: the survey-to-survey flux scale cannot be calibrated, with these
covariates and from isolated targets, to the precision this test needs.** Calibration is the
limiting systematic of the two-beam FASHI × ALFALFA blending test. This is the plan's
pre-stated negative and is reported as such. It also closes plan 97's question: its C2 failure
is at least partly a calibration failure, and a calibration good enough to separate the two
explanations has not been demonstrated.

**One more finding, from writing the tests (synthetic data only).** On the small synthetic sky
used in the unit tests, injected blends pushed 16–24% of C1's injected targets beyond the S/N
range of the calibration sample. There, both forms recovered a planted β = 1 erratically
(plan 97's form: −0.07 to 1.24; v2: −0.70 to 1.33). A linear spline extrapolates its last
segment, and a quadratic its curvature. On real data plan 97's C1 passed tightly (1.249 ± 0.038),
and plan 98 never reached C1. Any later attempt should still check, before running, what
fraction of injected targets lands outside the calibration range.

**Not done, deliberately.** Each further calibration form tried after this result would be
another fork chosen with C3′'s answer in hand. Two candidates are left as ideas, not runs:
- a per-region calibration (the misfit transfers poorly across RA strips, which a flux scale that
  varies by FASHI cube or ALFALFA drift strip would produce; untested);
- restricting the test to the middle S/N deciles, where C3′ is satisfied (this changes the
  estimand and the power, and needs C0 rerun).

Either would need its own frozen plan, written knowing that two calibration attempts have
already failed on these targets.
