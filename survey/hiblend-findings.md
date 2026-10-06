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
