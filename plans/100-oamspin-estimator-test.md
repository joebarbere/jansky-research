# Plan 100 — `oamspin`: can the "twisted light" estimator see black-hole spin in EHT 2017 data?

Status: **plan, controls frozen 2026-10-17 before any image has been made from the EHT data.**
The 2017 M87\* uvfits tarball has been downloaded once to verify it (GATE-0 data check) and was
deleted unopened beyond its file listing. Changes after the first real run go in
`survey/oamspin-findings.md` with the reason (CLAUDE.md, slice pattern). Idea: `ideas.md` L12.

## Context

Tamburini, Thidé & Della Valle 2020, MNRAS 492, L22 (doi:10.1093/mnrasl/slz176;
arXiv:1904.07923v4, read in full 2026-10-17) report an orbital-angular-momentum (OAM) asymmetry
in the EHT 2017 images of M87\* and infer spin a = 0.90 ± 0.05 at i = 17° (abstract). Their
method (Sect. 2(ii)–(iii), SM Sect. 4.2, p. 8):

- Two intensity images one day apart (epoch 1: 5+6 April 2017; epoch 2: 10+11 April; Sect. 2(ii))
  are treated as two planes along the line of sight, separated by the Earth's motion.
- The phase P is recovered from the Transport-of-Intensity Equation, ∂I/∂z = −(1/2π)∇·(I∇P)
  (SM Eq. 11a), solved by finite differences on 250 × 250 images.
- The field √I·e^{iP} on "concentric circles of pixels" (SM p. 6) is Fourier transformed in
  azimuth; q is the ratio of the m = +1 to m = −1 amplitudes (SM p. 8; Table 1 caption).
- q₁ = 1.417 ± 0.049, q₂ = 1.369 ± 0.047, mean 1.393 ± 0.024 → a = 0.904 ± 0.046 via KERTAP
  ray-traced rings over a = 0.5–0.9 (Table 1; p. 4).

GATE-0 (ideas.md L12, PR #345): no independent re-test found among 148 citing records. Read by
hand for this plan, three properties of the estimator decide what a test must do:

1. **A real brightness has |c₊₁| = |c₋₁|.** Any q ≠ 1 comes from P alone.
2. **P is driven by I₂ − I₁.** Two identical images give P = 0 and q = 1 exactly, so q measures
   the day-to-day image difference. For a source the paper itself calls stable to a few per cent
   (Sect. 2(ii)), that difference is dominated by imaging noise and calibration, not propagation.
3. **P scales as 1/Δz, and Δz is a convention.** The paper says the phase "recurs modulo λ"
   (p. 4). If q changes when Δz is rescaled, q reports the convention.

Two further points from Table 1 (p. 4), checked by hand: q = 1.391 at (a = 0.9, i = 17°) and
1.392 at (a = 0.5, i = 46°), so **q alone does not separate spin from inclination** — the quoted
a follows from the adopted i; and "Δq = 0.0012 corresponding to Δa ≃ 0.06" does not match the
table's slope (~0.24 in q per unit a gives Δa ≈ 0.005; our arithmetic).

**This slice tests the estimator, not the paper's number.** The authors' code and input images
are not public (GATE-0), so a mismatch with 1.39 could not separate "claim wrong" from "different
pipeline". What we can establish is whether the stated estimator, applied to EHT-like images, can
tell a spinning black hole from a non-spinning one at all.

## The estimator (frozen; implemented once, then not tuned)

`oamspin.q_estimator(I1, I2, pixel_uas, dz)`:
1. Both images on a common 250 × 250 grid (pixel set from the image field of view), each
   normalised to unit total flux.
2. TIE (Eq. 11a) with the standard FFT Poisson solution: ∂I/∂z ≈ (I₂ − I₁)/Δz, I = (I₁+I₂)/2,
   solved for P with a Tikhonov floor of 10⁻³·max(I) on I (the usual TIE regulariser; stated,
   not tuned).
3. u = √I·e^{iP}. Ring radius r_pk = peak of the azimuthally averaged I; annulus r_pk ± 5 µas
   (≈ the 2017 ring width scale); for each pixel circle in the annulus, an azimuthal DFT of u
   with 360 samples (bilinear interpolation); amplitudes summed over circles.
4. q = |c₊₁| / |c₋₁|.

Δz: the Earth's displacement along the M87\* line of sight over one day from the JPL ephemeris
(astropy), as the paper describes — the "fiducial" Δz. Unit tests (offline): an identical pair
gives q = 1 to machine precision; a planted Laguerre–Gauss ℓ = +1 field gives q ≫ 1 and its
mirror q ≪ 1; a real brightness with P ≡ 0 gives q = 1.

## Imaging (frozen)

`ehtim` (PyPI, GPL-3.0; installs and imports on Python 3.12, checked 2026-10-17) with one stated
regularised maximum-likelihood recipe for all data and all simulations: 2017 Stokes-I uvfits
(CyVerse doi:10.25739/g85n-f134, ODC-PDDL), 120 µas field, 64-pixel image upsampled to 250,
Gaussian prior of 40 µas FWHM, regularisers {simple entropy, TV} with fixed weights, three
imaging rounds with blurring between. The recipe is written into `oamspin.image_uvfits` and used
unchanged for every input. If the EHT's own published imaging scripts are confirmed public and
runnable, they are added as a second, independent imaging path; that is decided at step 1 and
recorded, before step 3.

## Frozen controls and interpretation rules

All simulations use the **exact uv coverage, times and per-baseline noise** of the real 2017 day
they stand in for (`ehtim` observe_same), so the null shares the real data's sampling.

- **N0 — identical pair (analytic).** q = 1 exactly. A unit test; also run on each real image
  against itself as a pipeline check.
- **N1 — non-spinning null.** A static, non-rotating crescent: a thin ring with a brightness
  asymmetry I(φ) ∝ 1 + A cos(φ − φ₀), **no phase structure and no spin**. Its radius, A and φ₀
  are taken from our day-095 lo-band image (ring radius ~21 µas expected from the EHT's 42 µas
  diameter — from memory, checked against EHT 2019 Paper I at step 1) **before q is computed on
  any real image**; station-gain error amplitudes are taken from EHT 2019 Paper III with a
  locator at step 1. For each of the two day pairs and both bands, 200 realisations: the same
  model observed on both days with independent thermal noise and gain errors, each imaged with
  the frozen recipe, q measured with the fiducial Δz. Output: the null distribution of q.
- **N2 — simultaneous-band pair.** The 2017 hi and lo bands of the same day were recorded
  simultaneously, so there is no physical Δz between them. q on (lo, hi) of each day.
- **N3 — Δz convention.** q on each real day pair with Δz × {0.5, 2, 10}.
- **N4 — calibrator.** 3C 279 2017 (CyVerse, HTTP 200 per GATE-0) through the same chain:
  q is reported; it has no spin interpretation.
- **R1 — the measurement.** q on the real pairs (095, 096) and (100, 101), each band, fiducial Δz.
- **P1 — positive control (step 2, optional).** If a public Kerr ray-tracer gives a spinning
  image with the same ring size and asymmetry, run it through the N1 chain; report q. Not
  required for the verdict.

Interpretation, fixed now:
- **Estimator cannot see spin** if the N1 95% interval of q contains 1.393, *or* N2 gives
  |q − 1| at least as large as the R1 departures, *or* N3 moves q by more than 0.05.
- **Estimator distinguishes spin** only if N1 is consistent with q = 1 (95% interval excludes
  1.3), N2 gives |q − 1| < 0.05, N3 moves q by < 0.05, and R1 departs from the N1 interval.
- Anything else is reported as **inconclusive** with the numbers. No outcome is written as a
  measurement or refutation of M87\*'s spin; the paper says what the estimator can and cannot do.

False-failure check before the run: on the N1 chain with noise turned off, q must equal 1 to
within 0.01 (if not, the imaging recipe itself injects asymmetry, and that is reported first).

## Deliverables

- `src/jansky_research/oamspin.py`: `tie_phase`, `spiral_spectrum`, `q_estimator`,
  `crescent_model`, `earth_dz_one_day`, `image_uvfits` and `simulate_pair` (`# pragma` for the
  ehtim/network parts), `run/_figure/_write_macros/_main` with `oamspinSyn*`/`oamspinReal*`
  macros through `preserve_live_macros` / `write_results`.
- Tests to the 85% floor; `results/oamspin_*.json` with `source`; `survey/oamspin-findings.md`;
  `papers/oamspin/` RNAAS note (software-cites `jansky-research`).

## Amendment A1 (2026-10-17, after the first real image, before any q)

Imaging is the EHT Collaboration's released `eht-imaging_pipeline.py` (2019-D01-02, commit
80d230e, `ttype='fast'`) instead of our hand recipe, which failed on the first image; the
estimator centres the annulus on the ring (equal-radius criterion, one common shift per pair); the
fiducial Δz uses physical units at M87 and the real-vs-null comparison is repeated across a Δz
sweep (post hoc). Details and the evidence for each in `survey/oamspin-findings.md`.

## Order of work

1. Estimator + offline tests (N0, planted ℓ = ±1, P ≡ 0). Imaging-recipe check on a model with
   noise off (the false-failure check). Decide on the EHT-scripts second path.
2. N1 null distribution (simulations only). 3. N2, N3, N4, R1 on real data. 4. Findings, GATE-2
   science review, note. Any control added after step 2 is labelled post hoc.

## Risks

- **Our TIE is ours.** The note says it tests the method as published, not the authors' code.
- **Imaging dominates.** That is the hypothesis under test, and N1/N2 measure it directly.
- **Contested topic.** Neutral framing; Tamburini et al. are cited for the method and the claim.
- **Compute.** 800 imaging runs at minutes each on a laptop CPU: ~1–2 days, resumable.
