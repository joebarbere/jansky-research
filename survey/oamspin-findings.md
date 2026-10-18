# oamspin findings — can the twisted-light estimator see black-hole spin? (plan 100)

Running log, newest section last. Plan: `plans/100-oamspin-estimator-test.md` (controls frozen
2026-10-17, PR #346). Every number below is in a committed `results/oamspin_*.json` unless marked.

## Step 1 — estimator (2026-10-17)

`oamspin.q_estimator` implements the plan's estimator (TIE SM Eq. 11a by the FFT/Teague method,
annulus spiral spectrum, q = |c₊₁|/|c₋₁|). Offline tests: identical images give q = 1 exactly; a
lopsided real brightness gives q = 1; planted ℓ = ±1 modes are recovered with their sign; the TIE
inverts a known phase exactly for uniform I and within 15% for a lopsided I (the curl-free
approximation of the FFT method, stated); the phase scales as 1/Δz. No EHT data touched.

## Inputs fixed before the first image (2026-10-17)

The plan left these to step 1 "with a locator"; they are set here before any image is made from
EHT data and are not changed afterwards.

- **Ring size.** EHT 2019 Paper I (arXiv:1906.11238): ring diameter 42 ± 3 µas (abstract;
  Table 1). The N1 crescent's radius, asymmetry A and orientation φ₀ are then measured from our
  own day-095 lo-band image, as the plan states, before q is computed on any real image.
- **Non-closing systematics.** EHT 2019 Paper III (arXiv:1906.11240): 2% in amplitude and 1° in
  phase (abstract; Sect. 8.4.5). Added in quadrature to the thermal σ of every simulated
  visibility.
- **Station gains in the null.** Amplitude gains random per scan at 10%; station phases fully
  random. Paper III Sect. 6.1 (p. 11) calls 10% "a lower limit" for the phased-ALMA SEFD
  uncertainty and gives no single network-wide figure, so 10% is a stated choice, not a
  measurement. The imaging's gain-free closure terms make the
  result insensitive to this choice; the amplitude term carries it.
- **Data terms of the frozen recipe.** Closure phases (weight 100), log closure amplitudes (100),
  visibility amplitudes (10) — the closure quantities are immune to station gains (Paper III
  Sect. 7), the amplitudes fix the scale. Not tuned after this point.
- **Total flux** for normalisation: images are normalised to unit flux before q, so the flux
  constraint only stabilises imaging; 0.6 Jy (EHT 2017 compact flux — from memory, and immaterial
  to q).

## Plan amendment A1 — after the first real image, before any q (2026-10-17)

No q has been computed on any real image when this was written.

- **Our hand recipe failed on the first real image.** As frozen it returned the Gaussian prior
  unchanged (objective NaN: ehtim's TV term divides 0/0 where the prior underflows to zero);
  with that fixed and the co-located AA–AP / JC–SM baselines flagged (1.14 Jy vs ≤ 0.42 Jy on all
  other baselines, day 095 lo), it made a ring but fit amplitudes poorly (χ²_amp ≈ 48; no
  self-calibration). Tuning it further would be fitting an imaging recipe to M87\*.
- **Replacement: the EHT Collaboration's own released pipeline**, `eht-imaging_pipeline.py` from
  github.com/eventhorizontelescope/2019-D01-02 (EHT 2019 Paper IV; GPL-3.0 per its header; the
  repo has no licence file, so it is downloaded at run time, not vendored), pinned to commit
  80d230e76c08edd31548e7cfe61dc7cf94300780. One change: `ttype = 'nfft'` → `'fast'` (the NFFT C
  library is not installed; the transform type is numerical, not a modelling choice). On day 095
  lo: χ² ≈ 0.8 in closure phase, log closure amplitude and visibilities; 51 s per image; vs the
  exact `'direct'` transform (4.7 min) the normalised images correlate at 0.998 (max pixel
  difference 6%). Every image — real and simulated — goes through the same `'fast'` path. This
  removes our own imaging freedom rather than adding to it; the plan anticipated the EHT scripts
  as an imaging path (step 1).
- **Ring centring, an estimator detail the plan left implicit.** The annulus must be centred on
  the ring, and RML images carry no absolute position. `ring_center` picks the centre about which
  the per-azimuth peak radius is most uniform (brightness cannot bias it; the azimuthal-profile
  peak was biased 1.5 px toward the bright side of a strongly lopsided test crescent, the
  centroid 18 px). Both images of a pair are shifted by one common offset, so the day-to-day
  difference is unchanged and identical inputs still give q = 1 exactly.
- **Units of Eq. 11a are not stated in the paper.** P scales as (pixel length)²/Δz, so the
  convention sets P's magnitude. Fiducial: x = physical transverse length at M87 (D = 16.8 Mpc, the
  paper's value, p. 3), z = the Earth's line-of-sight displacement in one day. Because that
  convention is arbitrary, the real-vs-null comparison is also made at each Δz of a sweep — a
  unit-free comparison. **Post hoc**: added here, after A1, before any q.

## False-failure check — FAILED (2026-10-17; simulations only, results/oamspin_falsefailure.json)

The plan's gate: a noiseless, calibrated, non-spinning crescent (radius 19.97 µas, width 15.9 µas,
A = 0.535, φ₀ = −2.15 rad, all measured from our day-095 lo image before any q; results/
oamspin_crescent.json) observed on the exact day-095 and day-096 lo coverage and imaged with the
EHT pipeline must give |q − 1| ≤ 0.01 at the fiducial Δz. **It gives q = 2.31.** Across the Δz
sweep the same pair gives q = 1.05, 2.31, 1.26, 0.17 at ×0.5, 1, 2, 10, and 0.36, 1.79, **1.369**,
5.73, 0.59, 0.71, 0.97 at ×10¹², …, 10¹⁸. The only difference between the two "days" is their uv
coverage.

Why: in physical units (x at M87, z = the Earth's one-day line-of-sight motion, 5.5 × 10⁵ km for
95→96) the TIE phase has rms ~10¹⁶ rad — wrapped ~10¹⁵ times — so e^{iP} is effectively random per
pixel and q is a noise statistic. P reaches ~1 rad only near Δz × 10¹⁶; q approaches 1 (0.97) only
at ×10¹⁸, where P ~ 0.01 rad.

Reading, before any real-data q: **for a source with no spin and no noise, this estimator returns
values spanning the paper's (q₂ = 1.369 is reproduced exactly at Δz × 10¹⁴) depending on a unit
convention the paper does not state, driven by day-to-day uv-coverage changes alone.** The plan
says such a failure is reported first; it is. N1 (with noise) still runs, as frozen, to give the
null's spread.

Two corrections made while setting this up, both before any real q: the model's orientation was
mirrored by the ehtim/FITS round trip (φ₀ → −φ₀; fixed with a vertical flip, now −2.141 vs −2.15),
and the Δz sweep was extended to 10¹²–10¹⁸ (post hoc) to reach the P ~ 1 rad regime.

## N1 — non-spinning null with noise (2026-10-18; results/oamspin_n1.json)

800 realisations (200 per day pair × band) of the crescent above, observed on the exact 2017
coverage with thermal noise, 10% gains, random station phases and the 2% / 1° systematics, each
imaged by the EHT pipeline (A1); 17.4 h wall clock. The 95% interval of q **contains the paper's
1.393 at every Δz in the sweep**:

| Δz × fiducial | median q | 95% interval | fraction q > 1.3 |
|---|---|---|---|
| 0.5 | 1.078 | 0.189 – 6.261 | 0.404 |
| 1 | 1.023 | 0.143 – 6.115 | 0.375 |
| 2 | 1.004 | 0.160 – 6.434 | 0.364 |
| 10 | 1.018 | 0.153 – 6.699 | 0.376 |
| 1e+12 | 0.936 | 0.155 – 7.968 | 0.338 |
| 1e+13 | 1.057 | 0.240 – 5.840 | 0.365 |
| 1e+14 | 0.995 | 0.429 – 2.667 | 0.084 |
| 1e+15 | 1.000 | 0.516 – 1.770 | 0.048 |
| 1e+16 | 0.997 | 0.670 – 1.690 | 0.059 |
| 1e+17 | 0.970 | 0.178 – 3.770 | 0.287 |
| 1e+18 | 0.979 | 0.354 – 2.853 | 0.214 |

By the plan's frozen rule ("cannot see spin if the N1 95% interval of q contains 1.393"),
**the estimator cannot distinguish a non-spinning source from the published value at any Δz**.
The narrowest interval (×10¹⁵, where the TIE phase is ~1 rad) is still 0.52 – 1.77. The
real-data controls (N2–N4) and the real pairs (R1) follow and are reported regardless.

## Real data: R1, N2, N3, N4 (2026-10-18; results/oamspin_real.json)

All eight 2017 images (4 days × lo/hi) and four 3C 279 images through the A1 pipeline; 9.7 min.

**R1 — the real day pairs** (q at Δz × 0.5 / 1 / 2 / 10 / 10¹⁵):
- 95-96 lo: 0.60 / **2.02** / 0.77 / 0.29 / 1.22 (fiducial inside the N1 95% interval: True)
- 95-96 hi: 0.29 / **1.61** / 2.24 / 1.12 / 0.75 (fiducial inside the N1 95% interval: True)
- 100-101 lo: 0.84 / **0.94** / 0.98 / 7.72 / 0.93 (fiducial inside the N1 95% interval: True)
- 100-101 hi: 0.79 / **0.21** / 3.08 / 2.19 / 1.13 (fiducial inside the N1 95% interval: True)

The two bands of a pair were recorded simultaneously yet give q = 2.02 vs 1.62 (95–96) and
0.94 vs 0.21 (100–101) at the fiducial Δz. Nothing here is a stable ≈ 1.39.

**N2 — simultaneous lo/hi of the same day** (no physical Δz):
- day 95: 1.97 / **0.92** / 0.07 / 3.00 / 2.13
- day 96: 1.45 / **1.04** / 0.95 / 0.63 / 0.28
- day 100: 1.88 / **0.64** / 3.56 / 0.16 / 1.58
- day 101: 1.11 / **0.63** / 1.15 / 1.40 / 1.22

**N3 — Δz rescaling:** in every R1 row q moves by far more than the 0.05 the plan allows
(e.g. 95–96 lo: 0.60 → 2.02 → 0.77 → 0.29 at ×0.5 / 1 / 2 / 10).

**N4 — 3C 279** (a blazar jet, no ring; same chain, run as frozen although the M87 pipeline's
fixed flux/prior/field are not tuned for it):
- 95-96 lo: 12.75 / **0.43** / 1.32 / 1.12 / 1.32
- 100-101 lo: 0.39 / **0.75** / 0.40 / 1.29 / 1.14

## Verdict by the frozen rules

The plan's three "cannot see spin" criteria are each met independently: the N1 95% interval
contains 1.393 (at every Δz); N2 departures from 1 are at least as large as R1's; N3 moves q by
far more than 0.05. **The published twisted-light estimator, applied to EHT 2017 sampling, cannot
distinguish a spinning black hole from a non-spinning one.** This is a statement about the method
as published and as we implemented it (our TIE and our ring centring; the authors' code is not
public), not a measurement of M87\*'s spin, which is untouched by this result.

