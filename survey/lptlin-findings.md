# lptlin — findings

Plan 92. GATE 0 complete and the first measurement made, 2026-08-21. Novelty pass in flight;
**no claim of a first should be made until it returns.**

## GATE 0 overturned this plan's own premise

Plan 92 was written with a scope reduction already applied: a CASDA query on 2026-08-20 found
only 2 EMU Q/U products at one LPT position and concluded "RACS carries no q/u at these
positions, only EMU does", so the uniform census "does not exist to be done".

**That generalised from a sample of one, and it is wrong.** Querying all sixteen catalogued
LPT positions:

- **14 of 16 have Q/U imaging in CASDA**, from VAST, EMU, and a dedicated ASKAP project
  named `ThreeLPTs` (440 Q/U products at ASKAP J1832−0911 alone, across 12+ SBIDs).
- **VAST does produce Stokes Q/U continuum images** (`image.q.VAST_1824-06.SB60182.cont.
  taylor.0.restored.conv.fits` and siblings) — the thing the plan assumed absent.
- The position originally sampled happens to be a RACS-era epoch, and RACS-low DR1 genuinely
  has no Q/U. The original finding was right about RACS and wrong as a generalisation.

The lesson is the cheap one: **a feasibility check on n = 1 is a feasibility guess.**

## Pulse-epoch coincidence: 4 of 7

A linear measurement is only meaningful in an observation that contains a pulse. Of the seven
adjudicated detections, four have both Q and U in the *same* scheduling block — all four VAST
epochs. The three RACS-era pulses (SB57929, SB20398, SB68311) have none.

## Why this survives the failure that killed plan 91

`lptspec` died because `taylor.1`, the MFS spectral term, cannot represent a source that
varies within the synthesis. This slice reads `taylor.0` — band-averaged flux — in Q and U,
exactly as `lptv` already does for I and V. A pulse present for a fraction *f* of the
integration is diluted to *f·I*, *f·Q*, *f·U* alike, so **dilution cancels in any ratio of
Stokes parameters from the same image**.

The cancellation is exact for the raw ratio and approximate after Ricean debiasing, because
noise does not dilute with signal: measured at *f* = 0.5, 0.2 and 0.05 the recovered fraction
is 0.02%, 0.12% and 2% low respectively. It errs **downward**, which cannot manufacture
polarization.

## The measurement

`results/lptlin_metrics.json`. Forced at the catalogue pixel in each Stokes image;
`L = sqrt(Q² + U²)` debiased by subtracting the noise in quadrature; leakage veto at 0.6% of
|I|, the same two-part criterion `lptv` applies to V.

| pulse | I (mJy) | L (mJy) | sigma | L/I | EVPA | detected |
|---|---|---|---|---|---|---|
| **ASKAP J1832−0911 / SB60804** | 249.97 | **27.21** | **33.2** | **10.9%** | 98° | **yes** |
| ASKAP J183950.5−075635 / SB62646 | 68.22 | 5.43 | 4.1 | 8.0% | 51° | no |
| ASKAP J175534.9−252749.1 / SB47253 | 24.93 | 0.87 | 3.1 | 3.5% | 0° | no |
| ASKAP J183950.5−075635 / SB62032 | 12.16 | 0.80 | 2.3 | 6.6% | 2° | no |

**One detection.** For ASKAP J1832−0911 the same epoch gives, from `lptv`, |V|/I = 5.6%
(V = 14.04 mJy), so this pulse is **more linearly than circularly polarized** — worth noting
for a class usually described as strongly circularly polarized. The total polarized fraction
is 12.3%, comfortably inside the physical bound (a `total_polarization` helper now checks
`sqrt(L²+V²)/I <= 1` explicitly; a value above 1 would mean the photometry or leakage is
wrong, not that the source is exotic).

The three non-detections sit at 2.3–4.1 sigma. Their L/I values (3.5–8.0%) are **not**
measurements and must not be quoted as such; they are what the debiased estimator returns
below the detection threshold.

## What this measurement is not

- **The 10.9% is a lower limit on the intrinsic linear fraction.** Q and U are band-averaged
  across 288 MHz at 887.5 MHz, so any appreciable rotation measure rotates the plane within
  the band and depolarizes it. J1832−0911 sits in the Galactic plane where |RM| can be large.
  The novelty pass is checking the |RM| that would cause ~50% band depolarization and any
  published RM for this line of sight.
- **The EVPA is observed, not intrinsic** — no ionospheric or Galactic Faraday correction is
  applied, and RM cannot be measured from a single band-averaged image.
- **The leakage floor is assumed**, taken equal to `lptv`'s I→V figure of 0.6%. A measured
  per-field I→Q/U leakage would be better. At L/I = 10.9% the detection clears it by a factor
  of 18, so the assumption is not load-bearing for the detection — but it is for any future
  claim near the floor.
- One epoch per pulse; no independent re-imaging check.

## Next

1. Novelty pass (in flight): whether J1832−0911's linear polarization is already published
   (its discovery paper, Wang et al. 2025 arXiv:2411.16606, reports X-ray and radio
   properties and may include polarimetry), and what linear fractions exist class-wide.
2. If novel, this is a short note: one detection plus three limits, with the depolarization
   caveat stated as a lower limit rather than buried.
