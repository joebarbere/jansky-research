# lptlin — findings

Plan 92. GATE 0 complete, measurement made, novelty pass returned — and the novelty pass
**reinterpreted the measurement**, 2026-08-21. Read the "what this actually measures" section
before quoting any number from here.

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

**One detection**, and the same epoch gives |V|/I = 5.6% (V = 14.04 mJy) from `lptv` — our
own measurement, not a published one. The total observed polarized fraction is 12.3%, inside
the physical bound (a `total_polarization` helper checks `sqrt(L²+V²)/I <= 1`; above 1 would
mean the photometry or leakage is wrong, not that the source is exotic).

**RETRACTED: "more linearly than circularly polarized".** An earlier draft of this file drew
that comparison from 10.9% against 5.6%. It is invalid. Faraday rotation depolarizes the
*linear* component in a band-averaged image and does **nothing** to the circular one, so the
two numbers are not on the same footing and cannot be ranked against each other. Corrected
below.

The three non-detections sit at 2.3–4.1 sigma. Their L/I values (3.5–8.0%) are **not**
measurements and must not be quoted as such; they are what the debiased estimator returns
below the detection threshold.

## What this actually measures: depolarization, not the source's linear fraction

The novelty pass supplied the missing number and it changes the interpretation. **The source's
rotation measure is published**: RM = +89.1 ± 0.1 rad m⁻² (SB55237) and +90.5 ± 0.1 (SB58609),
both from Wang et al. 2025 (arXiv:2411.16606) via RM synthesis on channelized data.

For the SB60804 band (744–1032 MHz), uniform band-averaging retains
`|sin(RM·Δλ²)/(RM·Δλ²)|` of the intrinsic linear polarization, with Δλ² = 0.0780 m².
Computed independently here rather than taken on trust:

| |RM| (rad m⁻²) | linear polarization retained |
|---|---|
| 24.3 | 50% |
| **89.1** | **8.9%** |
| 90.5 | 9.9% |
| 159.3 | 1.2% |

So an intrinsically 90% linearly polarized pulse would appear at **8.0–8.9%** in our image.
**We measured 10.9%.** The measurement is therefore consistent with a source that is
intrinsically ~90–100% linearly polarized and depolarized roughly tenfold by its own known RM
— which matches the *total* polarized fraction of 92 ± 3% that Wang et al. report for a
different epoch of this source.

Inverting naively gives an implied intrinsic fraction of 110–123%, i.e. unphysical, so the
uniform-weighting sinc model slightly over-predicts the depolarization here (real channel
weighting, flagging, and a possibly different RM at this epoch all enter). The agreement is
at the tens-of-percent level, which is all this method can support.

**Consequently 10.9% is not the linear fraction of ASKAP J1832−0911 and must never be quoted
as one.** It is a strong lower limit, and its real content is that the pipeline recovers
exactly the depolarized signal the published RM predicts.

A clean internal check on that reading: Caleb et al. 2024 report >90% linear for
ASKAP J1935+2148 at RM = +159.3 rad m⁻² in the same ASKAP band — where naive band-averaging
would retain 1.2%. They used channelized RM synthesis. The technique, not the source, is what
differs.

## Novelty: the gap is real, but our number does not fill it

Of nine LPTs with polarimetry in the literature, **seven already have a published linear
fraction** — GLEAM-X J1627 (88 ± 1%), GPM J1839−10 (~100%), ASKAP J1935+2148 (>90% bright,
~40% weak), ASKAP J183950.5−075635 (80% main, 90% interpulse), ASKAP J175534.9−252749.1
(>60%), ASKAP J174508.9−505149 (23–97%, elliptical), ILT J1101+5521 (51 ± 6%). Every one
comes from RM synthesis or a narrow sub-band.

**Only ASKAP J1832−0911 and CHIME J0630+25 lack a published linear fraction**, both having a
published RM. Wang et al. 2025 give this source a *total* polarized fraction (92 ± 3%) for one
epoch and two RMs from two others; their observation log has no polarization column at all,
and nothing whatever about SB60804. So the gap our target sits in is genuine — but a
band-averaged 10.9% does not fill it, because it is not the intrinsic quantity.

**The right way to fill it** is RM synthesis on channelized data for SB60804, if CASDA holds
per-channel images or visibilities for that block. That would either recover a ~90% intrinsic
fraction (confirming this reading) or return ~10.9% again, which would imply a much smaller RM
at this epoch than two months earlier — itself worth reporting.

## Other limits on this measurement

- **The EVPA is observed, not intrinsic** — no ionospheric or Galactic Faraday correction is
  applied, and RM cannot be measured from a single band-averaged image.
- **The leakage floor is assumed**, taken equal to `lptv`'s I→V figure of 0.6%. A measured
  per-field I→Q/U leakage would be better. At L/I = 10.9% the detection clears it by a factor
  of 18, so the assumption is not load-bearing for the detection — but it is for any future
  claim near the floor.
- One epoch per pulse; no independent re-imaging check.

## Next

1. ~~Novelty pass~~ — done; see above. The gap is real, our number does not fill it.
2. Check CASDA for channelized (per-channel or cube) Q/U products for SB60804. If they exist,
   RM synthesis turns this from a depolarization demonstration into the missing measurement.
   If they do not, the honest output is a methods note: band-averaged MFS Q/U cannot measure
   LPT linear fractions at Galactic-plane RMs, quantified, with this source as the worked
   example.
3. Open item the novelty pass could not close: a Galactic-foreground RM for this line of sight
   (l ≈ 21°, b ≈ −1°) from SPICE-RACS DR2 or the Hutschenreuter Faraday sky, which needs a
   catalogue cone-search rather than a text search.

## Is the RM-synthesis route actually open? Checked, 2026-08-21

Queried CASDA for what SB60804 holds at this position (38 products):

- **No channelized or cube Q/U products.** Only `cont.taylor.0` and `cont.taylor.1` in each
  of i/q/u/v, plus residual/weights/noise maps. So RM synthesis **from images is impossible**
  for this epoch — the depolarized band-average is all the image products can give.
- **The visibilities are there, and they are small.** Two beams covering this position,
  `scienceData.VAST_1824-06.SB60804.…beam26/27_averaged_cal.leakage.ms.tar`, at
  **~636 MB each** (`access_estsize` 636,380 kB) — already calibrated, channel-averaged and
  leakage-corrected.

So the route is open in principle and the data volume is not the obstacle. What it needs is
per-channel Q/U at one sky position, which for a *single* position does not require full
imaging: phase-rotate the visibilities to the source and average per channel, then run RM
synthesis on the resulting spectrum. That is a new capability for this repo (reading a
measurement set at all), not an increment of this slice, and it should be planned rather than
bolted on.

**Current honest state of plan 92:** a real detection, correctly interpreted as
depolarization-limited, with the path to the intrinsic value identified, costed (~1.3 GB for
two beams) and deliberately not started.
