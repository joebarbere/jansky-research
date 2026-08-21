# lptspec — findings

Plan 91. GATE 0 complete 2026-08-21 (recoverability **and** novelty). No image has been
fetched, and none should be until this gate is read.

## GATE 0, part 2 — novelty: CLEAR

No paper reports in-band ASKAP Taylor-term indices for any LPT pulse. Checked directly: an
arXiv full-text search for `taylor.1` returns no radio-astronomy hits, and none for
"long-period radio transient" AND "Taylor term"; the 2026 LPT review (arXiv:2601.10393)
discusses no such measurement. Every published LPT spectral index is either band-averaged /
cross-survey (GLEAM-X J1627 at alpha ~ -1.16 from MWA, Hurley-Walker+2022) or frequency
drift from dynamic spectra (Rose+2026) — different data, different measurement. One
near-miss worth naming: ASKAP J173608.2−321635 (arXiv:2405.13183) quotes an "in-band"
alpha = -3.1 +/- 0.2, but between THOR-GC at 1.23 GHz and VLITE at 339 MHz — two separate
telescopes, so cross-survey despite the label.

## GATE 0, part 1 — is a Taylor-term index recoverable at all?

Rose et al. (2026, Nature Astron. 10, 1166) report bursts that "drift in emission frequency".
ASKAP images continuum with a multi-frequency-synthesis Taylor expansion, so an in-band index
is available as `alpha = T1/T0` from `taylor.1` / `taylor.0` — images CASDA already serves
(98 `taylor.1` products confirmed at an LPT position, 2026-08-20). No sub-band re-imaging.

**The lever arm is the constraint, and it is short.** RACS-low spans 288 MHz at 887.5 MHz, so
the in-band frequency range is only about +/-16%. Fitting a slope across it gives a Taylor-1
image `sqrt(12) * nu0 / B ~ 10.7x` noisier than Taylor-0, hence

    sigma_alpha ~ 10.7 / (S/N)

so **S/N ~ 36 would be needed for sigma_alpha = 0.3** on that idealised calculation — against
the S/N ~ 5 that suffices for a detection. This is why the plan gated on an injection study
rather than assuming a detected pulse is a measurable spectrum.

**But the idealised number is roughly twice too optimistic, and the literature says so.**
Rashid et al. (2024, arXiv:2405.18978) pushed simulated point sources through MT-MFS — the
same Taylor-expansion synthesis ASKAPsoft uses — and found in-band indices "for SNR <~ 100
that have errors >~ 0.2, making them unreliable". The idealised formula gives 0.107 at S/N
100, so real MT-MFS carries about a **2x penalty** from deconvolution error, beam and
bandpass systematics, and spectra that are not exactly power laws. Applying it reproduces
the published behaviour (0.213 at S/N 100 against their ">~0.2") and moves the bar to
**S/N ~ 71 for sigma_alpha = 0.3**. That factor is a uGMRT simulation, not an ASKAP
calibration, so the real-cutout injection study must replace it with a measured value.

**Injection-recovery**, 40,000 trials per point, measuring bias as well as scatter (a ratio
with a noisy denominator is biased, and an error bar cannot show that):

| S/N | recovered (true −0.7) | scatter | bias |
|---|---|---|---|
| 300 | −0.700 | 0.036 | −0.000 |
| 100 | −0.700 | 0.107 | −0.000 |
| 50 | −0.701 | 0.214 | −0.001 |
| 20 | −0.702 | 0.538 | −0.002 |
| 8 | −0.711 | 1.371 | −0.011 |

The harness detects bias where a ratio estimator must have it — growing monotonically as S/N
falls — so "bias ~ 0 in our regime" is a measurement, not an artefact of a blind test. Below
S/N ~6 the `T0 > 5 sigma` guard itself selects upward fluctuations of T0 and *shrinks*
|alpha|: a selection effect, not noise.

## The verdict: 3 of 7 pulses are usable, and the motivating source is one of them

Across both legs of the `lptv` census (`results/lptspec_gate0.json`):

| pulse | leg | S/N (I) | sigma_alpha idealised | realistic | usable |
|---|---|---|---|---|---|
| ASKAP J1832−0911 / ASKAP-60804 | VAST | 262 | 0.04 | **0.08** | yes |
| ASKAP J183950.5−075635 / ASKAP-57929 | RACS | 214 | 0.05 | **0.10** | yes |
| **ASKAP J174508.9−505149 / ASKAP-20398** | RACS | 109 | 0.10 | **0.20** | **yes** |
| ASKAP J175534.9−252749.1 / ASKAP-47253 | VAST | 57 | 0.19 | 0.37 | no |
| ASKAP J183950.5−075635 / ASKAP-62646 | VAST | 36 | 0.30 | 0.59 | no |
| ASKAP J183950.5−075635 / ASKAP-62032 | VAST | 21 | 0.51 | 1.01 | no |
| ASKAP J165130.3−450520 / ASKAP-68311 | RACS | 15 | 0.71 | 1.42 | no |

Quoting the idealised column alone would have called five pulses usable. The published
MT-MFS behaviour says three. The difference is entirely the 2x penalty, and it is the
difference between a class-wide claim and three measurements.

The third row is the point: **ASKAP J174508.9−505149 is the Rose et al. source whose
frequency drift motivated this plan**, and its RACS pulse supports an index to +/-0.20 on the
published MT-MFS behaviour (+/-0.10 on the idealised floor) — enough to separate steep aged
plasma from flat coherent emission, though not to measure a subtle curvature. The measurement the paper's
qualitative claim invites is available in archival imaging.

## What this gate does NOT say

- **It is the idealised case**: Gaussian noise, no deconvolution error, no primary-beam or
  bandpass systematic, and a power-law source spectrum assumed across the band. Real
  Taylor-term alpha will be no better than this and probably worse, so a pulse failing here
  cannot be rescued by real data — but a pulse passing here is not thereby guaranteed.
- **Pulse dilution is not folded in.** `taylor.1` is fitted over the whole synthesis while
  the pulse occupies part of it, so the recovered alpha mixes the pulse with an empty
  integration. That is a separate systematic and must be quantified before any alpha is
  interpreted as intrinsic.
- A NaN guard was needed in the pulse selector: `float("nan")` does not raise and every NaN
  comparison is False, so a bare `< threshold` test admitted 319 non-detections as pulses on
  the first run. Fixed, and worth remembering — the same shape of bug would silently inflate
  any "N detections" count.

## Next

1. ~~Novelty pass~~ — clear, see above.
2. Fetch `taylor.0` and `taylor.1` cutouts for the three usable pulses (the CASDA
   machinery already exists in `stokesv.fetch_racs_cutout`; its filename mask is parameterised
   and `taylor.1` products are `restored.conv`), measure alpha, and quantify pulse dilution.
