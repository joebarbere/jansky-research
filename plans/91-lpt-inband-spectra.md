# 91 — `lptspec`: in-band spectra of LPT pulses from the Taylor-term images

Status: 📋 planned 2026-08-20. Sourced from Rose et al. 2026 (Nature Astron. 10, 1166).

## Context

Rose et al. report that ASKAP J174508.9−505149's bursts "drift in emission frequency,
potentially due to a longer beat period". Frequency structure within a single pulse is
otherwise almost unmeasured for the LPT class from imaging surveys, because everyone quotes
band-averaged flux.

**Feasibility, measured 2026-08-20 (not assumed).** A CASDA ObsCore query at that position
returns 292 products, of which **98 carry `taylor.1`** — e.g.
`image.i.RACS_1750-50.SB38330.cont.taylor.1.restored.conv.fits`. ASKAP continuum imaging is
a multi-frequency-synthesis Taylor expansion, so `taylor.1 / taylor.0` at a pixel is a direct
in-band spectral-index estimate across the ~288 MHz band. No sub-band re-imaging is needed,
and `stokesv.fetch_racs_cutout` already stages exactly these files —
`_racs_science_mask` matches `restored`+`conv` and is parameterised by Stokes, so the only
change is selecting the Taylor term.

The repo holds seven adjudicated LPT pulse detections (three RACS, four VAST), including
J1745−5051 at *I* = 21.6 mJy / 21.6σ in *V* — the brightest, and the same source Rose et al.
characterise.

## The thing that will make this wrong

**Taylor-term spectral indices are unreliable at low S/N and this is the whole ballgame.**
`taylor.1/taylor.0` is a ratio of two noisy quantities; its error blows up non-linearly as
S/N falls, and the estimator is biased, not merely noisy. Published guidance puts usable
α at S/N of order 50–100 per pixel, and most of these detections sit well below that.

GATE 0 therefore requires an **injection study before any science claim**: inject sources of
known α across the observed S/N range into real RACS Taylor-term cutouts, recover α, and
publish the recovered-vs-injected relation. If α is only recoverable for one or two of the
seven detections, that is the result — quote those and give honest limits for the rest. Do
not report an α for a pulse the injection study says is unrecoverable, and do not let a
band-averaged non-detection masquerade as a flat spectrum.

Second systematic: `taylor.1` is fitted over the full synthesis, while a pulse lasts a
fraction of the snapshot. The recovered α therefore mixes the pulse with the (empty) rest of
the integration. Quantify this dilution before interpreting any α as intrinsic.

## Deliverables

- `src/jansky_research/lptspec.py` — Taylor-term α estimation with per-pixel error
  propagation, the injection-recovery harness, and an offline synthetic fixture.
- `results/lptspec_metrics.json` — per detection: taylor.0 flux, taylor.1, α with
  uncertainty, injection-derived recoverability flag, dilution factor.
- `survey/lptspec-findings.md` + figure: recovered vs injected α as a function of S/N, with
  the seven detections marked at their measured S/N.

## GATE 0

1. Injection study first, as above; if nothing is recoverable, stop and record the negative.
2. Novelty: check whether any LPT discovery paper already quotes in-band α from ASKAP
   Taylor terms (several report band-averaged spectral indices from separate surveys, which
   is a different measurement).
3. Confirm the pulse-dilution correction is computable from the committed `duration_s` and
   the published pulse widths; if not, scope the claim to "band-averaged over the snapshot".

## Related

`lptv` (the detections and the CASDA machinery), `stokesv` (`fetch_racs_cutout`,
`_racs_science_mask`), plan 90 (`lptduty`, same epoch table).
