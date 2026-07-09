# Cross-slice audit — the fixed-threshold sensitivity confound (2026-07-09)

## The confound

A slice reports a raw trend "detection / occurrence / duty-cycle **rises** as [distance decreases]
or [source flux / brightness increases]", computed by applying a **fixed detection threshold** to a
signal whose amplitude depends on that same binning variable. Because the signal scales with the
variable (flux ∝ 1/d², or a brighter population clears a fixed SNR cut more often), the fixed
threshold *manufactures* an occurrence trend that is largely a detection-sensitivity artifact, not
intrinsic. The fix is a **sensitivity null model**: normalize the signal to a common reference value
of the binning variable, re-threshold, and report the corrected trend — the residual bounds the
intrinsic part.

This was found and fixed in two slices; this note records the repo-wide audit for the rest.

## Fixed (the confound genuinely bit)

- **`skr`** — raw SKR proximity duty-cycle near/far **3.33×**; 1/r² null (`distance_correct_flux`) →
  **1.39×**. Reframed from "proximity law" to a bounded near-null.
- **`junodam`** — raw Jovian-DAM range-quartile near/far **196×** (330× on the same p90 detector);
  1/r² null (`sensitivity_corrected_active`) → **2.2×**. The celebrated "~180×" is a
  threshold-amplified inverse-square visibility effect, not an intrinsic rise.

## Audit result for every other merged slice: NO further at-risk slices

Every other occurrence/detection-style slice either **already models the selection function** or is
**not a threshold-vs-covariate trend**. Verified 2026-07-09:

**Already handled (selection modelled explicitly):**
- **`svsbi`** — the beaming-fraction posterior is inferred with the detection threshold *inside the
  SBI simulator*: `forward_model` puts S = L/4πd² at each star's Gaia distance and detects at
  `DET_NSIGMA=5σ` AND above the leakage floor; SBC-validated. The threshold is in the simulator.
- **`fashienv`** — uses injection-validated **1/Vmax** for the flux-limit selection, and explicitly
  *dropped* its raw median-HI-vs-radius trend because it is selection-biased (deficient galaxies
  fall out of the flux-limited sample). Already made the honest move this audit looks for.
- **`vlass`**, **`drift`**, **`frblens`** — completeness/selection bounded by **injection–recovery**
  (inject a known flare / SETI signal / lensed pair, re-run the cut, measure the recovered fraction).
- **`ecallisto_census`** — occurrence *rate vs sunspot number*; the covariate does not set per-burst
  SNR, and the one selection effect (station coverage) is divided out (rate = N/coverage).
- **`stokesv`** / **`stokesv_discovery`** — a per-region leakage floor + an upper-limit census; no
  detection-rate-vs-distance trend is computed.

**Not applicable (amplitude / period / index / structure / limit / null, not an occurrence trend):**
`sourcecounts` (validated dN/dS — see below), `stacking` (stacked *mean flux* of sub-threshold
sources — immune to a per-source threshold by construction), `uss` (a flux-*scale* bias on spectral
index, not a detection threshold), `vlbi` (variability amplitude), `solarbursts` (exciter speed),
`offsets`, `pulsarspec`, `ppdot`, `hi`, `rmsky` / `rmstructure` / `rmdipole`, `peaked` / `southern`,
`windwaves` / `swaves` / `triangulate` / `type3synthesis`, `lpt` / `lptv`, `wdpulsar`, `frbwait` /
`frbperiod`, `torchdsp` / `torchfdmt`.

**Closest candidate, explicitly ruled out — `sourcecounts`:** its Euclidean-normalized differential
counts rise toward ~1 Jy (raw differential slope −1.91), and flux is the binning axis, so a flux
threshold does bias bins near the limit. But (i) it is a *validated reproduction* of the Hopkins
(2003) counts to **0.061 dex**; (ii) completeness is handled by cutting at **3.5 mJy** (above the
NVSS ~50%-complete point at ~2.5 mJy) with Eddington bias caveated; and (iii) the sensitivity-null
recipe — "renormalize the signal to a common reference of the binning variable" — is *structurally
inapplicable* when the signal **is** the binning variable. The correct tool there is a
completeness / Vmax function, which is a different fix and is already handled by the flux cut.

## Conclusion

The repo is methodologically consistent on this confound: it is corrected where it applies (`skr`,
`junodam`) and handled by an explicit selection model everywhere else it could arise. No further
sensitivity-null follow-up is warranted. Recorded so the cross-slice check need not be re-derived.
