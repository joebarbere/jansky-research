# Findings — sub-threshold radio stacking with injection-recovery (SDSS quasars × VLASS-SE)

`jansky_research.stacking` measures the *average* radio flux of an optically-selected population whose
members are individually below the radio detection limit, by image-plane stacking — with the
credibility-critical step of **injection-recovery** bias calibration (White et al. 2007; Karim et al.
2011). Catalogue + cutout level, reusing the project's verified VLASS CADC-SODA path.

## The two access problems solved

1. **VLASS Single-Epoch is on CADC** (`.se.` products alongside `.ql.`). SE is self-calibrated and
   deeper-CLEANed — the right substrate, since Quick-Look's residual-vs-restored flux bias makes raw
   QL stacks of sub-threshold sources unreliable.
2. **Efficient cutouts.** `Cadc.get_image_list(query, pos, radius)` returns **server-side cutout**
   URLs; the filtered SE Stokes-I `tt0` URL downloads a ~73×73 stamp (~0.18 MB) in ~2.7 s — so a
   few-hundred-source stack is tractable (~7 s/source incl. the query) instead of downloading
   multi-GB full tiles.

## Real run: the mean radio flux of optical quasars

> **SUPERSEDED (2026-08-25).** The narrative below (and the two binned tables after it) documents
> an *earlier* run — N=279, "2.5° cone", searched-peak photometry, the identity injection ratio —
> that does not match the committed evidence (N=236, radius 3°), and its estimator carries both
> round-9 blockers (searched peak; identity calibration). It is kept as the historical record;
> the current numbers are in the RESOLVED section at the end of this file and in
> `results/stacking_metrics.json`.

Target: **SDSS DR16 quasars** (Lyke et al. 2020; VizieR `VII/289/dr16q`) over a 2.5° cone at RA 180°,
Dec +25° (SDSS × VLASS overlap). Of 300 tried, **279** had a VLASS-SE 3 GHz Stokes-I image (93%) and
entered the median stack:

| quantity | value |
|---|---|
| sources stacked | **279** |
| stacked central peak | 45.4 µJy/beam (0.0454 mJy/beam) |
| annulus RMS | 9.2 µJy/beam |
| **stacked SNR** | **4.9** |
| injection-recovery ratio | **1.00** |
| **de-biased mean flux** | **45.4 µJy/beam** |

So the median stack of individually-undetected SDSS quasars yields a **~4.9σ central source at
~45 µJy/beam** — far below the VLASS single-source limit (~0.7 mJy QL, ~0.36 mJy SE), recovered only
by stacking. The mean radio flux of $\sim$tens of µJy is consistent with the radio-quiet-quasar
population. The injection-recovery ratio of **1.00** is itself a result: the VLASS-SE flux scale is
**unbiased** for these centred sub-threshold sources (unlike Quick-Look would be) — vindicating the
choice of SE and confirming the de-biasing step is honest, not a fudge.

## Magnitude-binned: the radio--optical trend recovered

Binning the sample into three equal-count bins of SDSS $i$-band magnitude and stacking each (with its
own injection-recovery) turns the single number into a **trend** (279 quasars, 93 per bin):

| $i$-band bin | median $i$ | mean radio flux | SNR |
|---|---|---|---|
| bright | 18.81 | **81.4 µJy/beam** | 4.9 |
| mid | 19.95 | 60.4 µJy/beam | 4.0 |
| faint | 20.92 | 43.7 µJy/beam | 2.9 |

The mean radio flux **rises monotonically with optical brightness** — the optically-brightest third is
$\sim$1.9$\times$ radio-brighter than the faintest. This is the expected radio--optical luminosity
correlation (more optically-luminous quasars are radio-brighter on average), recovered by stacking from
*individually-undetected* sources, with each bin separately bias-calibrated (injection-recovery ratio
$\approx$1 throughout). The faintest bin is only 2.9$\sigma$, so the trend's faint end is marginal, but
the monotonic ordering across the three bins is clean. This is exactly the step a single stacked number
cannot provide.

## Redshift-binned: a non-monotonic observed-flux run

Re-binning the **same** 279 cutouts by SDSS redshift (each bin separately injection-recovered) costs no
extra data and gives the run of *observed-frame* mean flux with redshift:

| $z$ bin | median $z$ | mean radio flux | SNR |
|---|---|---|---|
| low | 0.91 | 57.7 µJy/beam | 3.7 |
| mid | 1.92 | **64.0 µJy/beam** | 4.4 |
| high | 2.58 | 41.7 µJy/beam | 2.8 |

Unlike the magnitude trend, this is **not monotonic** — the observed mean flux peaks in the middle
($z\simeq1.9$) and the highest-$z$ bin is the faintest *and* only 2.8$\sigma$. We deliberately **do not
over-read** this: an observed-frame flux folds the luminosity distance, the $K$-correction, and any
intrinsic luminosity evolution together, and the bins are wide. It is reported as a **demonstration
that the binning machinery generalises to any catalogue property** — not as an evolutionary
measurement. (The magnitude trend, where brighter optical → brighter radio is unambiguous, is the clean
result; the redshift run is the honest "it also runs on $z$, and here's what the data say" companion.)

## Honest assessment & caveats

- **Methodology + a calibrated measurement, not a discovery.** The contribution is a tested
  stacking-plus-injection-recovery pipeline demonstrated end-to-end on public data, and one calibrated
  population-mean flux.
- **A 4.5σ stack is a marginal-to-solid detection, not a high-significance one** — a larger sample
  (more sources, or co-adding fields) would tighten it; this run is a demonstration, not a survey.
- **SE coverage drops sources.** 6% of targets had no SE image and were dropped; the SE footprint is
  still filling in, so a different field would give a different N.
- **Clean-PSF injection-recovery** calibrates the flux-scale/snapshot bias (and finds it negligible
  for SE here) but does not model every deconvolution subtlety; confusion and bright-neighbour
  contamination are mitigated by the median but not eliminated.
- **Binning is a demonstration, not a measurement.** The magnitude- and redshift-binned trends show
  the machinery generalises and recover a clean radio--optical correlation; the redshift run is
  non-monotonic and observed-frame, so it is *not* read as luminosity evolution. A larger sample would
  tighten the marginal (faintest / highest-$z$) bins.

## Referee round on the style conversion (2026-08-24)

All 11 macros re-verified against `results/stacking_metrics.json`; no number moved. Fixed:

1. **MAJOR.** "The honest limitations **are** those of any catalogue-and-cutout stacking" had become
   "**Several** limitations **apply to any** ... stacking". Two things changed beyond tone: the
   grammatical subject moved from *this work's* limitations to the method class (while the list that
   follows still opens "the result depends on the Single-Epoch coverage at the target positions",
   i.e. this paper's result), and a definite claim of completeness became an explicit subset claim.
   That downgrade is not survivable here: the dropped-source selection has **no stated denominator**
   anywhere (`\stN` = 236 is the *stacked* count, and the number of DR16Q targets queried is in
   neither the paper nor the metrics), the injection-recovery ratio of exactly 1.0 is the sole
   support for the paper's central word *calibrated*, and the binned trends the abstract headlines
   rest on bins of SNR 2.6--4.7. Restored to a definite, this-work-owned sentence.
2. The abstract's global "not a new astrophysical claim" had become the last item of a comma series,
   so it read as qualifying only the third contribution rather than all of them. Given its own clause.
3. "**We do not over-read this:**" was deleted from the redshift paragraph. The operative limitation
   survived, but that clause was the only signal that what follows is a *refusal to interpret*.
   Restored.
4. A sentence split let a bin-level ratio generalise: "The optically-brightest quasars are ~2x
   radio-brighter on average" acquired a generic plural subject, where under the em-dash it was a
   gloss tethered to the two numbers just quoted (82.8/43.8 = 1.89 in these three bins of these 236
   objects). Re-tethered.
5. Deleting "This is the kind of population trend a single stacked number cannot show" was
   **correct** --- the point is made twice already and the deleted form was an unsupported
   impossibility claim. But it left the paragraph *ending* on "the monotonic ordering across the
   three bins is clean", and terminal position is emphatic. Checked against the committed bins
   (error = flux/snr, since `snr = peak/rms`): adjacent steps are **1.00 sigma** and **0.60 sigma**,
   end-to-end **1.60 sigma**. The ordering is monotonic in the central values only; no adjacent pair
   is resolved. Reworded to say so.

**Also fixed (pre-existing):** the magnitude paragraph pointed the reader at
"(Figure~\ref{fig:stack}, **right**)" where the caption puts the magnitude panel in the **middle**,
and the redshift paragraph points at "right" too. Both pointed at the same panel.

## Full referee round (2026-08-25): MAJOR REVISION, 15 findings, two BLOCKERs

The measurement is probably real (the annulus RMS scales as N^-1/2 across a factor of 3 in N to
1.4% — a genuine internal check that passes), but the title claim and every quoted number rest
on two defects the referee proved rather than suspected.

**BLOCKER 1: the injection-recovery ratio is an algebraic identity, not a measurement.**
`injection_recovery` adds the same PSF plane to every cutout and differences two sigma-clipped
medians: sigma_clip's cenfunc/stdfunc are shift-equivariant, so the clip mask is unchanged and
the recovered value ≡ the injected value — ratio ≡ 1.0 for ANY input (verified analytically and
numerically on a hostile cube: heteroscedastic gains, interlopers, NaN edges — 1.0 exactly at
every amplitude and FWHM). The title's "Injection-Recovery Calibration", the abstract's
"de-biases the result", "the SE flux scale is unbiased (unlike Quick-Look would be)", and "the
same stack without the step would quote a biased flux" are all unsupported — the identical code
returns 1.0 on Quick-Look or pure noise; \stPeak and \stDebiased are the same number; and two
unit tests assert conditions the identity makes unfailable (the "test can lock a defect in"
pattern).

**BLOCKER 2: every flux and SNR is a 3-pixel SEARCHED peak** — the exact stokesv lesson written
in vlass.py's own line-1022 docstring. Monte Carlo of the code's own geometry: the searched
peak on pure beam-correlated noise reads +1.57×RMS (positive 99% of the time), so the headline
4.5σ is ~3.6σ honest, 43.5 µJy is ~28–32, and the faintest bin's 2.6σ is ~1.4σ. The committed
JSON carries the fingerprint twice: the count-weighted bin means (61.6/59.9 µJy) exceed the
full-sample 43.5 by the 1.42× a pedestal predicts (bins have 1.75× the RMS). And the
calibration leg measures the CENTRE pixel while the science leg searches — the two legs don't
even measure the same statistic.

**MAJORs:** no off-source control stack exists anywhere (an annulus is blind to a
centre-common pedestal — the control leg would turn the referee's simulation into the paper's
own measured null); "quasars that lack an individual VLASS detection" is not implemented (no
flux cut, no catalogue cross-match — the sample is every DR16Q row with a cutout); no
denominator/footprint/target list committed (n_queried, radius — code default 3.0° vs the
findings doc's "2.5° cone" — max_sources=300 row-limit truncation: the effective sky area is
undefined and "regenerates from a clean checkout" is not true); the magnitude bins are
APPARENT magnitude over z 0.9–2.6 (the "expected radio–optical luminosity correlation ...
recovered" is inseparable from a luminosity–distance effect, end-to-end 1.6σ, and the "~2×"
becomes 3.2× under the pedestal correction); the findings doc's main narrative documents a
DIFFERENT earlier run (N=279/SNR 4.9/45.4 µJy vs committed 236/4.5/43.5) including the only
stated denominator ("of 300 tried, 279"); the paper says "mean" throughout while the estimator
is a 3σ-clipped MEDIAN (for a skewed radio-loudness distribution these differ by a large
factor — illustration: true mean 158, clipped-median stack 38 — and the offline fixture gives
every source identical flux, so mean ≡ median by construction and no test can see it).

**MINOR/NIT:** White/Karim/Lindroos are cited as plain text and missing from the reference
list (main.bbl has 7 entries); √N noise sentence wrong for a median (1.25σ/√N); 4.44-vs-4.5
rounding inconsistency; mixed mJy/µJy units hide the bins-above-full-sample smoking gun;
\stSource generated but unused (field coordinates never reach the reader); macros
mode-dependent and un-namespaced (guards wired, defence-in-depth holds); the local
arxiv-submission tarball predates the honesty edits.

**Status: fixes pending.** The single change: replace the searched peak with a forced
central-pixel measurement, re-derive all seven fluxes/SNRs, and add the off-source control
stack — then either retract the calibration claim or build an injection test that can fail.

**Status: RESOLVED (2026-08-25).** One real re-run (`--ra 180 --dec 25 --radius 3`), fetching
science AND 30″-offset control cutouts (598 stamps total). SE coverage has grown since the
committed run: 300/300 queried targets now have an SE image (was 236), so every number moved
with the sample as well as the estimator; the committed denominator chain is
300 queried → 300 with cutout → 19 individually detected (5σ per-cutout test, now implemented,
not asserted) → **281 stacked**, all in `results/stacking_targets.csv`.

1. **BLOCKER 2 (searched peak) fixed**: `measure_stacked_flux` reads the FORCED central pixel;
   the searched max survives only as a labelled diagnostic, and a committed test asserts forced
   photometry goes negative ~half the time on pure noise while the searched peak is positive
   >90% of the time. Re-derived headline: **40.8 µJy/beam at 4.9σ** (annulus RMS 8.3). The
   searched max equals the forced value on the science stack — a real centred source. The
   referee's predicted "28–32 µJy at ~3.6σ" was computed for the old N=236 stack; on the
   current 281-source stack the forced value stands at 40.8.
2. **BLOCKER 1 (identity calibration) fixed by retraction + replacement**: the paper now states
   the shift-equivariance identity outright and retracts the "injection-recovery calibration"
   framing (title changed). The replacement injects per-cutout sub-pixel offsets (σ = 0.3 px,
   a documented assumption) at the MEASURED amplitude over 8 draws: ratio **0.923 ± 0.024** —
   a number that can differ from 1 and does. A test proves the new version fails (ratio < 0.8)
   for a badly-centred population. The ratio is reported as a conditional systematic, not
   folded into the headline.
3. **The off-source control exists**: 298 control stamps through the identical pipeline give a
   forced flux of **−17.5 ± 8.2 µJy/beam (−2.1σ)** — no positive centre-common pedestal; the
   paper notes the mildly negative level is consistent with zero at ~2σ and, if real, samples
   the CLEAN-residual environment 30″ from real sources (conservative direction).
4. **The sample cut is implemented** (`individually_detected`): 19 sources flagged ≥5σ in their
   own cutout are excluded and counted.
5. **The magnitude trend is restated honestly**: bright third 72.0 ± 14.0 (5.1σ), middle 41.5
   (2.8σ), faintest **9.1 ± 14.5 (0.6σ — a non-detection)**. No bright/faint ratio is quoted
   (the denominator is consistent with zero), and the run is framed as observed-frame machinery
   (apparent magnitude over z ≈ 0.9–2.6), not the radio–optical luminosity correlation. Same
   for redshift: 47.3 / 45.5 / 12.2 (0.9σ) µJy/beam.
6. **"Mean" → clipped median throughout**, with a log-normal fixture (`flux_scatter_dex`) whose
   committed test shows the stack recovers the median while the mean is ~2.6× larger; the noise
   sentence uses 1.25σ/√N.
7. White/Karim/Lindroos are real `\citep` entries (Crossref-verified; main.bbl resolves them);
   `\stSource` and the field/radius/max_sources reach the reader as macros; units are µJy/beam
   throughout; synthetic runs can no longer clobber the real figure (guarded + tested); the
   stale earlier-run narrative above is banner-marked SUPERSEDED; arXiv package regenerated
   (0 errors).
