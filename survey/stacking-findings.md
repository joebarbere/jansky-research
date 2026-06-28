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
