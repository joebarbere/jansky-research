# Findings — solar type III drift → exciter speed (e-Callisto)

`jansky_research.solarbursts` fits a type III radio burst's frequency drift in an e-Callisto dynamic
spectrum and inverts it, via the Newkirk coronal density model, to the exciter (electron-beam) speed,
reusing `jansky.solar.density_from_plasma_frequency` + `newkirk_radius`. This is the real-data +
recover-a-known leg.

## Event selection (the gating step)

The **Monstein SGD burst lists** (`solarradio/data/BurstLists/2010-yyyy_Monstein/`) catalogue every
e-Callisto burst with type, quality, and frequency range. A flare interval is a type III **storm** —
many overlapping bursts with no single drift (our first attempt, the 2011-08-09 X6.9 flare, smears into
an artificially slow, incoherent ridge). The fix is to pick **isolated, quality-1** type III bursts; the
tool's fit **R²** then empirically flags a coherent single drift (≈0.9) versus a storm (<0.3) — a useful
indicator, not a rigorous segmentation.

Four isolated quality-1 events at Birr (BIR, 10–90 MHz) were analysed. **Only the 2011-09-14 event is a
clean fit (R² = 0.90); it anchors the result.** The other three formally land in the canonical 0.1–0.5 c
band but their fits are marginal (R² < 0.55), so their individual speeds are not reliable:

| date (UT) | drift (MHz/s) | R² | speed (harmonic, 1× Newkirk) | fit |
|---|---|---|---|---|
| **2011-09-14 11:50** | **−3.3** | **0.90** | **0.14 c** | **clean** |
| 2011-07-05 10:54 | −9.8 | 0.48 | 0.27 c | marginal |
| 2011-07-16 13:10 | −17.9 | 0.52 | 0.45 c | marginal |
| 2011-09-19 07:43 | −2.6 | 0.41 | 0.11 c | marginal |

The 2011-09-14 drift (−3.3 MHz/s) is *slower* than typical: the Alvarez & Haddock (1973) and Zhang et al.
(2018, A&A 618, A165, median −6.94 MHz/s) relations give ~−6 to −9 MHz/s at mid-band, so this is a
relatively slow event (≈25th percentile) — a fair test of the pipeline on a weak drift, not a typical burst.

## Recover-a-known: the 2011-09-14 burst (the cleanest, R² = 0.90)

A single isolated type III at BIR, ridge spanning **10–79 MHz**, 66 channels (61 after robust
sigma-clipping), drift **−3.3 MHz/s**, coherent height–time track (**R² = 0.90**). Mapping the ridge
through Newkirk gives heliocentric heights **1.8–2.4 R⊙** and an exciter speed:

| emission mode | density model | speed | heights |
|---|---|---|---|
| fundamental | 1× Newkirk | **0.086 c** | 1.41–1.80 R⊙ |
| harmonic | 1× Newkirk | **0.137 c** | 1.76–2.41 R⊙ |
| harmonic | 4× Newkirk | **0.272 c** | 2.33–3.63 R⊙ |

The recovered speed (**~0.09–0.27 c** across the assumptions) sits in the established type III range
(0.1–0.5 c; Reid & Ratcliffe 2014). The headline 0.137 c (harmonic, 1× Newkirk) is a **peak-time** speed
and is directly consistent with the **peak-time** mean of 0.17 c from LOFAR type III imaging (Reid &
Kontar 2018) — no extra de-biasing assumed. The synthetic fixture confirms the inversion is
**algebraically self-consistent** (injected 0.2/0.3/0.4 c recovered within 10%) — a round-trip code
check, since the forward fixture and the inverse use the same Newkirk mapping; it does *not*
independently validate the density model against the real corona.

(Heights and `r_lo`/`r_hi` are from the sigma-clipped ridge, not the raw band edges.)

## Honest assessment & caveats

- **The two model knobs span a factor ~3.** Emission mode (fundamental vs harmonic) and the Newkirk
  fold (1× quiet → 4× active region) move the speed from 0.09 to 0.27 c — so we quote the **grid**, not
  one number. Ground-based 20–80 MHz bursts are often harmonic; 1× Newkirk is the quiet-corona default.
  The Newkirk model is best constrained within ~1–3 R⊙; the 4× / 3.6 R⊙ case extrapolates near its edge.
- **Peak time, not onset.** The ridge is built from each channel's *peak* intensity, which traces the
  bulk of the electron beam, not its leading edge; peak-time speeds run ~15–30% below front-of-beam
  speeds (Reid & Kontar 2018). Radio-wave scattering and projection further suppress apparent speeds at
  <50 MHz (Kontar et al. 2017) — both biases act *downward*, so the true front speed is a little higher
  than our peak-time number, not lower.
- **R² is an empirical flag, not a rigorous test.** On the 2011-08-09 flare *storm* the pipeline gives
  R² < 0.3 and an incoherent drift, so the metric usefully flags "no single burst" — but a proper
  discriminator would segment bursts on the dynamic spectrum, and R² is sensitive to the point spacing.
- **Uncalibrated, single-station, ionosphere-limited.** e-Callisto data are arbitrary digitiser units
  (no flux calibration) — this is a purely morphological drift analysis. Only Birr (BIR) was used, so
  there is no second-station confirmation. Ridge channels at ≲12 MHz may be affected by ionospheric
  refraction near the ground-based cutoff and add uncertainty to the high-altitude end of the fit.
- **A recover-a-known, not a survey or a discovery.** Four hand-picked clean bursts validate the tool
  end-to-end on public data; a blind drift-rate catalogue would need automated burst segmentation and
  RFI flagging across many stations.
- **Reproducible:** `python -m jansky_research.solarbursts --recover` regenerates the 2011-09-14 result,
  the metrics JSON, the dynamic-spectrum + height–time figure, and the macros from the public FITS.

## Stale systematics grid (found 2026-08-23, partially fixed)

The Results paragraph quotes a harmonic/fold grid --- "$0.086\,c$ (fundamental, $1\times$) through
$0.137\,c$ (harmonic, $1\times$) to $0.272\,c$ (harmonic, $4\times$)" --- whose middle value
contradicted `\sbSpeedC` = **0.1347** three lines above it, for the *same* grid point. The grid
comes from a superseded run that differs from the committed evidence on every axis:

| quantity | this findings file (earlier run) | `results/solarbursts_metrics.json` |
|---|---|---|
| R^2 | 0.90 | **0.811** |
| drift | -3.3 MHz/s | **-2.55** |
| heights, harmonic 1x | 1.76--2.41 R_sun | **1.622--2.572** |
| ridge channels kept | 61 | **62** |
| speed, harmonic 1x | 0.137 c | **0.1347** |

**Fixed:** the harmonic-1x point now cites `\sbSpeedC`, so that value is regenerable and cannot
drift again.

**NOT fixed, deliberately:** the flanking values (0.086, 0.272) and the "3.6 R_sun" figure in the
Discussion are still hand-typed from the superseded run. They cannot be recomputed from committed
evidence, because `exciter_speed` needs the raw ridge (frequencies and times) and only its summary
is in the metrics. The conclusion is unaffected --- the bracket is still ~0.09--0.27 c, inside the
canonical band --- but the grid is not auditable. The fix is to emit the full harmonic x fold grid
into `results/solarbursts_metrics.json` on the next real run and macro-ise all three points.

### Resolved 2026-08-24 (real re-run): the grid is emitted, and the discrepancy is explained

`speed_grid` now computes all three systematics points from the same fitted ridge as the headline
number, so the middle point IS `\sbSpeedC` by construction; the raw ridge is committed as
`results/solarbursts_ridge.csv` (66 channels), so both the headline and the grid are recomputable
from evidence. The paper's grid sentence, abstract bracket and Discussion radius now cite macros.

**The stale-grid mystery is solved, and it was a parameterization, not code drift.** The committed
metrics (r2 0.811, drift -2.55, 0.1347) are reproduced byte-for-byte by `pad_s=10.0` (the `run()`
default); the findings file's "superseded run" (r2 0.90, drift -3.3, 0.137) is exactly what
`pad_s=5.0` gives, and `--recover` had 5.0 pinned in it --- so the one command this file said
"regenerates the 2011-09-14 result" in fact produced different numbers. `--recover` now pins 10.0,
matching the committed evidence, and a re-run against the live archive confirms the pipeline is
deterministic and the archive unchanged.

The committed grid (pad 10): **0.0813 / 0.1347 / 0.2448 c**, a factor of 3.01. Two prose
corrections fell out: "All three lie within the canonical 0.1--0.5 c band" was false (the
fundamental point is below 0.1, and was at 0.086 in the old prose too --- pre-existing); and the
Discussion's "3.6 R_sun" outer radius was the pad-5 value, where the committed grid reaches
**4.01 R_sun**, which extrapolates *beyond* the Newkirk model's well-constrained range, not "near
its edge". Both now stated honestly.

Pad sensitivity, for the record: 5 -> 10 moves r2 0.897 -> 0.811, drift -3.259 -> -2.55, speed
0.1368 -> 0.1347. Both are defensible windows; the committed choice is the wider one the paper's
numbers have always cited.

## Full referee round (2026-08-24): MAJOR REVISION, 18 findings, one BLOCKER

The referee reproduced the committed headline exactly from the committed ridge CSV (0.1347 /
62 / 0.811 / 1.622-2.572 and all three grid points) -- committing the ridge was what made the
audit below possible.

**BLOCKER: Figure 1's plotted line is not the fit the caption describes.** The right panel fits
all 66 ridge points unweighted (slope -> 0.1154 c, R^2 = 0.411 on the points shown, and the
clipped 4.97 R_sun outlier is drawn unmarked) under a caption asserting R^2 = 0.811. A reader
measuring the slope off the axes gets a speed 14% below the headline.

**MAJORs:** the abstract binds the 10.0-78.94 MHz span and the -2.55 MHz/s drift to the
62-channel fit, but the kept channels span 25.44-78.94 MHz and the drift is the unclipped
66-point polyfit (robust drift: -3.27, a 28% difference -- and the drift feeds the "25th
percentile" caveat); `_robust_linfit` is NOT converged at its hard-coded n_iter=3 -- the speed
ranges 0.111-0.147 c and R^2 0.68-0.94 purely over the iteration count, the returned slope and
mask come from different iterations (refitting on the reported 62 gives 0.1469 c), and
convergence is at 55 points / 0.1173 c (the innerrc parameter-on-a-bound shape, in a loop
bound); no uncertainty anywhere and 40377.5 km/s is six significant figures on a quantity with
a +/-13% undeclared analysis systematic; the event-selection function (the paper's stated crux)
has no committed evidence and the three rejections were adjudicated at pad_s=5 (superseded) at
an unconverged n_iter -- commit a candidates file incl. the storm control at the headline
parameterisation; the recover-a-known cannot fail (target band spans 5x, model grid spans 3x,
and all three REJECTED events also land in-band 0.11-0.45 c) -- the storm control is the test
that can fail and it is the one not committed; harmonic=2 is the assumption that puts the
headline in-band and is justified nowhere in the manuscript (R^2 does not prefer it: 0.852
fundamental vs 0.811 harmonic); "directly consistent with the 0.17 c peak-time mean" compares
one event to a 31-burst ensemble mean derived through a different density model, with the
paper's own grid showing that choice moves answers 3x.

**MINOR:** the results JSON omits pad_s/snr_threshold/n_iter/clip sigma (pad_s is the exact
parameter whose invisibility caused the stale-grid incident) and the ridge CSV has no
provenance header; both quoted band extremes rest on single isolated channels, one of which
(78.94 MHz at t=367.0) drifts the WRONG WAY for a type III and sets both f_hi and r_lo; the
abstract says the fundamental point is "at ... the lower half" of the canonical band when it is
19% below its lower edge (the previously-retracted formulation's shape); "25th percentile" is
unsourced and computed from the contaminated drift; "15-30% below front-of-beam" -- the cited
numbers give 15% (peak) and 25% (back), not 30%; the findings file's first two-thirds present
superseded pad-5 numbers as current; the untracked arxiv-submission still carries the retracted
"squarely in the established range" abstract.

**NIT:** RECOVER_EVENT docstring says R^2 ~ 0.9 and still carries pad_s 5.0 while --recover
hard-codes 10.0; \sbTruth/\sbRatio render as -- (namespace as \sbSyn*); left panel shows
spectrum to ~110 MHz vs "10-90 MHz" prose; kontar2017 pages field carries an article number.

**The one change that matters most: iterate the robust fit to convergence and report the
headline as an interval over the analysis choices** -- that converts the weakest thing in the
paper into its actual contribution.

**Status: fixes pending** (all fixable from the committed ridge + one candidates re-run into a
scratch dir).
