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
