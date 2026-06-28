# Findings — southern peaked-spectrum catalogue (GLEAM-X + RACS), real data

`jansky_research.southern` selects peaked-spectrum (GPS/CSS-candidate) and ultra-steep-spectrum (USS)
sources in the southern sky by fitting a log-parabola to a multi-band SED built from **GLEAM-X DR2**
(20 in-band sub-bands, 76–227 MHz; VizieR `VIII/113`) and the three **RACS** bands (887.5, 1367.5,
1655.5 MHz; `J/other/PASA/38.58`, `41.3`, `42.38`). Unlike the northern `peaked` slice — which had to
use TGSS 150 MHz as an *upper limit* — the GLEAM-X in-band points let `fit_log_parabola` **measure**
the turnover frequency $\nu_\mathrm{pk}$.

## The unit trap (caught up front)

GLEAM-X DR2 integrated fluxes (`Fint###`) are in **Jy**; RACS fluxes (`Ftot`) are in **mJy**. Mixing
them would put GLEAM-X $10^3\times$ too low and fake a steeply rising spectrum into RACS for *every*
source. `fetch_gleamx` converts Jy→mJy (×1000). (Same class of bug as the VLASS µJy/mJy trap.)

## Real run: a 3° cone at RA 30°, Dec −30° (high galactic latitude, GLEAM-X footprint)

**1545 GLEAM-X∩RACS matched sources.** A naive log-parabola-peak selection over the SED gave **246
"peaked" (16%)** with a median $\nu_\mathrm{pk}$ piling at the GLEAM band edge — not a result, two
selection effects, each with a principled fix (mirroring the northern slice):

1. **GLEAM (2′) ↔ RACS (15–25″) resolution mismatch.** Sources extended at GLEAM's 2′ beam lose flux
   at RACS resolution, so the spectrum drops GLEAM→RACS and fakes a high-frequency turnover. Fix: a
   **compactness cut** using the GLEAM-X integrated/peak ratio (`Fintwide/Fpwide` $<1.2$); non-compact
   sources are flagged `extended`, not peaked. (246 → 197.)
2. **Steep-source low-frequency flattening.** An ordinary steep source with mild low-frequency
   flattening fits a shallow parabola peaking spuriously at the GLEAM band edge. Fix: require the
   **optically-thick rising side** — the measured GLEAM in-band index $\alpha_\mathrm{lo}>-0.1$ (the
   GPS/SSA signature), so a real peaked source rises below its turnover. (197 → 90.)

**After both cuts:** **90 peaked candidates (5.8%)**, **59 USS** (candidate high-z radio galaxies,
$\alpha<-1.2$ throughout), **28 extended** rejected, **median measured $\nu_\mathrm{pk}\approx211$
MHz** — i.e. turnovers in the few-hundred-MHz range, the intermediate population this method targets.

## Honest assessment

A reproducible, maximal-reuse multi-band **measured-turnover** selector — the upgrade the northern
slice's upper-limit method could not provide. It is a **scout catalogue + methodology** result, not a
discovery: 90 peaked candidates over 28 deg² is ~2.5× the surface density of the careful Bayesian
RadioSED II (Kerrison+2025, ~1.2 deg$^{-2}$ over Stripe 82), so some over-selection remains — expected,
since a least-squares parabola is more permissive than a Bayesian SSA/FFA model and GLEAM-X DR2 exposes
**no per-band errors** via VizieR (a nominal 10% per-band error is used for the fit weights).

## Recover-a-known: measuring the Callingham (2017) turnovers

`validate_callingham` is the validation the northern slice could not do — it *measures* the published
turnover rather than bounding it. For each Callingham peaked source (VizieR `J/ApJ/836/174/pkfreq`,
which gives a measured $\nu_\mathrm{pk}$) it fetches a small GLEAM-X+RACS cone, runs the pipeline, and
compares. On **40** sources (all with coverage):

| quantity | value |
|---|---|
| recovered as peaked | **30/40** |
| median $\lvert\Delta\log_{10}\nu_\mathrm{pk}\rvert$ (measured vs published) | **0.12 dex** (~30%) |
| within a factor of two (0.3 dex) | **90%** |
| recovered by published $\nu_\mathrm{pk}$: 72–250 / 250–500 / 500–2000 MHz | 23 / 6 / 1 |

So the method recovers 75% of Callingham peaked sources and reproduces their literature turnover to a
median ~0.12 dex — the measured-turnover headline the southern data enable. Recovery climbs with
$\nu_\mathrm{pk}$ (sources peaking well below the GLEAM band fall outside the rising-side window, as
designed).

### Limitations / next steps
- **No per-band GLEAM-X errors** → unweighted-ish fit; `chi2_red` is only nominal. A proper Bayesian
  SED (per-band rms from the images) would tighten purity.
- **Flux-scale tie** GLEAM-X↔RACS (~10% each) and the GLEAM↔RACS frequency gap (227 MHz → 887 MHz)
  mean $\nu_\mathrm{pk}$ in that gap is interpolated; a sampled point in the gap (e.g. SUMSS 843 MHz,
  or RACS-low at 888 MHz which we already use) constrains it.
- **Compactness cut is blunt** (GLEAM 2′ only); a RACS-resolution morphology cut would be cleaner.
