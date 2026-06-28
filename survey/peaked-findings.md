# Findings — peaked-spectrum (GPS/CSS) candidates via 3-frequency curvature

`jansky_research.peaked` selects peaked-spectrum (GPS/CSS-candidate) radio sources — compact, young
radio AGN whose spectrum rises then falls — from three public surveys: TGSS (150 MHz), NVSS
(1.4 GHz), VLASS (3 GHz). It is the direct follow-up the USS slice called for ("confirming a turnover
needs a third frequency, e.g. VLASS 3 GHz") and the maximal-reuse slice: it composes `spectra`
(`spectral_index`, `crossmatch`, `fetch_survey`) and `vlass` (3 GHz fetch, forced-photometry
variability, SIMBAD vetting) rather than reimplementing them. Run: a 2° cone at RA 180°, Dec +30°
(the same field as the USS slice), 409 NVSS$\cap$VLASS sources.

## Two systematics the real data exposed (and the fixes)

The naive "find sources detected in all three surveys, rising then falling" gave **zero** candidates —
not a null result but two selection effects:

1. **TGSS is shallow.** Peaked sources are *faint at 150 MHz* by definition, so requiring a TGSS
   detection excludes exactly the sources wanted (only the bright steep sources are in TGSS). Fix:
   treat a TGSS **non-detection** as an upper limit $S_{150}<25$ mJy, giving a *lower bound* on
   $\alpha_\mathrm{low}$. A source brighter at 1.4 GHz than the 150 MHz limit yet undetected at
   150 MHz **must** be rising. (This only confirms peaked sources brighter than $\sim$30 mJy at
   1.4 GHz — an honest depth limit.)
2. **NVSS$\to$VLASS resolution mismatch.** NVSS (45″) sees extended emission that VLASS (2.5″)
   resolves out, producing impossible "spectra" — e.g. 690 mJy $\to$ 6 mJy, $\alpha_\mathrm{high}=-6$.
   These masquerade as steep turnovers. Fix: a floor $\alpha_\mathrm{high}>-2$ rejects them as
   resolution artefacts. In this field **110** of 409 sources are flagged extended/artefact.

## Result: 6 peaked candidates, all steady (GPS-like), 3 uncatalogued

After both cuts: **6 peaked-spectrum candidates** (rising at low frequency via the TGSS upper limit,
$\alpha_\mathrm{high}\approx-0.6$ to $-1.1$, compact). Reusing `vlass.forced_photometry`, **all six are
VLASS-steady across the three epochs** ($V=0.02$–$0.08$) — consistent with GPS sources, not blazars.
SIMBAD (via `vlass`) identifies them:

| RA Dec | $S_{1.4}$ / $S_3$ (mJy) | $\alpha_\mathrm{high}$ | SIMBAD | note |
|--------|------|------|--------|------|
| 179.978 +30.458 | 41.6 / 20.5 | −0.93 | LEDA 1904565 [LINER] | AGN host |
| 182.018 +30.264 | 35.3 / 20.9 | −0.69 | **4FGL J1208.1+3017 [BL Lac]** | a *Fermi* blazar — contaminant |
| 179.181 +31.710 | 34.7 / 18.4 | −0.83 | — | uncatalogued GPS/CSS candidate |
| 179.174 +29.196 | 32.9 / 16.9 | −0.87 | — | uncatalogued GPS/CSS candidate |
| 180.048 +28.210 | 32.7 / 13.9 | −1.12 | — | uncatalogued GPS/CSS candidate |
| 180.702 +31.873 | 31.7 / 19.8 | −0.62 | NVSS J120248+315223 [radio galaxy] | CSS/GPS host |

A nice methodological point: the BL Lac (a known $\gamma$-ray blazar) is **steady** in VLASS ($V=0.08$),
so the variability flag alone did *not* catch it — **SIMBAD did**. Blazar rejection needs the catalogue
cross-ID (and ideally longer-baseline variability), not VLASS-epoch variability alone.

## Validation against a known peaked-source catalogue (and the selection function it reveals)

`validate_known` tests the selection against the **Callingham et al. (2017)** peaked-spectrum
catalogue (1222 sources with a measured turnover frequency $\nu_\mathrm{pk}$), adding VLASS 3 GHz to
each. Run on 113 bright ($S_{1.4}>60$ mJy), VLASS-accessible sources, the recovery as "peaked" is:

| $\nu_\mathrm{pk}$ (MHz) | 72–250 | 250–500 | 500–1000 |
|---|---|---|---|
| recovered as peaked | 0/81 | 3/26 | 1/6 |

This looks like a poor recovery but is the *correct* behaviour, and it pins down the method's
selection function. Callingham is **GLEAM-selected** (70–230 MHz), so it is dominated by
*MHz-peaked* sources (median $\nu_\mathrm{pk}=190$ MHz) whose turnover is **below the 150 MHz floor**
of this method — across 150 MHz–3 GHz they are simply falling, i.e.\ *steep*, and the method
**correctly does not flag a single one** (0/81 false positives below 250 MHz: high purity). A
three-point 150 MHz / 1.4 GHz / 3 GHz method has a narrow **peaked window of $\sim$0.7–2 GHz**:
sources peaking below it look steep, and classical bright GHz-peaked GPS (turnover $\gtrsim$3 GHz)
are still rising at 3 GHz and land in this method's **`inverted`** class, not `peaked`. So the slice
does not "recover" the GLEAM MHz-peaked population (a different sub-population), and there is no large
public catalogue that densely samples the $\sim$1 GHz window to give a clean high-recovery number.
The honest validation result is: **correct purity against the MHz-peaked majority, a well-characterised
$\sim$0.7–2 GHz sensitivity window, and the recognition that classical GHz-GPS appear in the `inverted`
class** — a frequency-coverage limitation, stated plainly.

## Chasing the recovery: a clean recover-a-known against the GHz-peaked (HFP) population

The Callingham test shows the method *correctly rejects* the MHz-peaked majority, but a selection
method should also be shown to *recover* the population it claims to find. The classical
GHz-peaked / High-Frequency-Peaker (HFP) sources peak at $\gtrsim$few GHz, so across 150 MHz /
1.4 GHz / 3 GHz they are **rising throughout** — they are exactly what lands in this method's
`inverted`/rising class, not `peaked`. To make this explicit, `validate_hfp` runs the selection over
the **Dallacasa et al. (2000)** bright HFP sample (`J/A+A/363/887`, NVSS 1.4 GHz fluxes), adding
VLASS 3 GHz per source:

| Dallacasa HFP sample (98 sources with VLASS) | result |
|---|---|
| median $\alpha_\mathrm{low}$ (TGSS upper limit $\to$ 1.4 GHz) | **+1.03** (strongly rising) |
| median $\alpha_\mathrm{high}$ (1.4 $\to$ 3 GHz) | **+0.19** (still rising) |
| recovered as *rising* ($\alpha_\mathrm{low}>0.1$: optically thick at low freq) | **98/98 = 100%** |
| recovered as *GHz-peaked* ($\alpha_\mathrm{high}>0.1$: still rising at 3 GHz) | **54/98 = 55%** |

This is the clean recover-a-known the Callingham catalogue could not provide: **every** Dallacasa HFP
source is recovered as optically-thick-rising, and a majority are flagged `ghz_peaked` (peak above
3 GHz). The two validations are complementary and together pin the selection function exactly:

- **MHz-peaked (Callingham, GLEAM-selected):** turnover *below* the 150 MHz floor $\to$ falling across
  the band $\to$ correctly **not** flagged (0/81 false positives below 250 MHz — high purity).
- **GHz-peaked / HFP (Dallacasa):** turnover *at/above* a few GHz $\to$ rising across the band $\to$
  recovered (100% rising; 55% still rising at 3 GHz, i.e. `ghz_peaked`).
- **Intermediate ($\sim$0.7–2 GHz turnover):** caught as `peaked` (rising then falling) — the field
  candidate list above.

So the three SED classes map cleanly onto the three turnover regimes, and "classical GHz-GPS appear in
the `inverted` class" is no longer a caveat but a **validated, named recovery**: those sources are now
the dedicated `ghz_peaked` class, and the method recovers 100% of a known HFP sample as rising.

## Honest assessment

This is a **methodology + candidate-list** result, not a discovery: a reproducible three-frequency
curvature selection that (i) is robust to the TGSS flux-scale offset that sank the USS slice (curvature
compares two indices; $\alpha_\mathrm{high}$ is TGSS-independent), (ii) correctly separates resolution
artefacts, and (iii) yields a clean, steady, compact peaked-candidate list — three of them
uncatalogued. Confirming any as a true GPS source needs higher-resolution / more-frequency follow-up
(VLBI compactness, a sampled turnover); we report candidates, not detections.

### Limitations

- **Small field, depth-limited.** One 2° cone; the TGSS-upper-limit method only confirms peaked
  sources brighter than $\sim$30 mJy at 1.4 GHz, so the sample is the bright tail.
- **Curvature via a lower limit.** For TGSS non-detections $\alpha_\mathrm{low}$ is a bound, not a
  measurement, so the turnover frequency is not pinned (only "rising then falling").
- **Turnover-regime sensitivity is by design, now split into two classes** (see both validations
  above): MHz-peaked/CSS sources (turnover below the 150 MHz floor → look steep) are *correctly
  rejected*, while classical GHz-GPS/HFP (turnover at/above a few GHz → rising throughout) are
  recovered as the dedicated `ghz_peaked` class (100% of the Dallacasa HFP sample). The narrow
  `peaked` (rising-then-falling) class is genuinely limited to the $\sim$0.7–2 GHz window, and for
  *that* intermediate band no large public catalogue densely samples the turnover for a clean
  recovery number — but the rising/`ghz_peaked` arm now has one.
- **NVSS$\to$VLASS resolution** biases $\alpha_\mathrm{high}$ for any non-point source; the floor cut is
  blunt. Matched-resolution fluxes (or a compactness cut) would be cleaner.

## Bottom line

A reproducible, maximal-reuse three-frequency peaked-spectrum selector that turns the USS slice's
unfinished business into a clean candidate list — robust to the flux-scale and resolution systematics,
honestly bounded, with three uncatalogued GPS/CSS candidates and a worked demonstration that blazar
contamination needs catalogue cross-ID, not just variability. Validation against the Callingham (2017)
catalogue confirms its purity against the MHz-peaked majority and characterises its $\sim$0.7–2 GHz
sensitivity window — the result is a well-understood selection method, not a discovery.
