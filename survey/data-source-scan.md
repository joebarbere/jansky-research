# Survey — under-mined radio data sources, non-traditional datasets, and GPU leverage

A scan (2026-06) of where an independent researcher with a workstation (16 GB GPU, ample disk)
and a reproducible Airflow-on-Podman pipeline could add something noteworthy. Verdicts are honest:
for an amateur with no telescope time, the realistic unit of contribution is a **value-added
catalogue, a vetted candidate list, or a reproducible methods/benchmark result** — not a discovery.

## Decision

**Next slice: a VLASS Epoch 1+2+3 variability / transient candidate catalogue** (see
`plans/13-vlass-variability.md`). Rationale: catalogue-level (no raw-visibility calibration), a
canonical Airflow batch job, and the published "census of variables" stops at Epochs 1–2 while
Epoch 3 is public and not yet incorporated — a real priority window. The 16 GB GPU is useful later
for ML classification of candidates.

## Citizen-science radio data

| Source | Access | Calibration reality | Verdict |
|--------|--------|---------------------|---------|
| **Radio JOVE** | radiojove.net/archive; MASER CDF | Mostly uncalibrated "SkyPipe units"; no bulk API | **LOW** — the Nancay Decameter Array already did Jovian-DAM statistics over 26 calibrated years; JOVE is a noisier lens on a known result |
| **e-Callisto** (solar bursts, 150+ stations, 20+ yr) | e-callisto.org; `ecallisto-ng` | Uncalibrated ADU at most stations; severe, heterogeneous RFI | **MEDIUM–HIGH** — a complete multi-cycle Type II/III burst occurrence census with completeness corrections does not exist; QC is most of the work; CESRA groups are circling |
| HamSCI/Grape, SuperSID | pswsnetwork, Stanford | GPS-disciplined (Grape) / relative | space-weather, not astrophysics |
| BRAMS/RMOB meteor | brams.aeronomie.be | PI-team priority / heterogeneous | LOW–MEDIUM |
| Amateur pulsar / H-line | — | small aperture | LOW for novel science (contribute GPU to Einstein@Home / NANOGrav PSC instead) |
| **Breakthrough Listen** open data | seti.berkeley.edu/opendata; `blimpy`/`turboSETI` | Professional-grade (GBT/Parkes/MeerKAT) | **HIGH** for non-SETI secondary science (pulsar scintillation, single-pulse re-search) — the best-quality bulk-accessible set |

## Large, under-mined archives (catalogue/image level — Airflow-friendly)

- **ASKAP / CASDA**: RACS three bands (887.5 / 1367.5 / 1655.5 MHz, ~2–3M components each) + **VAST**
  (variables/slow transients, ~500k light curves). TAP via `astroquery.casda`. A homogeneous
  multi-band RACS spectral catalogue and the new RACS-low Stokes V (circular-pol star/pulsar
  candidates) are largely unexplored. **Open.**
- **VLASS** (2–4 GHz, 2.5"): three epochs 2017–2024; CIRADA Quick-Look component catalogues (~3.4M
  components/epoch) free from CADC/CIRADA over HTTP/VO. Published variability stops at E1–E2.
  **Open — chosen.**
- **LOFAR LoTSS DR2/DR3** (144 MHz): DR3 = 13.7M sources over 88% of the northern sky (Jan 2026);
  Stokes V uniquely catches coherent stellar emitters. Catalogues free; raw is 18.6 PB (not local).
- MeerKAT/SARAO (raw visibilities → calibration barrier), FAST (proprietary year; teams have priority),
  GLEAM-X DR2, CHIME/FRB baseband (140 bursts).

## Non-traditional / repurposed datasets

- **Astrogeo VLBI image database** — 139k images of 21k compact sources from decades of geodetic
  VLBI; the largest public VLBI image archive. Free multi-decade dual-band AGN monitoring
  (jets, flux variability, VLBI–Gaia optical-core offsets). **HIGH.** astrogeo.org/vlbi_images
- **Wind/WAVES + STEREO/WAVES** — ~30 yr calibrated solar/interplanetary radio on NASA CDAWeb;
  existing burst catalogues incomplete (Type IV). ML re-classification + CME/X-ray cross-match is
  publishable. **HIGH.**
- **Gaia DR3 × ICRF3 (geodetic VLBI)** radio–optical offsets — both public; **HIGH.**
- Planetary PDS (Cassini RPWS, Juno Waves, Voyager PRA via MASER) — **MEDIUM.**
- GNSS receivers as L-band solar-burst monitors (IGS RINEX) — **MEDIUM** (strongest bursts only).
- Dead ends for astrophysics: SuperDARN/riometers, DSN telemetry, KiwiSDR beyond HF solar,
  generic WebSDR/RFI monitors.

## Signal-from-noise methods and where the GPU helps

- **Image-plane stacking** at multiwavelength-prior positions (eROSITA eRASS1 / Gaia-CRF3 AGN ×
  RACS/VLASS, incl. radio non-detections) — productive, amateur-tractable. Pitfalls that decide
  credibility: snapshot/CLEAN flux bias (~1.4× for sub-threshold sources; Bonnassieux et al. 2023),
  confusion masking, beam mismatch, probabilistic cross-match (NWAY / likelihood-ratio).
- **Classical source finders** (Aegean / PyBDSF / Selavy) remain the baseline (Hydra II, Boyce 2023).
  **ML source finders** (YOLO-CIANNA) beat them but need per-survey retraining.
- **DL denoising / deconvolution** (AIRI, R2D2) is for raw-visibility imaging, **not** post-processing
  already-imaged survey FITS — limited gain for archival cutouts.

**GPU verdict (asymmetric):**
- *Enabling*: incoherent **de-dispersion / single-pulse & FRB search** on filterbanks (Heimdall++,
  Bifrost; 10–100× over CPU); **ML inference** on image cutouts (~500 img/s vs ~10 on CPU).
- *Idle*: image stacking, cross-matching, spectral fitting — network-I/O and CPU bound.

## Shortlist that came out of this (ranked for our setup)

1. Cross-survey spectral + variability value-added catalogue (extends `spectra`) — GPU idle.
2. **VLASS E1+E2+E3 variability / transient catalogue — chosen** — GPU helps classify.
3. Multiwavelength stacking (eROSITA/Gaia × RACS/VLASS) — GPU idle.
4. GPU archival single-pulse / FRB re-search (extends `driftsearch`) — GPU critical.
5. Astrogeo VLBI or Wind/STEREO WAVES mining — high upside, new toolchain.

Full cited agent syntheses were generated 2026-06 and condensed here; key resources: CASDA
(`astroquery.casda`), CADC/CIRADA VLASS, LoTSS portal, Astrogeo, NASA CDAWeb, Breakthrough Listen
open data, `ecallisto-ng`, `blimpy`/`turboSETI`, Heimdall++/Bifrost.
