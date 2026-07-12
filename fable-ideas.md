# fable-ideas.md — deep opportunity re-scan (2026-07-05)

A fresh, ground-up hunt for new/innovative research slices, run as **12 parallel deep-research
agent sweeps** (FRB/time-domain, LPTs/coherent transients, catalogue cosmology/statistics,
solar/heliospheric, planetary/exoplanet, pulsars/PTA, HI/spectral-line, rooftop-station science,
a deliberately speculative wild-card vetter, and a GPU/ML coordinator with six sub-sweeps:
foundation models, SBI, DSP ports, RFI segmentation, FAST-FREX, solar dynamic-spectra ML),
each cross-checked against this repo's 26 merged slices and the 2026-07 opportunity scan.
It **supersedes the shortlist in `survey/opportunity-scan-2026-07.md`** (whose Tier-1 items are
now largely executed: `stokesv_discovery`, `lpt`, `rmstructure`, `torchfdmt`, `junodam`).

**Every idea below is written to be turned into a `plans/NN-*.md`** (Context / Deliverables /
Approach / Verification) by a planning model without re-doing the survey: each entry carries the
gap, the novelty evidence (closest papers, arXiv IDs), the exact data products + access path, a
method sketch naming the reusable modules, a recover-a-known validation, the biggest risk, and
its GATE-0.

## ⚠️ Standing caveat: this scan could not fetch primary sources

The session that produced this file ran behind an egress policy that blocked **all** astronomy
hosts (arxiv.org, ADS, Zenodo, CASDA, CDS, PDS, lofar-surveys.org… — 403 at the proxy CONNECT
stage). Web *search* worked; page/PDF *fetch* did not. Consequences:

1. Every literature claim below rests on **cross-corroborated search snippets** (each claim was
   required to match across ≥2 independent search hits), not full-text reads. arXiv IDs are
   internally consistent but **unverified against the actual abstract pages**.
2. **No data-access path was live-verified in this scan.** Paths marked "verified" inherit that
   status from the repo's earlier scans (`survey/opportunity-scan-2026-07.md`, findings files).
3. Therefore **every idea has an implicit GATE-0**: from the workstation (normal network), do a
   direct ADS/arXiv full-text pass on the cited "closest papers" to confirm the gap is still
   open, and re-verify the data URL, before writing the plan. Ideas whose open/closed verdict is
   *fragile* say so explicitly.

## Workstation profile (design ideas around this)

- **GPU:** AMD RX 7600 XT 16 GiB (gfx1102). PyTorch-ROCm **empirically working** (pinned wheel;
  see the GPU addendum in `survey/opportunity-scan-2026-07.md`): FFT 9×, gather-style
  dedispersion 24× vs CPU. Rule: **pure PyTorch runs; anything shipping CUDA kernels does not**
  (Heimdall, dedisp, AstroAccelerate, CuPy-locked tools, RAPIDS).
- **Disk:** ~275 GB free. **Owner explicitly OK with multi-day/multi-week jobs** — ideas that
  need "compute nobody bothers to spend" (all-pairs statistics over 13.7M-source catalogues,
  long SSL training, week-long forced-photometry sweeps) are in-scope and several below exploit
  exactly that.
- **Auth available:** CASDA (OPAL + `~/.casda_pw`, SODA verified), free-registration archives OK.
  CASDA TAP endpoint is `https://casda.csiro.au/casda_vo_tools/tap`.
- **Rooftop station** (Philadelphia): H-line receiver hardware in acquisition, meteor-scatter
  station planned; first light targeted late summer 2026 (Perseids Aug 12–13). Station ideas are
  hardware-gated — see the Station track at the end.

---

## Corrections & closed doors (checked 2026-07; cite these to avoid dead ends)

Findings from this scan that **invalidate or foreclose** earlier framings — a planning model
should treat these as hard "do not plan" flags unless the GATE-0 re-check overturns them:

- **"First ML burst detector on OVRO-LWA" is GONE.** The OVRO/NJIT team published a real-time
  YOLOv8 type III/IIIb detector, ApJ 1003, 57 (2026), arXiv:2603.25446. What survives: a
  retrospective public *catalogue*, a *type II* detector, type I statistics, open weights
  (F13/F14 below). Update the corresponding entry in `survey/opportunity-scan-2026-07.md`
  (GPU avenue #3) when next touched.
- **CHIME Cat 2 "obvious" population angles are closed**: fluence/DM/scattering debiasing via
  587k injections (arXiv:2606.26334), 80-repeater uniform rate statistics + DM drift
  (arXiv:2605.08410), repeater-vs-one-off spectral split (arXiv:2601.16048), ML repeater
  classification (arXiv:2509.02645, 2512.08308), polarization dichotomy test (arXiv:2401.17378 +
  companion). The *wait-time/duty-cycle* census survives (F2).
- **Saturn SED (lightning) census is closed** — Fischer et al. 2025 (JGR, 10.1029/2024JA033560)
  covers the complete Cassini SED set incl. polarization; Voyager-era closed by Imai et al. 2026
  (10.1029/2025JE009079). Do not plan; only a narrow proximity-duty-cycle companion might
  survive a full-text check.
- **CGM magnetization via RM stacking on Mg II absorbers is closed** — Van Eck/Malik et al.,
  A&A 2026, arXiv:2605.16924 (SPICE-RACS DR2, 612 sightlines). The *group/cluster-halo tracer*
  variant survives (F19).
- **Total-intensity radio source-count dipole is crowded** (arXiv:2509.16732 PRL joint
  NVSS+RACS+LoTSS-DR2 at 5.4σ; arXiv:2509.18689; RMP colloquium arXiv:2505.23526). The **RM
  dipole** (F1) and **polarized-flux dipole** (F18) axes are the open differentiators.
- **RFI segmentation "new architecture on the standard benchmark" is saturated** (~10+ entries
  on the Mesarcik LOFAR/HERA Zenodo benchmark 10.5281/zenodo.6724065 since 2017). Survivors:
  cross-observatory generalization, uncertainty calibration, an open pure-PyTorch inference
  package with public weights (folded into F6).
- **IPTA DR3 is NOT public** (mid-2026; ipta4gw.org says "next few years"). Use NANOGrav 15yr
  (Zenodo 7967584/8060824/8092346) and PPTA DR3 (CSIRO DAP 10.25919/w0nw-jt05, 10.25919/23wj-1d69).
- **Gaia DR4 is NOT out** (ESA roadmap: 2026-12-02). Nothing below depends on it; F21 gets
  sharper when it lands.
- **NANOGrav just closed two obvious PTA angles** (June 2026): per-pulsar solar-wind density
  from chromatic noise (arXiv:2606.28571/2606.28554) and profile-change event scans
  (arXiv:2604.05453). A solar-index *cross-check* of their released series is open but
  overlap-fragile — treat as a watch item, not a slice.
- **LPTs are discovered by image differencing, not folding searches** (GLEAM-X J1627,
  GPM J1839-10, ASKAP J1935+2148 all image-domain). A GPU FFA (F6) targets classical pulsar
  reprocessing, *not* LPT discovery — don't motivate it that way.
- **TESS/Kepler-flare × VLASS/RACS temporal coincidence is claimed** (MNRAS 516:540;
  arXiv:2506.21169; arXiv:2408.14612).
- **FRB DMs as a solar-wind probe: killed by arithmetic**, not literature — expected
  ΔDM ~10⁻³ pc cm⁻³ vs per-burst errors 10⁻²–1; a null is worthless. Revisit only at
  ~10⁻⁴ pc cm⁻³ DM precision.
- **Blind VLASS UCD search is closed** (arXiv:2506.21169: 14,915 Gaia UCDs × 3 VLASS epochs,
  zero BD counterparts). The RACS **Stokes-V** two-epoch companion survives (F22).
- Other kills (reasons in one line): sub-second lensed-FRB echoes (CHIME did it in baseband:
  Leung+2022/Kader+2022); ionospheric-RM-vs-geomagnetic from survey catalogues (pulsar-based
  work does it better); GW-O4 afterglow limits from survey epochs (no suitable event → null
  worthless); comet OH/ISO forced photometry (out of band / ≪µJy); asteroid occultations in
  snapshots (expected events ≈ 0); catalogue-level ORC hunting in LoTSS DR3 (image-level
  searches already active: arXiv:2506.08439, 2510.01999); OH megamaser hunting (needs raw
  visibilities the teams own); MeerTRAP statistics (no public unified catalogue); amateur
  meteor **head-echo** classification (physically inaccessible to forward scatter — reframed in
  S5); Saturn-lightning re-analysis (above); "PulsarNet"/"ffancy" (do not exist — phantom tools).

---

# Tier 1 — highest novelty ceiling, plan-ready (ranked)

## F1. The first rotation-measure dipole/anisotropy test (SPICE-RACS DR2)

- **Gap.** The cosmic-dipole-anomaly literature is entirely source-count/flux based; **no
  published test uses Faraday rotation measures as the tracer** (checked against
  arXiv:2509.16732, 2509.18689, the RMP colloquium 2505.23526, and Mittal & Lewis 2605.27520 —
  none mentions RM). SPICE-RACS DR2 (~2.5–3.4×10⁵ RMs over 87.5% of sky, arXiv:2605.16917) is
  the first catalogue large enough.
- **Data.** SPICE-RACS DR2 RM table, CSIRO DAP `csiro:64891` (~5 GB, no auth) — **already
  fetched and used by the merged `rmstructure` slice**, so data risk ≈ 0.
- **Method.** Reuse `rmstructure.py`'s per-|b| Galactic-floor machinery to form per-source
  extragalactic-RM residuals at |b|>30–45° (DEFROST regime, arXiv:2605.13605); HEALPix-bin;
  fit a dipole in mean |RM| or RM² (likelihood/Bayesian); compare direction to the CMB dipole
  and to the Böhme+ source-count dipole; null distribution from footprint-preserving RA
  scrambles (the footprint (Dec≤+49°) is the dominant systematic — the scramble test is
  load-bearing).
- **Recover-a-known.** Inject a known dipole into a mock all-sky RM catalogue with the real DR2
  footprint/noise; confirm amplitude+direction recovery (mirrors `rmstructure`'s synthetic-screen
  validation).
- **Risk.** No clean kinematic-expectation amplitude exists for RM (interpretation of a positive
  is hard — frame as an isotropy test with honest nulls); dipole field moves fast → GATE-0
  includes a same-week ADS re-search.
- **Fit.** Pure NumPy/healpy, CPU, days. Tooling: `rmstructure.py` ~70% reuse.

## F2. CHIME/FRB Catalog 2: uniform repeater wait-time & duty-cycle census (time-sensitive)

- **Gap.** Cat 2 (arXiv:2601.09399; 4,539 bursts, 83 repeaters) has closed debiasing, rates,
  spectra, polarization, ML classification (see Corrections) — but **no paper computes one
  uniform statistic (Weibull clustering k + exposure-corrected activity-window/duty-cycle) across
  all ~83 repeaters from the public table**. All periodicity results are single-source campaigns
  (20240209A: 2502.11215; 20240114A ~112.9 d; 20220912A: 2604.09098).
- **Data.** `chime-frb.ca/catalog2` CSV/FITS + exposure/sensitivity maps (no auth). **GATE-0:**
  the site 503'd on 2026-07-01 and was unreachable from this sandbox — confirm recovery and
  mirror everything into `data/` immediately.
- **Method.** Extend `frbstats.py` (Weibull) + `frbperiod.py` (Rayleigh Z², LS) into one
  hierarchical pass: per-repeater k, exposure-corrected periodogram, population k-distribution,
  stated burst-count completeness cut. The Cat-1 slice's honest lesson (exposure-blind FAPs are
  not rigorous — `survey/period-findings.md`) is the design driver: use the catalogue's own
  exposure maps.
- **Recover-a-known.** Reproduce **both** published anchors from the same dataset: 20180916B's
  16.35 d and 20240114A's ~112.9 d periods, before trusting the other ~81.
- **Risk.** Highest scoop risk in this file — CHIME/FRB is the natural author. Move fast or not
  at all.
- **Fit.** CPU, days–weeks. `frbstats`/`frbperiod` ~80% reuse.

## F3. SBI population inference for the RACS Stokes-V coherent-emitter class (first anywhere)

- **Gap.** Neural simulation-based inference has hit pulsars (arXiv:2312.14848, 2412.04070),
  magnetars (2503.11875), FRB selection functions (2606.26334), and per-source QU-fitting
  (VROOM-SBI, 2605.27538) — but **no SBI population inference exists for circularly-polarized
  stellar/coherent emitters**. The repo already owns the two hard inputs: a validated forced-V
  measurement (`stokesv.py`) and a per-field leakage-floor model (≈ the selection function).
- **Data.** RACS-low/mid V via CASDA (verified access pattern); SRSC (VizieR `J/other/PASA/41.84`)
  as labelled positives; Gaia CNS5/DR3 M-dwarf parent sample.
- **Method.** Forward model: draw a population (luminosity function slope/break, beaming
  fraction, distances from the Gaia parent sample) → fold through RACS sensitivity + the
  per-beam leakage floor → NPE (the pure-PyTorch `sbi` package, v0.26.x — ROCm-safe; avoid
  LtU-ILI's TF backend) on (LF slope, beaming fraction, break) conditioned on the real
  detections/non-detections from the `stokesv_discovery` target list.
- **Recover-a-known.** Inject synthetic V populations into real RACS cutouts via
  `measure_circular_pol`; verify NPE posterior coverage (simulation-based calibration).
- **Risk.** Model misspecification (beaming/duty-cycle degeneracy at small N) — report posterior
  widths honestly; the deliverable is the first calibrated posterior, not a tight number.
- **Fit.** GPU (NPE training is exactly what the 16 GiB card is for). New dep: `sbi`.
  `stokesv.py`/`stokesv_discovery.py` heavy reuse. Also seeds two siblings: SBI on the
  RM-structure-function turbulence parameters (extends F1/`rmstructure`) and SBI on the LPT
  population using the published GLEAM-X injection-recovery completeness (arXiv:2509.06315) —
  pick one, don't do all three at once.

## F4. Radio counterpart survey of the 56 optically-selected white-dwarf-pulsar candidates

- **Gap.** Pelisoli et al. 2025/26 (MNRAS 540, 821; arXiv:2505.04693) published 56 AR Sco-like
  candidates from Gaia+WISE light curves (26 previously uncharacterized). **No systematic radio
  search of the list exists** — the one confirmed WD pulsar found this route (J1912−4410,
  arXiv:2306.09272) was followed up object-by-object. A detection = a new white-dwarf pulsar;
  the non-detection table is honest value regardless.
- **Data.** Candidate table from the paper's supplementary/VizieR (**GATE-0: locate the
  machine-readable 56-row table**); RACS-low/mid I+V via CASDA; VLASS QL2 (CADC/CIRADA, no auth)
  for Dec>−40 targets.
- **Method.** `measure_circular_pol` forced I/V photometry per epoch per candidate +
  `classify_emitter` leakage vetting + VLASS QL2 cone checks; report detections, |V|/I, and
  upper limits.
- **Recover-a-known.** J1912−4410 is *in* the candidate list and is radio-detected — the pipeline
  must re-find it (same pattern as the GJ 65 control in `stokesv_discovery`).
- **Risk.** Pelisoli's team is the natural competitor (may hold MeerKAT time); table access
  unverified.
- **Fit.** CPU + CASDA I/O, ~a week wall-clock (resumable, like the M-dwarf run).

## F5. Lensed-repeater test: recurring time-delay patterns in CHIME Cat 2 (theory-only niche)

- **Gap.** Theory fully specifies how a strongly-lensed one-off FRB masquerades as a "repeater"
  with a **fixed pattern of mutual delays** (Dai & Lu 2017, 10.3847/1538-4357/aa8873; Li+2018,
  Nat. Comm. 9:3833; review arXiv:2412.01536). **No observational catalogue-level search has
  been published.** CHIME's own lensing searches are intra-burst/microsecond baseband work
  (different regime). Cat 2 is the first dataset big enough.
- **Data.** Same as F2 (per-burst TOAs, DMs, morphologies, exposure functions) — bundle with F2.
- **Method.** Per repeater: all-pairs TOA delay histogram (GPU-trivial), search for recurring
  delay differences; require DM agreement at the measurement floor + morphology consistency;
  FAP from **exposure-aware** TOA scrambles (CHIME's transit cadence aliases — the scramble is
  mandatory); injection tests map the detectable (delay, magnification-ratio) region.
- **Recover-a-known.** Injected lensed pairs recovered at stated FAP; known non-lensed
  high-count repeaters yield clean nulls.
- **Risk.** Expected yield ≈ 0 (optical depth ~10⁻⁴) — the paper *is* the first empirical upper
  limit on the lensed fraction; frame that way from day one.
- **Fit.** CPU/GPU, days. `frbperiod` scramble/FAP machinery ~60% reuse.

## F6. `torch-dsp`: coherent dedispersion + RFI-excision kernels in pure PyTorch (JOSS + science leg)

- **Gap.** Verified per-algorithm (July 2026): **no PyTorch/JAX coherent dedispersion exists**
  (dspsr/CDMT are CUDA-only); **no PyTorch SumThreshold/spectral-kurtosis RFI kernels exist**
  (the field's freshest tool, `jess` — AJ Jan 2026 — is CuPy/CUDA-locked; IQRM/AOFlagger are
  CPU); **no GPU FFA of any kind exists** (riptide is CPU; its 2021 GPU issue #2 is unanswered).
  Prior art fence: `PyTorchDedispersion` (GitHub, 1 star) covers *incoherent* dedispersion only;
  this repo's own `torch-fdmt` covers FDMT.
- **Data (validation oracles).** CHIME/FRB public **baseband release** (DOI 10.11570/23.0029,
  140 known-DM FRBs, HDF5, 2.56 µs) for coherent dedispersion; the vendored Crab `.fil` +
  `jansky.rfi`'s CPU SK/SumThreshold as the RFI oracle; riptide-CPU as the FFA oracle.
  **GATE-0:** confirm one small CHIME baseband file is fetchable from CANFAR without auth.
- **Method.** Coherent dedisp = chirp-filter multiply via `torch.fft` (rocFFT-backed — the 9×
  FFT benchmark on this card de-risks it); SK/SumThreshold = block moment/threshold reductions
  (embarrassingly tensor-parallel); FFA is the stretch goal (budget for the same
  batched-gather vectorization wall `torch-fdmt` hit — honest framing: portability first,
  speed second). Package + benchmark + one real recover-a-known per kernel, JOSS-style.
- **Recover-a-known.** Re-dedisperse a CHIME baseband FRB to its published DM/structure; match
  `jansky.rfi` masks on synthetic + real RFI; re-find the Crab period with the FFA.
- **Risk.** Engineering, not novelty: ROCm/gfx1102 remains officially unlisted — pin and re-test
  torch versions (documented pattern in `torchfdmt-findings.md`).
- **Fit.** The flagship GPU slice; multi-day benchmarks welcome. Extends the `torch-fdmt` paper
  arc into a coherent-DSP suite; each kernel also unlocks later archival slices (F10, F24).

## F7. LPT catalogue v3 + the first multi-epoch Stokes-V forced photometry at all LPT positions

- **Gap.** (a) `data/lpt_sample.csv` (13 objects) predates ≥3 2026 discoveries: ASKAP
  J142431.2−612611 (arXiv:2603.07857), ASKAP J165130.3−450520 and J170036.6−445758 (VASTER,
  arXiv:2606.20067) — a ~30% undercount, and the Rea+2026 review's "no population synthesis"
  flag still stands. (b) The merged counterpart cross-match was **Stokes-I only** (VLASS/LoTSS);
  **nobody has done forced V photometry at the LPT positions across RACS low1/low2/mid** —
  RACS-low2 Paper VIII's blind V catalogue (arXiv:2606.16182) did not target them. A persistent
  V counterpart discriminates coherent-emission models; systematic V limits sharpen the
  burst-only character.
- **Data.** Discovery-paper tables (new rows); CASDA RACS V epochs at ~16–17 positions.
- **Method.** Extend `lpt.py` schema + re-run `crossmatch_counterparts` and the population
  statistics (does N=16–17 move the ~78-min binary-boundary test's p?); then
  `measure_circular_pol` per epoch per position with leakage vetting; report V detections/limits
  + any inter-epoch handedness/sign changes.
- **Recover-a-known.** SRSC positive controls (176/176) + GJ 65 re-detection before touching LPT
  positions; new CSV rows spot-checked against discovery papers (the process that caught the
  review's own data-file typo).
- **Risk.** Likely all-limits outcome at RACS depth (low duty cycles) — the table is the paper.
- **Fit.** CPU + CASDA, days. Near-total reuse of two merged slices (`lpt`, `stokesv`).

## F8. FASHI DR2 (156,411 HI sources): the environment statistics nobody has run

- **Gap.** FASHI **DR2** dropped ~June 2026 (arXiv:2606.31539; 156,411 sources, 19,500 deg²,
  z<0.09 — 4× DR1) with only a *global* Schechter HIMF published. The single environment paper
  (arXiv:2510.22902) used 230 DR1-era group galaxies (~3% of DR2 volume, no cluster-scale
  densities). Open: environment-split HIMF, HI-deficiency vs clustercentric radius (the
  Solanes+2001 curve at 40× the sample), void-vs-wall gas fractions (pre-FASHI ALFALFA-era only).
- **Data.** FASHI DR2 table (VizieR / China-VO mirror — **GATE-0: confirm the catalogue ID**);
  SDSS/DESI group catalogues (Tempel+2017 / Lim+2017 on VizieR); Douglass+2023 void catalogues
  (VizieR, three void-finders).
- **Method.** Position+velocity cross-match → clustercentric radius / local-density terciles /
  void membership; censored (survival-analysis) gas-fraction statistics; Schechter fits per bin.
  Report per void-finder (algorithm dependence is a known trap).
- **Recover-a-known.** Reproduce the declining HI-deficiency profile inside R200 for a
  well-sampled cluster, and the ALFALFA-era "void galaxies are gas-richer at fixed M*" result,
  before any new claim.
- **Risk.** FAST 2.9′ beam confusion at low z (flag blends); the FASHI team is fast — but has
  signalled global-HIMF interests, not cluster work.
- **Fit.** CPU, days. Extends `hi.py`/catalogue muscle memory; biggest *fresh-data* surface in
  this file.

---

# Tier 2 — strong slices, mostly value-added or one-gate-from-ready

## F9. Voyager 2 PRA: modern re-derivation of the Uranus & Neptune radio rotation periods

Lamy+2025 (Nat. Astron. 9, 658) moved Uranus's period by **28 s** using 11 yr of HST UV aurora —
the 1986 *radio* value (17.24±0.01 h) underlying System III was never reanalysed with modern
statistics; Neptune's 16.11 h has no independent check in 28 yr; the Cecconi+2017 PRA
refurbishing (arXiv:1710.10471) covers only Jupiter/Saturn. **Method:** extract burst/flux time
series from the PDS-PPI VG2-PRA Uranus/Neptune encounter volumes (**GATE-0: pin dataset IDs**),
run `frbperiod`'s LS/Rayleigh machinery with honest few-cycle uncertainties; three-way compare
(1986 value, Lamy+2025, this work). **Recover-a-known:** consistency with ~17.24 h/~16.11 h is
itself the pipeline validation. **Risk:** days-long flybys → wide posteriors; the honest paper
may be "modern uncertainties say the radio data cannot distinguish" — still citable.
**Fit:** CPU, small data.

## F10. Parkes Transient Events II: a giant-pulse/heavy-tail census across 363 pulsars in 1.5 GB

PTE-II (arXiv:2508.14403; sqlite DB, **165,592 single pulses from 363 pulsars**, raw segments
preserved) is brand new and its own paper is infrastructure-only. Only a handful of pulsars have
published single-pulse energy-tail fits. **Method:** per-pulsar energy histograms → lognormal vs
power-law-tail model comparison; rank tail index vs Ė (ATNF); exposure-normalize via
injection-recovery (the `stacking` discipline). **Recover-a-known:** Vela / B1937+21-class known
heavy-tail emitters. **Risk:** heterogeneous 1997–2001 sensitivity — normalize or report as
lower bounds. **Fit:** perfect for the GPU (`fdmt`/`singlepulse` reuse for re-extraction),
tiny footprint.

## F11. JBO glitch catalogue: the post-2018 population statistics (≥120 unanalysed glitches)

Catalogue now at **664 glitches / 207 pulsars**; the last population paper (Basu+2022,
arXiv:2111.06835) stops at end-2018 (543/178). Recent bimodality work is single-pulsar
(arXiv:2502.20017). **Method:** scrape `jb.man.ac.uk/pulsar/glitches/gTable.html` (access
pattern already verified in earlier scans); per pulsar with ≥5 glitches, BIC-selected
exponential-vs-mixture waiting-time + size-distribution classification; headline = which pulsars'
classifications *change* with post-2018 data. **Recover-a-known:** Crab/Vela exponential
waiting times; J0537−6910 behaviour. **Risk:** JBO team publishes updates every ~3–4 yr —
check for an in-press paper first. **Fit:** CPU, days; `ppdot`/`lpt` compilation pattern.

## F12. Five-archive simultaneous type III event at cycle-25 max (completes `type3synthesis`)

Exactly one precedent exists (2019-04-15: PSP+STEREO-A+Wind+e-Callisto+EOVSA, arXiv:2306.01910)
— pre-OVRO-LWA, off-max. All five archives are now simultaneously public (e-Callisto qkl tree,
OVRO-LWA spectrograms, Wind/WAVES, PSP/RFS, SolO/RPW). **Method:** time-boxed hunt (say 3
weeks of search effort) for one burst crossing all bands, then run the existing drift-to-distance
ladder end-to-end; deliverable = mutual consistency of five distance/speed estimates on one real
event (or an honest bounded null). New code: only an OVRO-LWA spectrogram reader. **Risk:** the
hunt can fail — pre-commit the time box. **Fit:** CPU; near-total reuse of 4 merged slices.

## F13. OVRO-LWA metric type II census × LASCO CME catalogue (the detector lineage skips type II)

Every type-II-adjacent census is cycle-24, ascending-phase-only, or N=10 case studies
(arXiv:2512.21846; RSTN Solar Phys. 2024), and the published OVRO-LWA detector (2603.25446) is
type-III-only. The repo's own TODO already wants a slower-drift template. **Method:** type II
detector (slow-drift ridge + band-split heuristic) over the OVRO-LWA Level-1 spectrograms
(2024-04→now; portal `ovsa.njit.edu/lwadata-query`, **GATE-0: verify access/format**) and/or
e-Callisto; cross-match LASCO CME v2 + GOES flares; occurrence vs cycle phase (SILSO).
**Recover-a-known:** the established type II–fast-CME association fraction. **Risk:** ~2 yr
baseline limits cycle-phase claims; RFI at 13–40 MHz. **Fit:** CPU/GPU, `solarbursts` +
`ecallisto_census` coverage-correction reuse. Sibling (same data, pick one first): **type I
noise-storm census** — also absent from the literature on this archive.

## F14. First *external* benchmark on FAST-FREX + an open-weights burst classifier (ROCm)

Verified: every published FAST-FREX user shares authors with the dataset/DRAFTS team; DRAFTS
itself is CUDA-locked. A zero-author-overlap, pure-PyTorch benchmark with released weights is
still open — and doubles as the ROCm demonstrator. **Data:** FAST-FREX (Science Data Bank,
doi:10.57760/sciencedb.15070; 600 bursts + 1,000 RFI FITS; **GATE-0: download works from a US
connection; total GB unknown**). **Method:** ResNet/U-Net (plain torchvision) vs the paper's
RaSPDAM baseline; publish weights + training recipe. **Recover-a-known:** the paper's own
metrics as the oracle. **Risk:** modest ceiling (benchmark paper); check no external benchmark
landed in the interim. **Fit:** GPU training, days.

## F15. SSL/anomaly sweep of a survey nobody has swept (RACS continuum, or LoTSS DR3)

Verified map (July 2026): SSL/SOM/anomaly sweeps exist for MGCLS, EMU-pilot, LoTSS DR1/DR2,
VLASS-QL-partial — **none for RACS continuum, none for LoTSS DR3** (5-month window), none for
EMU-main. The recipe is established and ROCm-safe: BYOL/DINO features → Astronomaly-Protege
active anomaly ranking (arXiv:2411.04188, 2602.15930). **Data:** LoTSS DR3 cutout API (no auth)
or RACS cutouts. **Method:** ~10⁵–10⁶ cutouts, SSL pretrain (or linear-probe a public DINOv2 —
arXiv:2409.11175 shows generic ViTs already transfer at F1 0.72–0.88), Protege loop, human-vet
the top-N; deliverable = ranked rare-morphology candidate list + the feature model.
**Recover-a-known:** known ORCs/rings/GRGs in-footprint must rank highly. **GATE-0:** read
STRADAViT (arXiv:2603.29660) + DR3 paper full-text for "in prep" claims; check the unverified
HF checkpoint `ISSA-ML/stradavit-base`. **Risk:** LOFAR/ASKAP teams could be mid-flight;
morphology novelty requires careful de-artifacting (pairs well with F26). **Fit:** the
long-GPU-run slice (owner-approved); weeks OK.

## F16. Broadband technosignature EIRP limits from survey Stokes V (empty haystack cell)

All EIRP-limit literature is narrowband (BL GBT, MeerKAT commensal arXiv:2103.16250, OVRO-LWA
narrowband 2606.04304); Lenc+2018 showed V isolates artificial emitters but nobody converted
survey-V non-detections into technosignature limits. **Method:** forced V at ~10⁴ Gaia CNS
(100 pc) positions in RACS-low/mid (`measure_circular_pol` at scale); 3σ V limits →
EIRP = 4πd²·S·Δν; place on the haystack axes vs narrowband limits; exclusion list from
`stokesv_discovery`. **Recover-a-known:** known V stars + satellite artifacts must appear.
**Risk:** referee pushback on broadband-transmitter priors — frame strictly as parameter-space
cartography; leakage sets the floor. **Fit:** CASDA-I/O-bound, resumable, ~week; near-total
`stokesv` reuse. RNAAS/arXiv-able.

## F17. e-Callisto as an accidental 15-year RFI observatory: the megaconstellation trend

The only RFI-trend precedent stops in 2019 — the year Starlink began (Pérez+2020, SoPh 295:11);
dedicated-instrument satellite-RFI studies (LOFAR 2307.02316, SKA-Low 2506.02831) don't touch
archival spectrograph records. **Method:** robust per-station/per-channel daily occupancy metric
(quantile-based, solar-burst-immune) for 10–20 configuration-stable stations, 2012–2026; trend
the documented unintended-emission bands; attribute via Space-Track TLE pass windows.
**Recover-a-known:** reproduce Pérez+2020's Spanish-site increase; FM/DAB fixed transmitters as
controls. **Risk:** station hardware changes masquerade as trends (differential in-station bands
mitigate). **Fit:** the Airflow ingest (`ecallisto_pipeline`) already streams this archive —
infra reuse is total. Null is still a citable spectrum-management result.

## F18–F20. Cosmology/statistics cluster (siblings of F1; do after it, same data/muscles)

- **F18 — Polarized-flux dipole (SPICE-RACS DR2).** Last done with NVSS (Tiwari & Jain 2013/15);
  never with a modern RM/pol catalogue. Reproduce the NVSS-era result first (recover-a-known),
  then DR2. Risk: sparse counts after S/N cuts; scan-pattern leakage systematics.
- **F19 — CGM RM² stacking on DESI group/cluster halos.** The Mg II angle is closed
  (2605.16924); the halo-mass/group tracer is not. Re-implement their annulus-GRM subtraction
  (recover-a-known = reproduce their Mg II excess), then stack behind DESI DR1 groups by
  richness. Risk: small N behind massive groups; don't spin a null as contradiction.
- **F20 — First real-data application of the ⟨RM²×g⟩ estimator** (Zhang & Lidz, arXiv:2512.06584
  — SKA forecast only). Likely an upper limit at DR2 size (30× fewer RMs than forecast) — say so
  up front. Hardest estimator of the three; stretch goal.

## F21. Secular aberration drift from public ICRF3×Gaia DR3 (cleanest quick win)

Reproduce the ~5 µas/yr Galactic-centre glide (Klioner+2021; refinements 2503.03389) from the
*public* cross-match with a VSH/glide fit — an independent-pipeline reproducibility note that
extends the merged `offsets` slice with near-zero new access. Recover-a-known *is* the result.
Risk: Gaia AGN PM noise floors — quote the discrepancy honestly. Gets 10× sharper with DR4
(2026-12). **Fit:** CPU, days; `offsets.py` ~90% reuse.

## F22. RACS-V two-epoch census of the Gaia UCD sample (companion to the VLASS null)

arXiv:2506.21169 closed the VLASS-I blind search (0/14,915 UCDs). The V-selected, lower-frequency
two-epoch companion is open — UCD emission is coherent/highly-V (new VLITE detection 2512.11120;
Driessen review 2606.27706 flags exactly this lever). **Method:** `measure_circular_pol` at the
Gaia UCD positions in RACS-low1/low2. **Recover-a-known:** WISE J0623 (arXiv:2306.15219) +
Pritchard UCD detections. **Risk:** likely near-zero detections — the limit table is the paper.
**Fit:** merges naturally with F16's sweep (same cutout I/O — do them in one pass).

## F23. Cassini SKR: post-2013 + Grand-Finale proximity census (port `junodam` to Saturn)

Fischer+2015 stops early-2013; nobody has run the RPWS-flux dual-period census through the 2017
proximal orbits, and nobody has asked "does the junodam ~180× proximity law hold for SKR".
**Data:** PDS-PPI `CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0` (60-s key params; **GATE-0: pre-2008
volume ID**); Horizons geometry. **Method/validation:** re-derive the known 10.8/10.6 h split
first; then range-binned duty cycles + LS through 2017. **Risk:** SKR "proximity" needs
magnetic-latitude weighting (real modelling step). **Fit:** CPU; `junodam` + `frbperiod` reuse.

## F24. Apertif TD DR2 single-pulse reprocessing (torch-fdmt's real-data leg 2)

DR2 (1,666 pointings, ~0.48 GB/pointing 1-bit Stokes-I PSRFITS) was FRB-searched (AMBER/ALERT,
arXiv:2406.00482) but has **no published archive-wide single-pulse RRAT/known-pulsar census**.
**GATE-0:** ASTRON tape-staging request turnaround (free helpdesk). **Method:** `fdmt.py` +
`singlepulse.py` over staged pointings; ATNF position-matched folding; redetections + candidate
trains. **Recover-a-known:** ≥1 known pulsar/RRAT per staged field. **Risk:** staging stalls;
modest ceiling — frame as tooling validation at scale. **Fit:** GPU; disk-aware batching.

## F25. Merged GPPS+CRAFTS pulsar catalogue (the `lpt` pattern at N≈1,000)

751 GPPS + >200 CRAFTS pulsars exist only as per-paper tables; no homogeneous provenance-typed
merge with Z-heights, luminosity-function placement, MSP fraction vs ATNF. Scrape
`zmtt.bao.ac.cn/GPPS/GPPSnewPSR.html` (live bulk table) + paper tables; `ppdot`/`pulsarspec`
reuse. Recover-a-known: GPPS's own faint-luminosity-tail claim. Risk: transcription errors at
N≈1,000 — keep the `lpt` flag discipline. **Fit:** CPU, compilation slice.

## F26. LoTSS DR3 artifact archaeology: exhaustive offset-vector ghost/duplicate census

Nobody publishes exhaustive small-separation pair statistics on mega-catalogues (QA sections
sample). GPU-chunked all-pairs within 10′ over 13.7M sources; cluster offset vectors in the
bright-source frame and mosaic-tile frame; ghost families = repeated offset/PA excesses.
Deliverable: "DR3 contains <X% duplicates above S" + a flag list — citable by every DR3 user
including F1/F15. Recover-a-known: injected duplicates + documented deblending modes.
**GATE-0:** read the DR3 QA section first, scope against it. **Fit:** the pure
compute-for-days GPU idea; pure torch.

## F27. Blind moving-source search across VLASS E1/E2/E3 (radio proper motions without Gaia)

All published radio PMs start from a known object (Gaia-anchored: arXiv:2409.18466). GPU
all-pairs E1–E2 linkage (0.3–5″/yr) + E3 collinearity confirmation + flux/morphology cuts,
*then* exclude Gaia/WISE counterparts — the survivors are optically-dark movers (Y dwarfs,
high-PM pulsars) or a first surface-density limit. Recover-a-known: UV Ceti (3.4″/yr) must fall
out blind. Risk: VLASS astrometric floor (0.25–0.5″) and variable-AGN false pairs — collinearity
+ scramble tests are load-bearing. **Fit:** GPU all-pairs, weeks OK.

## F28. WALLABY DR2 pair: BTFR/angular momentum vs environment + cube stacking at DESI positions

Still open (confirmed: no DR3 until ~2029; the 126 kinematic models are the ceiling; team papers
since are morphometrics/dark-sources, not BTFR-environment). Two legs: (a) BTFR + Fall-relation
residuals vs field/density (N=126 — pre-register the null framing); (b) sub-threshold spectral
stacking of WALLABY cubes at DESI-z positions, velocity-aligned (adapt the `stacking`
injection-recovery harness from continuum to cubes; recover catalogued sources first).
**GATE-0:** fetch the kinematic products (CSIRO DAP DOI 10.25919/7w8n-9h19 — team page, not
obscore) and check DESI footprint overlap. **Fit:** CPU, ~500 MB/field cubes.

---

# Tier 3 — quick wins, fillers, watch items

- **F29. Cycle-25 prediction test (fastest honest note).** Compagnino & Zuccarello 2021
  (arXiv:2103.13699) made a falsifiable cycle-25 radio-loud/halo-CME rate prediction from
  cycle-24 correlations; nobody has checked it against the real cycle now past max. LASCO CME
  v2 + SILSO + RSTN lists; near-zero new code (`ecallisto_census` statistics). Recover-a-known:
  refit their cycle-24 coefficients first.
- **F30. MALS DR3 Galactic HI-absorption demographics** (3,640 features / 19,130 sightlines,
  arXiv:2504.00097, own portal `mals.iucaa.in`): covering fraction vs |b|/ℓ, optical-depth
  distribution vs 21-SPONGE. Low ceiling, zero drama.
- **F31. Fermi 4FGL-DR4 × ATNF from-scratch positional+Ė/d² cross-match** → a small vetted
  "radio pulsar in unassociated gamma-ray ellipse" list. Must recover all ~294 known gamma
  pulsars first. `ppdot` reuse.
- **F32. MeerKLASS 2019 IM cube × DESI** (WiggleZ-only precedent, 7.7σ, arXiv:2206.01579):
  independent-tracer check. Foreground cleaning is the hard part — follow the published
  pipeline, frame conservatively.
- **F33. LoTSS DR3 faint-count injection spot-check** — an independent completeness harness
  (from `stacking`) vs the DR3 paper's own counts; also the natural *migration of plan 29* to
  DR3. Reproducibility note, not discovery.
- **F34. PINT quick-look toolkit** (residuals → normalized PSD → step-function/small-glitch
  matched filter → EFAC/EQUAD triage) on NANOGrav 15yr + PPTA DR3; NGC 6440E vendored.
  Methods/JOSS-adjacent; gate any "candidate glitch" against chromatic-noise nulls.
- **F35. e-Callisto fine-structure census** (spikes/J/U-bursts): network-scale sub-2-s detector;
  independent-instrument comparison to the LOFAR spike-pair result (Nat. Comm. 2026). Report
  cadence-limited completeness honestly.
- **F36. Juno/Waves HOM occurrence + moon-induced-emission incidence table** (same CDFs as
  `junodam`, re-banded 0.3–3 MHz; Louis+2023 encounters as positive controls). Crowded
  LESIA space — scope narrowly as the occurrence complement to their physics papers.
- **F37. Ground (Nançay JunoN) × Juno simultaneous Io-DAM census** — extends `junodam`'s
  "Earth boxes don't organise Juno's sky" finding with an independent ground stream.
  **GATE-0:** file-level access to the NDA/JunoN public products.
- **F38. RCW 103 (6.67-hr ultra-slow magnetar) radio counterpart check** — one-target
  `stokesv`/`lpt` appendix; clean expected-null closing a never-stated gap.
- **F39. EDGES raw-data averaging/flagging robustness check** (scoped: hold their calibration
  fixed, vary day-selection/flagging; recover their averaged spectrum first). NOT a
  recalibration — say so loudly. High effort; park unless keen.
- **Watch items (do not start):** DSA-110 NSFRB LPT-recovery test (their paper names
  CHIME J1634+44 as planned future work — be ready if they don't); NANOGrav solar-wind ×
  F10.7/SILSO cross-check (gate on their data release not already containing it); CRAFT/ASKAP
  scintillation-bandwidth population table (read Scott+2025/2510.05654 first — gap unverified);
  MWA GPM candidate cross-validation (their survey paper in prep will supersede).

---

# Station track (hardware-gated; sequence with the build, not the archive work)

From the station-science sweep — ranked by time-to-citable-result. Precedent checked; the
common finding: the amateur literature stops at qualitative detections, so **rigor (uncertainty
budgets, public-survey ground truth, geometry correction) is the novelty axis**, matching this
repo's house style.

- **S1. LNA-first vs filter-first noise-figure/overload A/B** — bench (VNA + Y-factor +
  out-of-band interferer); folk knowledge nobody has published measured. RNAAS-able methods
  note; do during commissioning. *Updated 2026-07-11 for the integrated Discovery feed:* the
  station's sealed feed fixes its own ordering, so the A/B is a pure methods note on a dedicated
  ~$85 discrete-parts buy (SPF5189Z + standalone 1420 filter + SAWbird+ H1 as the bench-measurable
  integrated reference), anchored by an on-sky ground/cold-sky Tsys measurement of the Discovery
  feed itself (not bench-injectable — antenna-integrated input). See `plans/77`.
- **S2. Open H-line pipeline with recover-a-known vs LAB/HI4PI** — per-pointing calibrated
  spectra reproducing survey profiles to stated accuracy; the CHART/campus-telescope genre
  (arXiv:2307.11173, 2404.17893) ships no tested, CI'd pipeline. `hi.py` LAB reader extends.
  Also the calibration substrate for S3/S4.
- **S3. Annual ±142 kHz H-line Doppler sinusoid with a real uncertainty budget** — the flagship;
  best precedent found is a two-date proof of concept, never 12 months with an error budget.
  The limiting term is gain/bandpass drift (weekly cal cadence + possibly the Dicke-switched
  second dish). Month-12+ deliverable; residual-vs-HI4PI per pointing is the falsifiable form.
- **S4. Urban drift-scan strip vs HI4PI: beam-forward-modelled residual/systematics note** —
  distinct from the many "look, a map" write-ups only if the beam convolution is modelled before
  attributing residuals; one well-characterized declination strip suffices.
- **S5. Meteor shower flux (Perseids/Geminids) vs Global Meteor Network optical flux** —
  geometry-corrected single-station radio rates against GMN's calibrated km⁻²hr⁻¹ product +
  RMOB; the honest ML add-on is ping/aircraft/RFI spectrogram classification (head-echo work is
  physically out of reach — dropped). First shower data ~4 weeks after station start.
- **S6. Meteor influx → sporadic-E coupling (Wallops digisonde, ~200 km)** — daily lagged
  correlation of station meteor counts vs public GIRO foEs + Madrigal TEC, controlling season/
  tides/geomagnetics; literature is professional-radar and seasonal-scale only. ≥6-month
  campaign; TLE-based satellite-glint excision is a by-product.
- **S7. Two-dish additive-interferometer solar diameter at 1420 MHz** — citable if the
  fringe-spacing error budget is honest; **GATE-0 is analytic**: required baseline for the ~0.5°
  disc before buying hardware. Sequence last (per `station/interferometry.md`).
- **S8. Calibrated 21-cm quiet-Sun monitor vs F10.7/RSTN** — only worth it flux-calibrated (Cas
  A/Cyg A transit); otherwise blog-grade. Flare-coincident burst = opportunistic bonus, not a
  deliverable.

---

# Suggested first moves (opinionated)

1. **F1 (RM dipole)** — data already on disk, tooling merged, field-defining question, zero
   archive risk. The natural next `plans/38`.
2. **F2+F5 (CHIME Cat 2 bundle)** — the moment `chime-frb.ca` serves the tables, mirror and go;
   most time-sensitive pair in the file.
3. **F6 (torch-dsp)** — the GPU flagship; extends the repo's strongest software arc
   (`torch-fdmt` → coherent-DSP suite) with three verified-empty niches.
4. **F4 (WD-pulsar radio sweep)** — highest single-object discovery ceiling per unit effort once
   the 56-row table is in hand.
5. **F8 (FASHI DR2)** — the freshest large dataset with textbook-shaped open questions.
6. Keep **F29** in a back pocket as the one-day honest-note palate cleanser.

Cross-checks before any plan is written: (i) the per-idea GATE-0; (ii) the standing full-text
novelty pass (top of file); (iii) `survey/opportunity-scan-2026-07.md` for infrastructure
corrections (CASDA TAP endpoint, e-Callisto tree moves) that still apply.
