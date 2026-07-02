# Opportunity scan 2026-07 — corrected assumptions & new research avenues

A deep re-survey (2026-07-01) of the project's assumptions and the open-data landscape, run as four
parallel agent sweeps: (1) live access verification of 14 candidate sources, (2) a 2025–mid-2026
literature-gap scan across 10 topics, (3) a repo-internal sweep of every self-flagged "future work"
item, and (4) an adversarial re-check of every LOW/NO-GO verdict in `data-source-scan.md` /
`new-findings-scan.md`. **This document is the seed for future `plans/` files** — each Tier-1/2
entry below is written to be turned into a plan (Context / Deliverables / Approach / Verification)
with its GATE-0 already identified.

## Headline: three project assumptions are now WRONG

1. **"Reliable no-auth data sources are largely used up" (TODO.md) — FALSE.** Since the June scan:
   LoTSS **DR3** (13.7M sources, 88% of the northern sky, no-auth FITS + live HTTP **cutout API**),
   CHIME/FRB **Catalog 2** (4,539 bursts / 3,641 sources / 83 repeaters — 8.5× Cat 1),
   **MeerKAT image/catalogue products** (MIGHTEE DR1, MGCLS, MALS DR2, SMGPS, MeerKLASS DR1 —
   no visibility calibration needed), VLASS **Single-Epoch images** for E1+E2 (unblocks stacking),
   **SPICE-RACS DR2** (~3×10⁵ RMs, 5× every previous RM catalogue combined), plus brand-new
   archives the old scan never listed: **OVRO-LWA public solar archive** (calibrated 13–87 MHz
   interferometric dynamic spectra, Aug 2025), **Apertif time-domain DR2** (0.8 PB, 1,666
   pointings), **Parkes PSRDA** (4.5 PB / 34 yr, formal data paper arXiv:2511.22702), and the
   **Juno/Waves flux-density CDFs** (Oct 2025, MASER doi:10.25935/6jg4-mk86).
2. **"16 GB GPU available for CUDA tools" — WRONG.** The workstation GPU is an **AMD Radeon
   RX 7600-class (Navi 33) with no ROCm installed** (verified `lspci`; `nvidia-smi`/`rocm-smi`
   absent). CUDA-only radio tools (Heimdall, GPU turboSETI) will not run out-of-the-box. Plans
   should be **CPU-first** (plan 28 already is, via `jansky.transients`); a ROCm/PyTorch-HIP setup
   is a separate infrastructure task if ever needed. Disk: **275 GB free**.
3. **"CASDA was the blocker for Stokes-V work" — RESOLVED, and the discovery half was never
   built.** `stokesv` (merged) is the *methods* paper. The two discovery angles the survey ranked
   "highest genuine-discovery ceiling" remain unclaimed in the literature (verified against
   RACS-low2 Paper VIII, arXiv:2606.16182, 2026-06): forced V photometry on a curated M-dwarf
   target list, and RACS-low1 vs low2 two-epoch V variability. CASDA auth + SODA cutouts are
   verified working. **Note: the CASDA TAP endpoint moved** to
   `https://casda.csiro.au/casda_vo_tools/tap` (old `/casda_tap/` 404s) — update scripts.

## Verdict flips vs the June scan

| June verdict | Re-check (2026-07) |
|---|---|
| MeerKAT: skip (calibration barrier) | **FLIPS** — 5+ public image/catalogue products; MEDIUM–HIGH for catalogue-level work |
| CHIME public = Cat 1 (~500) + 140 baseband | **FLIPS** — Cat 2: 4,539 bursts, 83 repeaters, exposure/sensitivity maps included |
| LOFAR: only catalogues usable | **FLIPS** — DR2/DR3 HTTP cutout API + LoLSS SIAP live; DR3 = 13.7M sources, no auth |
| VLASS stacking NO-GO (QL dirty-residual bias) | **PARTIAL FLIP** — E1+E2 Single-Epoch images public at CADC → two-epoch SE stacking is GO; E3 SE pending |
| FAST: skip (proprietary) | **PARTIAL FLIP** — archival releases public (GPPS lists, FAST-FREX FRB set on Science Data Bank); independent groups publish from it |
| Radio JOVE: LOW | HOLDS — but Juno/Waves public CDFs open a ground+in-situ sub-case (MEDIUM-LOW) |
| Planetary PDS: MEDIUM | HOLDS, sharpened — two concrete slices identified (Juno DAM census, Cassini SKR periodicity) |
| BL open data: HIGH + GB–TB gating | HOLDS — no lightweight index; secondary science actively published by outsiders |
| GNSS solar monitors: MEDIUM | HOLDS — still no public IGS burst catalogue |
| Amateur dish pulsar/H-line: LOW | HOLDS — archival compute on Parkes/FAST/BL is the right path |

## Live-verified data access (2026-07-01)

| Source | Status | Smallest tractable unit | Auth |
|---|---|---|---|
| LoTSS DR3 | **verified** | `lofar-surveys.org/public/DR3/catalogues/LoTSS_DR3_v1.0.srl.fits` (13.7M rows); cutout API `dr3-cutout.fits` | none |
| PSP FIELDS RFS (LFR/HFR) | **verified** (fetched) | daily CDF ~3.3 MB, `spdf.gsfc.nasa.gov/pub/data/psp/fields/l2/rfs_lfr/` | none |
| SolO RPW HFR/LFR-SURV | **verified** (fetched) | daily CDF ~57 MB, `spdf.../solar-orbiter/rpw/science/l2/hfr-surv/` | none |
| VLASS SE images (E1,E2) + CIRADA QL2 catalogue | **verified** | QL2 component FITS (~3.4M comps); SE via CADC `collection=VLASS`, `SEIP` products | none |
| NANOGrav 15yr par/tim + posteriors | **verified** | Zenodo 7967584 / 8060824 / 8092346 | none |
| EDGES raw + averaged | **verified** | two CSVs (8+21 KB) from Bowman 2018; raw = 290×~3 GB `.gsh5`, open HTTP | none |
| GLEAM-X DR1/DR2 | **verified** | VizieR `VIII/113` / `VIII/114` cone queries | none |
| EMU Pilot Survey | **verified** | CSIRO DAP `csiro:52997` (images + ~220k-source catalogue) | none |
| RMTable consolidated (Van Eck+2023) | **verified** | VizieR `J/A+A/671/A151` or CIRADA-Tools/RMTable FITS (~55k RMs) | none |
| JBO glitch catalogue | **verified** | HTML table `jb.man.ac.uk/pulsar/glitches/gTable.html` (scrape; no CSV) | none |
| WALLABY (CASDA) | **verified** | HI cube ~500 MB via datalink; **DR2 kinematic models are paper/team-page supplementary, not in obscore** | OPAL (free) |
| POSSUM / SPICE-RACS DR2 | open-unverified | pilot cubes in CASDA; **DR2 RM table access is a GATE-0** (see Tier 1 #2) | OPAL (free) |
| MeerKAT products (MIGHTEE/MGCLS/MALS/SMGPS) | open-unverified | SARAO archive UI (JS-rendered, no TAP); paper-supplementary CSVs | SARAO (free) |
| CHIME/FRB Cat 2 | open-**blocked now** | `chime-frb.ca` returning **HTTP 503 persistently**; CANFAR paths 403 | none |
| DSA-110 archive | open-unverified | JS-rendered table at `code.deepsynoptic.org/dsa110-archive/`; check `github.com/dsa110` for CSV | none |
| VCSS 340 MHz | **gated** | CIRADA bright-source catalogue only; full catalogue needs NRL collaboration | partial |
| e-Callisto curated burst lists | **moved** | `qkl/` spectrogram tree verified live (2000–2026); `e-callisto.org/Activity_Reports/` **404s** — the typed II/III/IV lists moved; GATE-0 hunt needed | none |

## Tier 1 — discovery-ceiling opportunities (each one plan-ready)

### 1. RACS Stokes-V discovery slice (forced target-list photometry + two-epoch V variability)
- **Gap (verified open):** RACS-low2 Paper VIII (arXiv:2606.16182) published the *blind* V
  catalogue (61 stars / 85 pulsars) but did **no** forced photometry on a curated M-dwarf/UCD
  input list (pushes below the blind 5σ threshold) and **no** RACS-low1(2020) vs low2(2025)
  two-epoch V variability census (~5-yr baseline, not even for its own 61 stars).
- **Data:** CASDA (verified; auth + SODA working; new TAP endpoint). Gaia DR3 + SIMBAD for the
  target list. **Tooling:** `stokesv.py` forced photometry + leakage floor reuse near-verbatim.
- **Amateur unit:** new circularly-polarized star candidates + a two-epoch V variability table.
- **Risk:** ASKAP team could publish either angle within ~a year — this is the most
  time-sensitive item. **GATE-0:** confirm RACS-low1 V images are servable via SODA like low2.

### 2. SPICE-RACS DR2 RM structure function vs Galactic latitude
- **Gap:** DR2 (arXiv:2605.16917, 2026-05) = 2.5–3.4×10⁵ RMs over 87.5% of sky — 5× every prior
  catalogue combined — released **without** a systematic Galactic structure-function analysis by
  |b|/longitude sector. Direct 10×-data upgrade of our `rmsky` (Taylor 2009, 37k RMs) slice; the
  team's DEFROST (arXiv:2605.13605) handles the extragalactic subtraction as a citable method.
- **Data GATE-0:** locate the DR2 machine-readable RM table (CASDA or PASA supplementary — the
  access sweep verified only DR1-era tables; do not confuse with VizieR `J/MNRAS/519/5723`, which
  is O'Sullivan+2023 LoTSS). **Tooling:** `rmsky.py` structure-function code extends directly.

### 3. LPT period–activity population catalogue (the missing P–Ṗ diagram for a new class)
- **Gap (explicitly flagged by the field):** the Rea, Hurley-Walker & Caleb 2026 review
  (arXiv:2601.10393) notes **no population synthesis exists** for the ~15 confirmed long-period
  transients. A systematic table (periods, activity fractions, spin-down limits, binary IDs,
  multiwavelength counterparts) + a P–Ṗ-style class diagram + VLASS-3-epoch/LoTSS-DR3-cutout
  counterpart checks is a value-added catalogue paper with zero new observations.
- **Data:** published tables + VLASS QL2 catalogue (verified) + LoTSS DR3 cutout API (verified).
- **Tooling:** `ppdot.py` diagram machinery + `vlass.py` cutout confirmation reuse. Very high
  honesty ceiling: compilation, not discovery claims. **GATE-0:** none — all inputs public now.

### 4. CHIME/FRB Catalog 2 repeater wait-time population (80 repeaters, one hierarchical model)
- **Gap:** the collaboration closed fluence/DM/scattering debiasing (arXiv:2606.26334) but a
  **uniform Bayesian wait-time / activity-window analysis across all 80 repeaters** (the
  FRB 20180916B periodicity question, asked of the whole population) is unpublished.
- **Data risk (live):** `chime-frb.ca` is returning 503 — GATE-0 is simply the site recovering
  (or the CANFAR mirror). **Tooling:** `frbstats.py` Weibull wait-times + `frbperiod.py`
  periodograms extend directly. **Risk:** CHIME-team speed is the main competition.

## Tier 2 — solid value-added slices (lower ceiling, high fit)

5. **e-Callisto/OVRO-LWA metric-band multi-cycle burst census** — the DH-band cycle-23–25
   comparison exists (arXiv:2409.02554); the **metric-band** census vs sunspot number does not.
   Our `ecallisto_census` statistic is built and validated; the real ingest showed raw-detection
   is coverage-limited (168-day probe → 5 events, both at 2014 cycle-max, results in
   `results/ecallisto_census_realdays.csv`). Two data paths, both GATE-0s: (a) find the moved
   e-Callisto curated burst lists (typed II/III/IV; `Activity_Reports/` 404s), (b) the new
   **OVRO-LWA public solar archive** (calibrated, Aug 2025+ — too short alone, ideal validation).
6. **WALLABY DR2 BTFR vs environment** — 126 public kinematic models (DR2, Feb 2025); the
   cluster-infall-vs-field BTFR-scatter question is unpublished. Extends `hi`. GATE-0: fetch the
   kinematic-model products from the team page (not in CASDA obscore).
7. **VLASS E1+E2 Single-Epoch stacking upgrade** — the old NO-GO was QL-specific; SE images are
   public at CADC. Re-run `stacking` on SE cutouts with the existing injection-recovery harness;
   quantifies the QL bias as a bonus result. GATE-0: SE cutout service granularity at CADC.
8. **Juno/Waves Jovian DAM occurrence census** — the Oct-2025 calibrated flux-density CDFs
   (MASER, 2016–2019+, 1-s cadence) are new; an Io-phase/CML occurrence census cross-checked
   against the Nançay 26-yr ground statistics reuses our whole CDF/dynamic-spectrum toolchain
   (`swaves`/`windwaves`). New domain, same physics muscle memory.
9. **Cassini RPWS full-mission SKR periodicity** — 60-s key-parameter CDFs on PDS (verified
   dataset ID `CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0`); Lomb–Scargle N/S dual-period census
   including the post-2012 seasons. Reuses `frbperiod` periodogram + CDF tooling.
10. **Pulsar-timing domain cluster** (4 backlog gaps, data vendored/verified): residual
    quick-look (PINT wrapper), red-noise PSD, small-glitch matched-filter scan, PTA noise budget
    — NGC 6440E vendored, NANOGrav 15yr Zenodo verified. Plus a **JBO live-catalogue glitch
    mixture-model** (post-2018 glitches unanalysed in print; HTML scrape).
11. **Plan 28 (BL single-pulse) & Plan 29 (LoTSS deep counts)** — both remain shovel-ready;
    29 should **migrate to DR3** (13.7M sources, better rms map). 28 is CPU-viable via
    `jansky.transients` (GPU caveat above); GATE-0 unchanged (find one small BL filterbank).

## Tier 3 — quick internal wins (from the repo's own future-work list)

- **Type II (slower-drift) template** for the e-Callisto detector — enables type II census.
- **Wind RAD1 + RAD2 combined** — extends the `windwaves` beam track toward 1 AU.
- **Simultaneous e-Callisto × Wind × STEREO single event** — the type3synthesis GATE-0 hunt.
- **DM→z pipeline with Macquart scatter** (backlog) — small, clean, fills the stale `fruitbat`.
- **Southern SED upgrade**: RACS-resolution compactness cuts + the 0.23–0.89 GHz gap-band point
  via CASDA (now unblocked); northern measured-turnover companion via LoTSS DR3 + VLASS.
- **SETI drift-search CLI + injection-recovery benchmark** on the vendored Voyager `.h5`
  (four specced-but-unbuilt uses of an already-registered dataset).
- **Blind V-population and multi-epoch leakage-floor V survey** (stokesv paper's own next step).

## Infrastructure corrections (apply before any new slice)

- CASDA TAP: use `https://casda.csiro.au/casda_vo_tools/tap`.
- CHIME/FRB: treat `chime-frb.ca` as flaky (503); mirror any fetched catalogue into `data/`.
- GPU: **superseded by the GPU addendum below** — PyTorch-ROCm is now installed and empirically
  verified on the card; CUDA-only tools remain walled (see addendum).
- WALLABY/SPICE-RACS/MeerKAT: catalogue products often live in **paper supplementary material
  or team pages, not the observatory archive** — every plan needs a data-location GATE-0, not
  just an archive name.
- e-Callisto quicklook tree verified at `soleil.i4ds.ch/solarradio/qkl/{YYYY}/{MM}/{DD}/`
  (1978–2026); the raw-FITS `/data/` subtree and Activity Reports have moved — re-scout paths.

## Ranked shortlist (what to plan first)

1. **RACS Stokes-V discovery slice** (Tier 1 #1) — highest ceiling, most time-sensitive, tooling ready.
2. **LPT population catalogue** (Tier 1 #3) — zero data risk, field-flagged gap, honest by construction.
3. **SPICE-RACS DR2 structure function** (Tier 1 #2) — 10× data upgrade of an existing slice.
4. **CHIME Cat 2 repeater wait-times** (Tier 1 #4) — gated only on the site recovering.
5. **Juno DAM census** (Tier 2 #8) — new domain, full toolchain reuse, low competition.

Key references: arXiv:2606.16182 (RACS VIII), 2605.16917 (SPICE-RACS DR2), 2601.10393 (LPT
review), 2601.09399 (CHIME Cat 2), 2605.08410 (80 repeaters), 2606.26334 (CHIME debiasing),
2602.15949 (LoTSS DR3), 2409.02554 (DH type II cycles), 2511.22702 (Parkes PSRDA),
2512.11964 (MeerKLASS), doi:10.25935/6jg4-mk86 (Juno/Waves flux density). Agent syntheses
generated 2026-07-01; condensed here.

---

# GPU addendum (2026-07-02) — ROCm on the RX 7600 XT: verified working; re-ranked GPU avenues

The scan's "no ROCm → CPU-first" correction is itself now corrected. **Empirical result (this
machine, 2026-07-02):** the PyTorch **2.12.1+rocm7.1** wheel (`pip install torch --index-url
https://download.pytorch.org/whl/rocm7.1`, venv `~/.venvs/rocm-test`) detects the card **natively
as gfx1102** — no `HSA_OVERRIDE`, and **no system ROCm install needed** (the wheels bundle the
runtime; `/dev/kfd` + `/dev/dri/renderD128` are world-writable on Fedora 44). Measured, with
GPU/CPU results verified matching:

| Kernel (radio-astronomy proxy) | CPU | GPU (RX 7600 XT, 16 GiB) | Speedup |
|---|---|---|---|
| matmul 4096² fp32 (ML inference) | 216 ms (0.64 TF) | 65 ms (2.1 TF) | **3×** |
| 2-D FFT 4096² complex64 (dynamic spectra / imaging) | 50 ms | 5.4 ms | **9×** |
| Dedispersion-style gather+sum, 1024 DM × 4096 chan × 8192 samp | 156 s | 6.5 s | **24×** |

Caveat kept honest: a web survey of the ecosystem (same date) found the *official* ROCm support
matrix still omits gfx1102 and older wheels (≤ torch 2.6) dropped it — the empirical wheel test
above is the ground truth for this machine, but pin the torch/rocm version that works and re-test
on upgrades. VRAM discipline matters: a naive 64-DM-batch gather OOM'd 16 GiB; batch to fit.

## What stays CUDA-walled (do not spec)

**Heimdall / dedisp / astro-accelerate / FREDDA** (no HIP ports), **YOLO-CIANNA** (custom
C+CUDA framework), **turboSETI GPU mode** (CuPy; `amd-cupy` exists but is Instinct-primary and
untested on RDNA3 — run turboSETI on CPU), **FETCH** (TensorFlow; TF-ROCm is fragile — prefer a
PyTorch classifier). The ROCm-viable rule of thumb: **if it's pure PyTorch, it runs; if it ships
CUDA kernels, it doesn't.**

## GPU-enabled avenues, ranked (a parallel track — does not displace the catalogue-work Tier 1)

1. **`torch-fdmt`: a pure-PyTorch Fast-DM-Transform package** — no maintained torch/JAX FDMT
   exists (verified; only Julia + a Python reference). Pure tensor ops → ROCm out of the box; the
   naive gather benchmark above (24×) de-risks the performance claim ("faster than CPU
   baselines", not "beats FREDDA"). Output: JOSS/RNAAS software paper + it upgrades plan 28's
   dedispersion leg to GPU on this very hardware. Composes `jansky.transients` as the
   correctness oracle (CPU reference implementation already tested).
2. **FAST-FREX FRB-classifier benchmark** — the public labelled FAST dataset
   (doi:10.57760/sciencedb.15070; 600 bursts + 1,000 RFI negatives, built *for* ML challenges);
   train a plain-PyTorch ResNet/UNet on ROCm, report precision/recall vs the paper's
   RaSPDAM baselines. 16 GiB VRAM is ample. Avoid custom-CUDA-op models.
3. **OVRO-LWA solar burst detector (first ML detector on the new archive)** — Ultralytics YOLO
   (pure PyTorch) fine-tuned on OVRO-LWA public dynamic spectra (Apr 2024+, cycle-25 max);
   nothing published on this archive yet; ties directly into the e-Callisto census (a calibrated
   cross-check instrument). Risk: sparse labels — seed-annotate, or transfer from the SWSC-2026
   e-Callisto YOLOv5 work (code availability unverified).
4. **Radio-morphology fine-tuning (Zoobot / STRADAViT) on LoTSS DR3 / VLASS cutouts** — Zoobot is
   PyTorch and explicitly gaming-GPU-certified; classify the VLASS variability-census sources'
   morphology as the science hook (connects an existing slice to interpretation).
5. **U-Net RFI segmentation on Apertif TD DR2** (weak labels from AOFlagger; test flag-vs-FRB
   preservation on the 24 FRB-positive pointings) — viable but the most crowded of the five.

Interaction with the main shortlist: the Tier-1 catalogue picks (Stokes-V, LPT, SPICE-RACS,
CHIME wait-times) are **GPU-idle and unaffected** — the GPU track runs in parallel. The single
best synergy: **#1 (`torch-fdmt`) + plan 28 (BL single-pulse)** become one slice arc — build the
tool, validate against `jansky.transients` synthetics, run it on a real public filterbank
(BL / Apertif / Parkes PSRDA), publish the benchmark.
