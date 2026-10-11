# jansky-research

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21482378.svg)](https://doi.org/10.5281/zenodo.21482378)
[![GitHub Sponsors](https://img.shields.io/github/sponsors/joebarbere?logo=githubsponsors)](https://github.com/sponsors/joebarbere)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-donate-ff5f5f?logo=kofi&logoColor=white)](https://ko-fi.com/joebarbere)

Amateur radio-astronomy research, end to end, on public data. A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky teaches radio
astronomy, this repo does original research with it.

**What it is for.** Most radio-astronomy data are public, and most published results are never
independently re-checked. This repo asks small, answerable questions of those archives (does a
published null reproduce, are a tool's error bars honest, is a claimed trend a survey artefact?)
and answers each with a tested tool and a short paper. The unit of work is a "slice": one
question, one tested tool, one real dataset, one write-up. The tools are pure NumPy/SciPy/astropy,
or pure PyTorch where device portability pays (the same code runs on CPU and on an AMD GPU).

**Why it is useful.**
- *Checkable:* every number in every paper is generated from committed results, and the whole
  chain regenerates from a clean checkout. 48 slices, each with its tool, tests and paper.
- *Honest about negatives:* nulls, inconclusive tests and retractions are written up the same
  way as positives; most of the results here are limits and cautions, and say so.
- *Reviewed adversarially:* each result passes a science-review gate and, before submission,
  presenter/referee rounds whose findings and retractions are kept in `survey/`.
- *Reusable:* the tools are tested modules (85% coverage floor) that run on a desktop, and
  the lessons about how analyses go wrong are recorded where the next slice will read them.

**Recent headlines.**
- *2026-10:* zdm's per-burst redshift intervals for fast radio bursts place host-identified
  bursts too far away: 95% coverage 0.81 against a floor of 0.86 on 36 bursts it was not fitted
  to (`dmzcal`, RNAAS note accepted after three referee rounds).
- *2026-10:* three more RNAAS notes refereed to acceptance: the Voyager 2 ice-giant rotation
  null (`vgpra`), the TGSS×NVSS ultra-steep-spectrum caution (`spectra`), and the independent
  reproduction of Breakthrough Listen's 3I/ATLAS null (`atlas3i`, also bound for arXiv).
- *2026-10:* a two-beam FASHI × ALFALFA test of H I beam blending, with controls frozen before
  the data, came out inconclusive and is published as such (`hiblend`).
- *2026-09:* the first blind, survey-scale radio proper-motion search across four VLASS epochs recovers
  UV/BL Ceti and finds no new movers (`vlasspm`).
- Release [v1.13.0](https://github.com/joebarbere/jansky-research/releases/tag/v1.13.0)
  (2026-10-10), archived on Zenodo.

**Up next** (from [`ideas.md`](ideas.md); each starts with a novelty and data-access check):
a Stokes-V census of Gaia ultracool dwarfs across two RACS epochs (plan 59); broadband
technosignature power limits from survey circular polarization (plan 53); an anomaly sweep of a
survey nobody has swept (plan 52); and, once the rooftop receiver is producing calibrated
spectra, the first self-collected data (plan 78). See [What's next](#whats-next).

**In this README:** [AI use disclosure](#ai-use-disclosure) · [Method](#method) ·
[Results](#results) · [What's next](#whats-next) · [Papers](#papers) ·
[Quickstart](#quickstart) · [The rooftop station](#the-rooftop-station) ·
[Relation to `jansky`](#relation-to-jansky) · [Layout](#layout)

## AI use disclosure

The tools, analyses, and paper drafts in this repo were developed collaboratively
with Anthropic's Claude, directed by me — an amateur, not a professional astronomer.
Each paper carries its own AI-use disclosure, and an AI/LLM is not an author.

The method section below describes real guardrails — every tool must recover a
known, published result before it's trusted with anything blind, an adversarial
science-review gate tries to knock each finding down, tests hold an 85% coverage
floor in CI, and every reported number regenerates from a clean checkout. But be
clear about what those guardrails are for: they make the work checkable, not
expert-reviewed. No professional radio astronomer has peer-reviewed these
results. That's why findings are framed as candidates, validations, and limits —
never discoveries — and why negatives are reported as negatives.

If you're qualified to review a slice, that's the most valuable contribution you
can make: each one is self-contained under `plans/`, `survey/`, and `papers/`, and
[issues](https://github.com/joebarbere/jansky-research/issues) citing a specific
result are very welcome.

## Method

Every result follows the same path, specified in advance in a `plans/NN-slug.md` file: a tested
tool (offline synthetic fixture, 85% coverage floor in CI), a run on real public data, an
adversarial science-review gate, then the write-up (`survey/<slice>-findings.md` and
`papers/<slice>/`). Findings are framed as candidates, validations, and limits, and every
reported number regenerates from the pipeline. Survey gaps that weren't picked are kept as a
backlog in `survey/candidate-gaps.md`.

## Results

Forty-eight slices, one line each; the long version of every row is `survey/<slice>-findings.md`:

| Slice | Tool | Outcome |
|-------|------|---------|
| FRB burst-statistics | `jansky_research.frbstats` | ✅ reproduced the CHIME repeater **width** result |
| Ultra-steep-spectrum hunt | `jansky_research.spectra` | ➖ a raw TGSS×NVSS USS cut is 17% pure / 14% complete against the reference's single-component USS sources (1/6, 1/7: disagreements at the threshold and sources outside the matched sample); no relative calibration offset; selection on the noisy index explains part of the candidates' steepening (17–81% by noise model; a fifth to a quarter with the robust core, leaving a ~2σ residual unexplained) |
| FRB repeater periodicity | `jansky_research.frbperiod` | ✅ recovered FRB 20180916B's **16.35-day** period (transit-collapsed Z²=23.9; coincidence probability ~10⁻², the honest headline) |
| SETI drift-search benchmark | `jansky_research.driftsearch` | ✅ benchmark + a real recovery: the carrier is located in the data (8419.297 MHz, −0.374 Hz/s) and recovered at S/N ≈ 997; the DC-spike caution stands (artifact reads 2×10⁵); an earlier draft's "Voyager null" withdrawn as a targeting error |
| HI rotation curve | `jansky_research.hi` | ✅ the inner Milky Way curve from 71 LAB sightlines, cross-validated against the VGPS: the edge fit agrees to one LAB channel (1.04 ± 0.21 km/s; the estimator's own choices span 2.1) with matched slopes (3.59 vs 3.57 km/s/kpc on common longitudes); non-Keplerian by ~80 km/s; the fitted slope depends on the inner radius cut, so the earlier "not flat" is a chord across an undulation |
| VLASS multi-epoch variability | `jansky_research.vlass` | ✅/➖ 703 deg² census: catalogue variability is **artifact-dominated**, but image-confirms **FK Comae Berenices** |
| Peaked-spectrum (GPS/CSS) selection | `jansky_research.peaked` | ✅/➖ three-frequency curvature selector with a field-local TGSS limit: 1 peaked + 1 GHz-peaked at 29.4 mJy, a count set by the limit (8 → 0 across 20–40 mJy); on the Dallacasa HFP candidates 53% rising, 32% falling — the earlier "100% recovery" was a flux cut the sample could not fail |
| Southern peaked-spectrum (GLEAM-X×RACS) | `jansky_research.southern` | ✅ multi-band curvature that **measures** the turnover ν_pk: 41 candidates over a 3° cone with the catalogue's own errors (49% of turnovers interpolated across the RACS gap); Callingham ν_pk to 0.091 dex; density 1.21× RadioSED II |
| Radio–optical offsets (ICRF3 × Gaia × MOJAVE) | `jansky_research.offsets` | ✅ reproduces the AGN offset **excess tail** (≫ Rayleigh) **and its alignment with the parsec-scale jet** (KS p=3×10⁻²²) |
| Pulsar radio spectral indices (ATNF) | `jansky_research.pulsarspec` | ✅ steep mean spectrum (−1.77 raw, −1.88 completeness-limited); MSP−normal offsets outside [−0.22, +0.26] excluded |
| Sub-threshold radio stacking (SDSS quasars × VLASS-SE) | `jansky_research.stacking` | ✅ forced-photometry stacking with an off-source control; median flux of undetected quasars 40.8 µJy at 4.9σ |
| Multi-decade VLBI variability (Astrogeo) | `jansky_research.vlbi` | ✅ **control-floor** method recovers OJ 287 & BL Lac; 12–13 of 14 non-controls above the floor across its jackknife (unselected medians 0.305 vs 0.193) |
| Solar type III exciter speed (e-Callisto) | `jansky_research.solarbursts` | ✅ drift→beam speed **0.117 c** (analysis spread 0.117–0.133 c; ~3× emission-mode/density bracket) on a clean isolated burst (Newkirk inversion) |
| Galactic Faraday rotation sky (Taylor+2009) | `jansky_research.rmsky` | ✅ plane enhancement **5.4 ± 0.6** (sky-block jackknife), inner-Galaxy sign antisymmetry at 2.8σ/5.0σ with alias-immune sign fractions |
| Pulsar P–Ṗ diagram (ATNF) | `jansky_research.ppdot` | ✅ MSP↔normal median fields 3.6 dex apart; 98.2% above the death line (92.5–99.6 across the valley); four named anchors gated |
| Inner-heliosphere type III (Wind/WAVES) | `jansky_research.windwaves` | ✅ beam speed 0.102 ± 0.011 c; across the emission-mode × density grid 0.054–0.215 c and a reach of 6–19.5 R⊙ set by the band edge, so no heliospheric regime is claimed; residual curvature unexplained |
| Interplanetary type III (STEREO/WAVES HFR) | `jansky_research.swaves` | ✅ per-sample fit gives a radial speed of 0.160 ± 0.014 c (harmonic, Leblanc); the reach of 0.19–0.38 AU is set by the HFR band edge and the emission mode, not the event; mode × density grid committed |
| 3D type III triangulation (STEREO-A+B DF) | `jansky_research.triangulate` | ✅ source localized in 3D (HEEQ lon 179.5° ± 7.6); triangulated distances sit a constant 13.4 ± 3.3 R⊙ outside Leblanc (a ~4° direction-finding bias, not a density enhancement), 3.4 R⊙ rms after the offset; the old r = 0.989 headline dropped (any monotone law scores > 0.95) |
| Euclidean source counts (NVSS) | `jansky_research.sourcecounts` | ✅ matches the Hopkins 2003 1.4 GHz counts (median ratio 1.01 ± 0.02, closed error budget); sub-Euclidean slope −1.90 ± 0.02 against −1.86 for the reference through the same estimator |
| RACS Stokes-V coherent emitters (ASKAP via CASDA) | `jansky_research.stokesv` | ✅/➖ selection validated on a synthetic field; I+V pinned to one obs_id: 9 of 15 targets valid (6 are single-epoch I non-detections); significant V in 4/9 (44%, Wilson 29–61%), threshold-robust; image I matches the catalogue only in the median |
| Type III occurrence census vs the solar cycle (e-Callisto × SILSO) | `jansky_research.ecallisto_census` | ✅/➖ the coverage correction rescues a broken raw correlation under network growth (r 0.71 → 0.98) and measurably fails under saturating confirmation (0.60); the linear arm is correct by construction; the real ingest is audited at 45 of 168 days |
| torch-fdmt: device-portable Fast DM Transform | `jansky_research.fdmt`+`singlepulse` | ✅ oracle-validated pure-PyTorch FDMT; real Crab DM recovered to 0.3% (giant pulse S/N 14); benchmarked fairly: FDMT-on-CPU beats brute-on-GPU 3.6× |
| RACS Stokes-V discovery: two-epoch photometry | `jansky_research.stokesv_discovery` | ✅/➖ recovers GJ 65 (BL+UV Cet); its V decline (10σ on image noise, 3.5–5σ with a flux-scale term) tracks I, so V/I is constant to 3.4%, a calibration signature; the others quiescent (median 5σ 0.83 mJy, beam-scale maxima from a peak search); no new detections |
| LPT population catalogue (16 objects, provenance-typed) | `jansky_research.lpt` | ✅/➖ the death-line count is structural (margins 42–10⁶×; the claim-carrying 6/6 is the no-companion subset); the ~78-min boundary not significant under any labelling (exact p, power 0.17 at the observed offset (0.06 at the published labelling)); caught a dropped digit in the review's own data file |
| RM structure functions (SPICE-RACS DR2) | `jansky_research.rmstructure` | ✅/➖ noise-debiased SF per \|b\| bin on the deduped DR2 goodRM sample; each plateau an upper bound; disc–halo ladder measured with per-bin errors |
| Jovian DAM occurrence census (Juno/Waves) | `jansky_research.junodam` | ✅/➖ the raw 330× proximity trend is threshold-amplified 1/r² sensitivity (censored near/far 2.7×, entangled with geometry); Io boxes are not detected to organise orbital occurrence (CI 0.53–1.59 admits both 1 and 1.6×); Io-B alone enhanced (1.67) |
| Cassini SKR proximity census (RPWS KEY60S) | `jansky_research.skr` | ✅/➖ the flux-series anchor matches the published 10.68/10.80 h to 0.05% (a certification, not a re-detection); raw near/far 3.3 ± 0.9 stands; the 1/r² decomposition is withdrawn as method-dominated (nulls 0.6–3.6), entangled with a 28° latitude span |
| OVRO-LWA type II burst census (F13) | `jansky_research.typeii` | ✅/➖ streamed all 768 observing days 2024–2026 in memory → a null: the 332 candidates are false-positive dominated, so a blind spectral type II census fails in this RFI-heavy band |
| PTE-II per-source giant-pulse census (F10) | `jansky_research.pte2` | ✅/➖ uniform heavy-tail test over 363 Parkes pulsars → a null: not identifiable without the detection floor (incidence 0 to ~23 of 136); detection-power limited; tails too steep for classic giant pulses; the Ė test is uninformative |
| JBO glitch waiting-time classification (F11) | `jansky_research.glitchpop` | ✅ gap-robust waiting-time classification of a pinned JBO snapshot (31 classifiable: 22 exponential, 8 quasi-periodic, 1 chance-consistent clustered); J0537, Vela and the out-of-sample Howitt anchor all recovered; no classification flip survives a deep bootstrap |
| Voyager 2 PRA ice-giant rotation periods (F9) | `jansky_research.vgpra` | ✅/➖ a controlled null: synthetic injections recover, neither real ice-giant period does; injections into the real series exclude a sinusoid at the historical period only at depths of 20% (Uranus) / 10% (Neptune) — the PRA total power is not a clean rotational sinusoid |
| e-Callisto megaconstellation RFI trend (F17) | `jansky_research.rfitrend` | ✅/➖ a burst-immune line-excess metric over 2012–2026: the two line-sampling stations disagree in sign and HUMAIN's "rise" is a reconfiguration staircase, so no Starlink attribution; the null excludes only a common rise ≳15 log-units |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `jansky_research.rmdipole` | ✅/➖ the RA-projected RM sky is consistent with isotropy once the top-1% residual tail is clipped (p 0.001 → 0.93): it excludes |p|/m > 0.35, nothing below, and is blind to a polar-axis dipole; the tail's dipole points at the declination edge, read as systematics |
| Uniform Cat-2 repeater timing census (CHIME/FRB) | `jansky_research.frbwait` | ✅/➖ anchor 20180916B re-found at 16.33 d and survives a transit-group null; the census median k is explained by transit censoring (mapped injection curve), so no population clustering claim; no new periods |
| Lensed-repeater search in Cat 2 (CHIME/FRB) | `jansky_research.frblens` | ✅/➖ 0/33 detections → lensed fraction < 0.37 at 95%, dividing by the summed injection efficiency (a source-count denominator would overstate the limit 4×) |
| torch-dsp: coherent dedispersion + RFI + FFA in pure PyTorch | `jansky_research.torchdsp` | ✅/➖ CHIME baseband burst re-dedispersed to its catalogue DM (conjugate-sign control 1.4 vs 4.0); SK and sequential-mode SumThreshold byte-identical to the CPU oracle, cross-device checked on ROCm; FFA 10.6× on GPU (single-session benchmark); the 2.1-s Crab file gives a period-search null with a measured 2σ/pulse injection limit |
| Radio survey of 56 WD-pulsar candidates (RACS+VLASS) | `jansky_research.wdpulsar` | ✅/➖ 0/49 RACS-covered candidates detected in I or V (median 3σ V limit 0.41 mJy) → persistently radio-loud fraction < 6.1% (95%); the AR Sco control re-detected in I in 2/5 epochs, never in V at 5σ; J1912−4410 undetected |
| How often are long-period transients on? (VAST archival) | `jansky_research.lptduty` | ✅/➖ 3 LPTs caught in 1.1–4.5% of snapshots (one common rate, LRT p=0.45); implied active fractions 0.04–0.74 — the class is not uniformly "mostly off"; 2 of 10 constraints void (+3 conditional on Ṗ); 6 discovery papers publish no reference epoch |
| Environment-split HI mass functions + two null tests (FASHI DR2) | `jansky_research.fashienv` | ➖ void − wall −0.088 dex, 1.1–2.5σ beyond a random-placement and a redshift × sensitivity label-shuffle null — not a detection; group − field +0.123 dex (3.5σ vs the shuffle) carries a FAST-beam blending signature (+0.16 dex more HI at fixed optical luminosity when a neighbour is inside the flux aperture, +0.03 just outside), not attributed; supersedes the DR1 −0.26 dex claim |
| SBI population inference for RACS Stokes-V emitters | `jansky_research.svsbi` | ✅/➖ with one detection no model parameter is measured: the LF break is a lower limit (log L* ≳ 13.9), f_beam is prior-located, only the slope median is stable |
| LPT catalogue v3 + Stokes-V forced photometry | `jansky_research.lptv` | ✅/➖ catalogue extended to 16; the ~78-min boundary is not significant (p = 0.52); RACS: 3 V detections, all of published epochs (zero new); VAST: 2 unpublished J1839−0756 pulses, one a candidate interpulse; uniform V-limit table (0.47 mJy median) |
| Inner Milky Way RC replication (plan 86) | `jansky_research.innerrc` | ✅ their ρ_DM arithmetic validates, but an interior refit of their own published curve gives 0.24 GeV/cm³ (1σ 0.16–0.31) — the published 0.107 is not uniquely determined by their curve |
| BL 3I/ATLAS GBT reproduction (plan 85) | `jansky_research.atlas3i` | ✅ the BL null reproduces across 1–12 GHz from the public archive (1.94M raw hits → 294 survivors → 0 confirmed) at a matching 99.2 mW L-band EIRP limit (up to ~3× weaker off the drift grid; completeness unmeasured); the survivors are a taxonomy of two-position-filter evasion |
| DR20 BHM radio-counterpart census (plan 88) | `jansky_research.dr20radio` | ✅ first radio census of SDSS-V BHM quasars: RACS south 3.95%, VLASS north 4.17–4.38% per epoch (4.67% in any epoch, raw); α −0.6 to −0.9 depending on which survey selects; north/south ratio 1.47 (1.39–1.61 from α), attributed to survey effects |
| Blind VLASS proper-motion search (plan 64) | `jansky_research.vlasspm` | ✅/➖ UV/BL Ceti recovered blind at 3.41 ± 0.04″/yr (three epochs, parallax-aware; system 3.34); zero new movers after calibrated errors and a compactness cut; 95% limit < 9.3e-5 per deg² on persistent compact optically dark movers (2.4–3.75 mJy, beyond 8 pc, 1.1–5″/yr) |
| Two-beam FASHI × ALFALFA beam-blending test (plans 97–98) | `jansky_research.hiblend` | ➖ inconclusive: 22,889 matched galaxies; the spectral null fails (β_null 0.26 ± 0.09) for an unidentified reason (the isolated-sample S/N misfit does not reproduce it; robust to 2–8° strip bootstraps); the injection test carries a +0.54 covariate-absorption bias plus a beam signal set by ALFALFA's unmeasured response to added flux, so its pass certifies nothing; a pre-specified recalibration stopped at a held-out gate that failed on χ² and on a tolerance a perfect calibration fails 2/3 of the time — no β is interpreted |
| Coverage audit of FRB DM→z posteriors on localized hosts (plan 99) | `jansky_research.dmzcal` | ➖/✅ negative: on 36 held-out spectroscopic hosts, zdm's p(z\|DM) fails a frozen coverage rule one-sidedly (95% coverage 0.81 vs floor 0.86, KS p = 1e-4; posterior medians ~1.85× the true z); not fixed by a constant host-DM change; the cause (population model, host selection, or z-dependent host DM) is not separable |
| Type III synthesis: corona → 0.4 AU (4 instruments) | `jansky_research.type3synthesis` | ✅/➖ unified drift-to-distance ladder to 0.384 AU (band edge) on three separate events; the same-event geometric check is additive (13.4 ± 3.3 R⊙ pointing offset, 3.4 R⊙ rms after it) and does not validate the model's log-log shape (slope 0.65) |

Most of these are recover-a-known validations and method work. Several are negatives, written up
as negatives: the USS cut proved 17% pure and 14% complete against the reference, the SETI "Voyager
detection" was a DC-spike artifact, the blind type II census is false-positive dominated, and the
megaconstellation RFI trend is systematics-limited. The review gate has caught at least one real
physics, citation, or statistics error in every slice before write-up, which is the argument for
keeping it.

## What's next

The idea backlog is [`ideas.md`](ideas.md) (a 2026-07 deep re-scan of the open-data
landscape, with dated updates since: 2026-09 releases and leads L1–L8; 2026-10 archival-mining
leads L9–L11, two demoted at GATE-0 and one blocked on data access). Every July entry has a plan
(`plans/38`–`93`) and most have been executed, so current work is depth rather than breadth:

- **Review.** Every paper passes mechanical triage (`scripts/triage_papers.py`), and a deeper
  presenter/referee round is working through the list — the Reviewed column in the papers table
  below is the record. Every round so far has come back with revisions, which is the argument
  for finishing it before submitting anything.
- **New slices.** Planned and waiting on their novelty and data-access checks: plan 59 (RACS
  Stokes-V census of Gaia ultracool dwarfs), plan 53 (broadband technosignature EIRP limits from
  survey Stokes V), plan 52 (self-supervised anomaly sweep of RACS or LoTSS DR3) and plan 51
  (an external FAST-FREX benchmark with an open-weights burst classifier). Plan 61 (Apertif
  single-pulse reprocessing) is re-scoped to a pulsar field and waits on an ASTRON account.
- **Publishing.** The submission queue and venue reasoning are in the papers section below.
- **The station.** Once the rooftop receiver produces calibrated spectra, self-collected data
  joins the public-archive slices (plan 78; `ideas.md` carries a station track for it).


### Running the full test suite locally

CI installs the CPU torch wheel (`uv sync --extra fdmt`), so the `fdmt`, `singlepulse` and
`torchdsp` tests all run there. A plain `uv sync` omits them: they `importorskip("torch")`
and the three modules report 0% coverage, which understates the suite by ~7 points
(89.9% without, **96.5%** with). To reproduce CI locally:

```bash
uv sync --extra fdmt   # once; ~200 MB CPU wheel from the pinned index
make cov
```

## Papers

Every slice is written up as an AASTeX paper under `papers/<slice>/`, authored by Joseph
Barbere with Claude credited in an AI-use disclosure and a `\software{}` citation (an AI/LLM is
not an eligible author). Each builds with containerized tectonic and takes its headline numbers
from a pipeline-written `generated/macros.tex`. Most are recover-a-known validations and method
papers; the framing column says which is which. The Reviewed column is the record of the
review campaign: a paper gets a mechanical triage pass (every paper passes it), then an
adversarial presenter/referee round; the Gate column links to the findings file with the full verdict.
The Tests column is the slice module's line coverage in CI (the 2026-10-11 run on `main`).
Rounds during the 2026-08 house-style conversion also re-checked each converted paper's claims
against its committed evidence, which is where several of the fixes below came from.

| Paper | `papers/…` | Framing | Tests | Gate | Reviewed |
|-------|-----------|---------|-------|------|----------|
| FRB burst statistics, validated on CHIME/FRB Cat 1 | `frbstats/` | validation (tested, reproducible tool) | 95% | **major revision** ([referee 2026-08-25](survey/findings.md)) | ✅ 2026-08-25 — the source is now the unit of analysis (width D = 0.56 at p = 1e-5; DM honestly marginal at 0.044); Figure 1 plots the waits the Weibull was fitted to (multi-source fixture added); γ carries the joint bootstrap error with Cat 1's slope as the cited comparand (0.46σ); the macro writer is guarded and the silent synthetic fallback raises |
| Recovering FRB 20180916B's 16.35-day period | `frbperiod/` | validation | 98% | **major revision** ([referee 2026-08-26](survey/period-findings.md)) | ✅ 2026-08-26 — one epoch per transit (19→14; Z²=23.9, P=16.32±0.09 by bootstrap); the abstract now prints the measured coincidence probability (0.004 uniform / 0.015 clustered null — the argument's real strength) instead of the disclaimed FAP; the null source's sensitivity is stated (twin power 0.32) and its numbers are in the JSON; per-burst MJDs committed; figure clobber hole closed |
| The inner Milky Way rotation curve from LAB HI, validated against the VGPS | `hi/` | validation | 94% | **major revision → resolved** ([referee 2026-08-31, RNAAS note](survey/hi-findings.md)) | ✅ 2026-08-31 — resampled 8 → 71 longitudes vs MG&D 2016's VGPS curve: agreement to about one LAB channel (+1.04, SEM 0.21; estimator choices span 2.12 km/s); slopes agree on matched sightlines (3.59 ± 0.76 vs 3.57 ± 0.71; the unmatched 2.52 vs 3.82 was a 2.2σ chord mismatch); the flatness claim **withdrawn** (1.79 ± 0.57 at N = 51, cut-dependent) and "Flat" dropped from the title; threshold bias +15.26 ± 0.56, not constant (sd 3.92) |
| A SETI drift-search benchmark + Voyager-1 recovered | `driftsearch/` | benchmark + real-data recovery + DC-spike caution | 94% | **major revision** ([referee 2026-08-26](survey/drift-findings.md)) | ✅ 2026-08-26 — the carrier is LOCATED in the data (never asserted): 8419.29696 MHz, −0.3741 Hz/s measured, recovered at S/N = 997 (committed, macro'd, regression-tested incl. a gated real-file test); the null withdrawn in print with the 0.92 MHz lesson kept as evidence; benchmark config committed; full p_detect matrix + noise-peak distribution + off-grid check; figure axis fixed |
| A TGSS×NVSS USS cut is neither pure nor complete | `spectra/` | cautionary negative | 97% | **RNAAS note accepted** ([referee rounds 2026-10-11](survey/uss-findings.md)) | ✅ 2026-10-11 — confirming rounds 2–4 on the note, then an end-to-end audit: selection explains 17–81% of the candidate offset depending on the index noise, a fifth to a quarter with the robust core (a ~2σ residual unexplained); 1/6 and 1/7 attributed to sources at the cut and three reference USS sources outside the matched sample; figure rebuilt with multi-component rows marked. ✅ 2026-08-25 — whole-field reference comparison committed (456 rows): population offset +0.002 ± 0.005 (no relative calibration offset; dec-edge refuted); candidate offset −0.11 vs −0.05 ± 0.03 from the committed selection-on-noisy-α model; cut scored 17% pure / 14% complete with the matched-sensitivity truncation quantified (α* ≈ −1.01; S₁₅₀ ≥ 102 mJy at −1.65); the limit-row disagreement flagged; both papers rewritten around the measured mechanisms |
| VLASS multi-epoch variability: a 703 deg² census + FK Com | `vlass/` | methodology + validation (recovers FK Com) | 92% | **major revision** ([referee 2026-08-24](survey/vlass-findings.md)) | ✅ 2026-08-24 — completeness re-derived with the fixed error model (99.5% at 10×, was "52% saturation"); the archival rejection is committed evidence (fixed-position photometry goes negative on the artifact) |
| Three-frequency curvature selection of peaked-spectrum sources | `peaked/` | methodology + two recover-a-known validations | 89% | **major revision** ([referee 2026-08-25](survey/peaked-findings.md)) | ✅ 2026-08-25 — HFP validation re-measured with real TGSS fluxes (53% rising, 32% falling; the 100% was a flux cut); field-local limit makes the count 1+1 with the 8/6/1/0 sweep committed; candidates published with errors + in-pipeline SIMBAD/NED vetting |
| Measuring the turnover: southern peaked sources from GLEAM-X + RACS | `southern/` | methodology + measured-turnover candidate list | 94% | **major revision** ([referee 2026-08-24](survey/southern-findings.md)) | ✅ 2026-08-25 — re-run with the catalogue's own errors + a significance-tested gate: 41 candidates (was 90), density ratio vs RadioSED II 1.21 (was ~2.5×, an error-model artifact), Callingham accuracy improved to 0.091 dex, full catalogue committed |
| The AGN radio–optical offset excess and its alignment with the jet (ICRF3 × Gaia × MOJAVE) | `offsets/` | reproduction (excess + jet alignment) + reproducible catalogue | 99% | **minor revision** ([referee 2026-08-25](survey/offsets-findings.md)) | ✅ 2026-08-25 — X rebuilt as a Mahalanobis distance (both catalogues' correlations fetched), so 24.0× now stands against an exact null; directional test on the sign-symmetric null (208/83, p=1.5e-13) with the upstream/disk component reported; error model measured (core scale 1.10, inflation sweep); the 3502-row catalogue ships |
| The steep radio spectra of pulsars from the ATNF catalogue | `pulsarspec/` | reproduction + MSP/normal comparison | 99% | **major revision** ([referee 2026-08-25](survey/pulsarspec-findings.md)) | ✅ 2026-08-25 — selection restated as coverage-dominated (764/2536 have S400) with the committed completeness cut (−1.77 raw → −1.88); the exclusion is now the CI [−0.22, +0.26] with 0.24 labelled 50%-power; per-source CSV, cut sweep, permutation p, and systematics floor committed |
| Forced-photometry stacking with an off-source control: SDSS quasars in VLASS-SE | `stacking/` | methodology + population median flux | 97% | **major revision** ([referee 2026-08-25](survey/stacking-findings.md)) | ✅ 2026-08-25 — every flux is now a forced central-pixel read (40.8 µJy at 4.9σ on the grown 281-source stack); the identity "calibration" is withdrawn in print and replaced by a jittered injection that can fail (0.923 ± 0.024); an off-source control stack exists (−17.5 ± 8.2 µJy, no pedestal); the detected-source cut is implemented (19 excluded); the faint bin is reported as the non-detection it is; "mean" → clipped median with a skewed fixture proving the difference |
| Multi-decade parsec-scale VLBI variability from Astrogeo | `vlbi/` | control-floor method + recover-a-known (OJ 287, BL Lac) | 99% | **major revision** ([referee 2026-08-25](survey/vlbi-findings.md)) | ✅ 2026-08-25 — Fassnacht & Taylor citation fixed; floor carries its jackknife (0.186–0.201, count 13→12; 10 at the control max; 12 epoch-matched) with the median-threshold selection function stated; amplitude comparison de-circularized (unselected 0.305 vs 0.193); CTA 102 resolved at 5 decimals (above by 7×10⁻⁵, z=0.0); per-source table + m_d committed |
| A solar type III exciter speed from an e-Callisto dynamic spectrum | `solarbursts/` | method + recover-a-known | 98% | **major revision** ([referee 2026-08-24](survey/solarbursts-findings.md)) | ✅ 2026-08-25 — fit converges (0.117 c, R²=0.935, spread 0.117–0.133 committed); figure draws the fit it captions; all four candidates + the storm control committed at one parameterization, and the gate is now falsifiable |
| The Galactic Faraday rotation sky from the Taylor et al. 2009 RM catalogue | `rmsky/` | reproduction (plane enhancement + sign antisymmetry) | 98% | **major revision** ([referee 2026-08-25](survey/rmsky-findings.md)) | ✅ 2026-08-25 — every error is now a 571-block sky jackknife: ratio 5.4 ± 0.6 (bootstrap said ±0.15), inner sign pattern at 2.8σ/5.0σ backed by alias-immune sign fractions; catalogue committed; Aitoff labels fixed; variants + correlated-fixture test added |
| The pulsar P–Ṗ diagram from the ATNF catalogue | `ppdot/` | reproduction + population classes | 99% | **minor revision** ([referee 2026-08-25](survey/ppdot-findings.md)) | ✅ 2026-08-25 — magnetar median reframed as threshold-set (sweep 12.7–14.1 committed; the robust claim is the 3.64-dex MSP↔normal separation); 2017 snapshot vintage disclosed; death-valley + flux-cut sweeps, four gated anchors, and the 2052-row CSV committed |
| Tracking an inner-heliosphere type III beam with Wind/WAVES | `windwaves/` | method + recover-a-known (6–19.5 R⊙ bracket; no regime claimed) | 92% | **major revision** ([referee 2026-08-25](survey/windwaves-findings.md)) | ✅ 2026-08-25 — grid committed (0.054–0.215 c, 6–19.5 R⊙) with the f²∝n degeneracy stated; headline rebuilt as a per-sample fit, 0.102 ± 0.011 c; ridge + burst epoch committed |
| Tracking a type III beam to 0.2–0.4 AU with STEREO/WAVES | `swaves/` | method + recover-a-known (genuinely interplanetary) | 100% | **major revision** ([referee 2026-08-25](survey/swaves-findings.md)) | ✅ 2026-08-25 — per-sample fit gives 0.160 ± 0.014 c (R²=0.986; the converged per-point clip collapsed to 0.031 c — measured, and retired); grid + ridge + epoch committed; title carries the 0.2–0.4 AU range |
| 3D triangulation of a type III source with STEREO-A+B direction-finding | `triangulate/` | method + independent geometric-vs-plasma distance cross-check | 97% | **major revision** ([referee 2026-08-25](survey/triangulate-findings.md)) | ✅ 2026-08-25 — longitude wrap-fixed (179.5° ± 7.6, was 168.9); the comparison reframed additive (13.4 ± 3.3 R⊙ constant, slope 1.12 — a ~4° DF bias, measured by the committed calibration grid; density enhancement retired); Leblanc reproduced to 3.4 R⊙ rms after the offset; epochs, directions, grids committed |
| Recovering the canonical 1.4 GHz Euclidean source counts from NVSS | `sourcecounts/` | reproduction (Hopkins 2003) | 99% | **major revision** ([referee 2026-08-25](survey/sourcecounts-findings.md)) | ✅ 2026-08-25 — fixed 0.2-dex grid: ratio 1.01 ± 0.02 with a closed budget (scatter 0.054 vs computed 0.052; χ²/dof 0.7 with the reference's own residual); cosmic variance moved to the normalisation (1.7–2.1%, computed); threshold profile, cut sweep, FoF sweep, per-bin table, and a closed-form estimator test committed |
| A type III drift-to-distance framework, corona to 0.4 AU | `type3synthesis/` | synthesis + same-event geometric check on the model distance | 99% | **major revision** ([referee 2026-08-25](survey/type3synthesis-findings.md)) | ✅ 2026-08-25 — re-run against the current siblings with a staleness-guard test; ladder carries the (harmonic × density) grids inline; the cross-check reframed additive (13.4 ± 3.3 R⊙, model reproduced to 3.4 R⊙ rms after one pointing bias; log-log shape NOT validated, slope 0.653); "geometrically validated" left the title; reach = the plasma band edge 0.384 AU |
| A reproducible RACS Stokes-V coherent-emitter pipeline, and the reach of single-epoch circular polarization | `stokesv/` | tooling + single-observation circular-polarization occurrence rate (4/9 valid targets) | 99% | **major revision** ([referee 2026-08-25](survey/stokesv-findings.md)) | ✅ 2026-08-25 — I+V pinned to one obs_id per target: the impossible ratios become interpretable I non-detections (9 valid of 15); V in 4/9 (44%, Wilson 29–61%), threshold-robust; every row committed with product identifiers |
| A streaming e-Callisto burst-ingest pipeline with cross-station coincidence QC | `ecallisto_pipeline/` | automation pattern + coincidence-vetted burst events (with its chance-coincidence rate) | — | **major revision** ([referee 2026-08-26](survey/ecallisto-findings.md)) | ✅ 2026-08-26 — fetch-by-filename + de-duplicated listings (the real day: 763 listed → 761 scanned, 2 failures counted, 2 candidates — BIR's is the known 11:50 UT burst, auditable from the committed t_peak_s); per-leg figures; DAG and CLI share one scan_file; the fixture gained contaminants that can defeat it (60/60-seed ensemble committed); event spans + closed-form chance rate quoted; the real-day null promoted to the abstract with its honest reading (13 of 15 stations saw nothing: the threshold, not coverage, binds) |
| A coverage-corrected type III occurrence census (method + recover-a-known) | `ecallisto_census/` | census statistic + three-arm validation (one arm can fail) toward a multi-cycle census | 100% | **major revision** ([referee 2026-08-26](survey/ecallisto-census-findings.md)) | ✅ 2026-08-26 — the real leg is an audited ingest failure (45 of 168 days ingested; 123 throttled with data listed, now recorded as missing, never zero); r=0.28 retired; the validation gained the arms that can fail (growth: raw 0.707→0.975; saturating confirmation: 0.604 with cov-corr −0.568, committed); linear arm labelled correct-by-construction; macros namespaced and code-written; k=0.03 stated |
| A pure-PyTorch Fast DM Transform (device-portable dedispersion) | `torchfdmt/` | tool + oracle validation + real Crab recover-a-known + CPU/GPU benchmark | 100% | **major revision** ([referee 2026-08-24](survey/torchfdmt-findings.md)) | ✅ 2026-08-25 — figure replotted with the statistic that produced the DM; 200-rep shift null committed (p=0.03 for the peak height; the positional p=0.0035 + boxcar S/N 14 carry the detection); benchmark re-measured in one invocation (24×, not the spliced 29×) |
| Two-epoch RACS-mid Stokes-V photometry of the nearest M dwarfs | `stokesv_discovery/` | method + GJ 65 recovery (decline consistent with a flux-scale offset) + peak-search upper-limit census | 100% | **major revision** ([referee 2026-08-12](survey/stokesv-discovery-findings.md)) | ✅ 2026-08-12 — 13 findings; the 10σ V decline now read as a flux-scale signature (V/I constant to 3.4%), and the quiescent limits as beam-scale peak-search maxima |
| A provenance-carrying LPT population catalogue | `lpt/` | cross-checked table + regenerable P–Ṗ statistics (novelty scoped vs the review's own diagram) | 95% | **major revision** ([referee 2026-08-26](survey/lpt-findings.md)) | ✅ 2026-08-26 — the paper now contains the catalogue (generated 16-row table, all discovery papers cited, DAS with DOI); figure regenerated from the corrected CSV with a derived title; the split test is truly exact (11,440 partitions, p=0.526 published labelling) with the labelling table primary (unknowns excluded: Δ=0.413, p=0.167) and measured power (0.17 at the primary offset; 0.06 ≈ size at the published one); 9/9 reframed as structural (margins 42–10⁶×, valley sweep, 6/6 no-companion carries the claim); "verified"→"cross-checked" with the five-defect disclosure; the CSV is cell-pinned by test |
| LPT catalogue v3 + Stokes-V forced photometry | `lptv/` | v3 catalogue (N = 16) + uniform forced-V limits at all LPT positions; RACS: 0 new detections (3 published epochs re-measured); VAST: 2 unreported J1839−0756 pulses, one a candidate interpulse | 96% | **major revision** ([referee round 3, 2026-08-15](survey/lptv-findings.md)) | ✅ 2026-08-15 — 3 rounds; all RACS detections are published epochs; interpulse downgraded to a candidate on polarization (18 ± 3% vs ≲10%); 6.4σ local → 5.4σ global; the earlier peak-search veto withdrawn |
| Cassini SKR proximity census | `skr/` | rotation-period anchor (0.05%) + occurrence-vs-range census: raw near/far trend 3.33 ± 0.93, the 1/r² decomposition shown method-limited (three nulls give 1.4/1.7/3.6) | 95% | **major revision** ([referee 2026-08-24](survey/skr-findings.md)) | ✅ 2026-08-24 — 1.39 withdrawn: three null models give 1.4/1.7/3.6, so no residual is quoted; raw trend 3.33 ± 0.93 |
| OVRO-LWA type II burst census | `typeii/` | slow-drift detector + in-memory streaming of the whole 2024–2026 archive → a null (a blind spectral census is false-positive dominated) | 96% | **major revision** ([referee 2026-08-24](survey/typeii-findings.md)) | ✅ 2026-08-24 — the confusable contaminant injected: false-positive rate ~38% against slow-drift background; match-rate deficit withdrawn as a coverage artifact |
| e-Callisto megaconstellation RFI trend | `rfitrend/` | notch-robust narrowband-UEM-line metric over 2012–2026 + coherence verdict with its measured power curve → a systematics-limited null | 93% | **major revision** ([referee 2026-08-26](survey/rfitrend-findings.md)) | ✅ 2026-08-26 — HUMAIN reported as a staircase (steps +7.5/−12/+15 at 2017.2/2018.9/2022.2, no regime rises; the down-step spans deployment); prieto2020 correctly attributed (DOI-verified); the null quotes its sensitivity (90% power only at ~15 log-units); equal-n test shows the pooled slope is amplitude-dominated (0.261 [0.224, 0.292]); Theil–Sen CIs, window-mirror check, shape discriminant, flank sweep with the ALMATY-matching row, and the splice-free write path all committed |
| BL 3I/ATLAS GBT 1–12 GHz reproduction | `atlas3i/` | independent pipeline over the public archive → the null reproduces (294 survivors, 0 confirmed) + the filter-evasion taxonomy (satellite downlinks defeat two-position filters by design) | 99% | **RNAAS note accepted; arXiv paper aligned** ([referee rounds 2026-10-11](survey/atlas3i-findings.md)) | ✅ 2026-10-11 — note gained main.tex's caveats (unmeasured completeness, band-dependent thresholds, ~3× off-grid penalty, no S/C/X sensitivity); cost figures labelled reconstructions. ✅ 2026-08-11 — 7 findings fixed |
| DR20 BHM radio-counterpart census | `dr20radio/` | first radio census of the first southern SDSS quasar spectra (RACS 3.95%, VLASS north 4.17/4.38% per epoch, chance-corrected) + a north/south ratio of 1.47 (1.39–1.61 from α) attributed to survey effects, with measured chance rates, a measured RACS-footprint bound, and carton cross-frequency controls | 98% | **major revision** ([4 referee rounds](survey/dr20radio-findings.md)) | ✅ 2026-08-12 — 4 rounds; α depends on which survey selects the sample (−0.6 to −0.9), so no single value is adopted; 4 systematics measured |
| An independent replication of the Sofue & Kohno inner Milky Way rotation curve, and the local dark-matter-density degeneracy | `innerrc/` | two-leg replication (their published tables + raw HI4PI) of Sofue & Kohno 2025's inner RC: their ρ_DM arithmetic reproduces, but their own curve prefers 0.24 GeV/cm³ (1σ 0.16–0.31, consensus inside) | 86% | **major revision** ([referee 2026-08-12](survey/innerrc-findings.md)) | ✅ 2026-08-12 — 19 findings, both blockers fixed |
| Voyager 2 PRA ice-giant rotation periods | `vgpra/` | blind Lomb-Scargle of the PRA flux + rotation-block bootstrap → a controlled null (synthetic recovers an injected period; neither real Uranus/Neptune period is recovered; injections exclude a sinusoid at the historical period only at depths of 20% / 10%) | 85% | **RNAAS note accepted** ([referee rounds 2026-10-11](survey/vgpra-findings.md)) | ✅ 2026-10-11 — four confirming rounds: the Aug fixes had never reached the note; rewritten with disjoint sub-bands, post-hoc controls and two post-result criterion changes disclosed. ✅ 2026-08-13 — 17 findings; null restated with FAP + selection function |
| PTE-II per-source heavy-tail test and its detection-power floor | `pte2/` | uniform per-source heavy-tail test over 363 Parkes pulsars → a non-identifiable null (giant-pulse incidence bracketed 0–~23 of 136; Ė test uninformative; tails too steep for classic giant pulses) | 78% | **major revision** ([referee 2026-08-13](survey/pte2-findings.md)) | ✅ 2026-08-13 — headline re-bracketed [0, 23]; test shown non-identifiable |
| JBO glitch waiting-time classification | `glitchpop/` | monitoring-gap-robust per-pulsar waiting-time classification of the live JBO catalogue (exp/quasi-periodic) + post-2018 increment (89 glitches; no classification flips); recover-a-known passes (J0537 & Vela quasi-periodic, gap excision required) | 90% | **major revision** ([referee 2026-08-24](survey/glitchpop-findings.md)) | ✅ 2026-08-24 — the title flip was MC noise (dissolves at 2×10⁵ boots); census pinned to a committed snapshot; "~184" was vintage arithmetic (measured: 89) |
| RM structure functions from SPICE-RACS DR1/DR2 | `rmstructure/` | method + error-model recover-a-known (jackknife reproduces field scatter) + definition-dependent DR2 disc–halo ratio + upper-bound plateaus | 96% | **major revision** ([referee 2026-08-26](survey/rmstructure-findings.md)) | ✅ 2026-08-26 — the quality-flag break claim withdrawn: the matched-pair comparison it lacked shows the shift preferred in only 37.5% of shared bootstrap replicates (both curves + intervals committed); every ladder error now a sky-block jackknife (bootstrap hid a factor ~7); block-size sweep committed with the fixture arbitrating 10° (1.1 calibrated, 1.7 ceiling); definitional variants run on DR2 (alt bins move the ratio +54%); CI asserts on the ensemble, not the outlier seed; per-leg figures in Galactic coordinates with the ladder panel |
| Jovian DAM occurrence from Juno/Waves | `junodam/` | census method + range trend shown to be inverse-square visibility + sensitivity-limited Io-box non-detection (Io-B alone enhanced) | 92% | **major revision** ([referee 2026-08-24](survey/junodam-findings.md)) | ✅ 2026-08-24 — test committed (CI 0.53–1.59: unity and 1.6× both admitted); Io-B alone enhanced (1.67); title softened |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `rmdipole/` | method + injection validation + isotropy null (tail-carried anisotropy disclosed as systematics) | 99% | **major revision** ([referee 2026-08-24](survey/rmdipole-findings.md)) | ✅ 2026-08-24 — headline restated as a limit (|p|/m > 0.35 excluded); injection unbiased over 20 seeds; clip shown non-circular |
| Uniform Cat-2 repeater timing census | `frbwait/` | anchor recovery + population k census + no-new-periods verdict | 94% | **major revision** ([referee 2026-08-24](survey/frbwait-findings.md)) | ✅ 2026-08-24 — population claim withdrawn: the censoring curve explains the median k; anchor survives the grouped null, the other two peaks do not |
| Lensed-repeater search in Cat 2 | `frblens/` | first empirical lensed-fraction limit + transit selection function + null-design lesson | 92% | **major revision** ([referee 2026-08-12](survey/frblens-findings.md)) | ✅ 2026-08-12 — 16 findings; limit corrected 4× weaker |
| torch-dsp: the coherent-DSP suite in pure PyTorch | `torchdsp/` | per-kernel oracle validation + real CHIME/Crab legs (CPU) + ROCm benchmarks | 92% | **major revision** ([referee 2026-08-25](survey/torchdsp-findings.md)) | ✅ 2026-08-25 — benchmark re-measured in one `--benchmark-only` session (FFA 10.6×, derived macros, hardware introspected); portability measured (cross-device kernel checks committed); chirp checked against the dispersion law + conjugate-sign trial (1.4 vs 4.0); Crab null given its noise level (5.0) and 2σ/pulse injection limit |
| A radio survey of the WD-pulsar candidates | `wdpulsar/` | AR Sco recover-a-known + systematic RACS/VLASS non-detection limit table | 96% | **RNAAS note refereed: minor, applied** ([referee 2026-08-15](survey/wdpulsar-findings.md)) | ✅ 2026-08-15 — note's inflated control efficiency restored (2/5 epochs, tied to f < 6.1%); main.tex peak-search survivors corrected. ✅ 2026-08-14 — re-measured forced; null real, bound stated |
| Implied active fractions of the long-period transients from archival ASKAP snapshots (RNAAS note) | `lptduty/` | per-snapshot detection probabilities + implied active fractions for the LPT class | 91% | **major revision** ([referee 2026-08-26](survey/lptduty-findings.md)) | ✅ 2026-08-26 — the weighting's true effect stated (0.02% at observed fluxes); counts fixed at the source (6 no-epoch, 2 published PEPOCHs); Kuiper p implemented + honest ×20 Bonferroni (verdict survives, 0.023); Ṗ folded into GATE-0 (8 stand, 5 unconditionally); the note now leads with implied f_active = p/((w+T)/P) (0.04–0.74; best limit <0.09 = GLEAM-X; one vacuous); \ldFActive carries its 95% interval and PEPOCH conditionality; discovery ephemerides cited; the macro producer is wired and lptduty builds via make |
| Environment-split HIMFs in FASHI DR2 against two null tests (RNAAS note) | `fashienv/` | random-placement + label-shuffle nulls, environment-restricted Vmax, counterpart-matched blending classes; void not detected, group offset shows a blending signature (supersedes the DR1 −0.26 dex paper) | 91% | minor revision at round 8, applied, not yet confirmed ([findings](survey/fashienv-findings.md)); package built | ✅ 2026-09-27 — 8 rounds (rounds 6–7 major, round 8 minor, all applied); group sigma dropped as occupancy-biased; blending favoured via HI at fixed optical luminosity |
| SBI for the RACS Stokes-V emitter population | `svsbi/` | first SBI framework for the M-dwarf Stokes-V population (beaming fraction barely updated from the prior; log L* a lower limit; SBC passes, with its blindness stated) + ROCm-trained NPE | 92% | **major revision** ([referee 2026-08-12](survey/svsbi-findings.md)) | ✅ 2026-08-12 — 20 findings; log L* withdrawn to a lower limit |
| No optically dark compact radio movers in four VLASS epochs (RNAAS note) | `vlasspm/` | first blind (Gaia-free) survey-scale radio proper-motion search: UV/BL Ceti recovered; calibrated error model; compactness cut calibrated on Gaia-star size noise (11 → 1 candidate); parallax-aware completeness; limit < 9.3e-5 deg⁻² | 99% | **ready** ([referee round 3](survey/vlasspm-findings.md)); package built | ✅ 2026-09-28 — 3 rounds; final accept with nits, applied |
| An inconclusive two-beam FASHI–ALFALFA test for H I beam blending (RNAAS note) | `hiblend/` | two-beam, same-galaxy test of H I beam blending with controls frozen before the data; negative control fails for an unidentified cause; the second attempt's gate failed on χ² and on a tolerance it could rarely pass | 98% | revised after referee round 4 ([findings](survey/hiblend-findings.md)) | ✅ 2026-10-10 — 4 rounds, all minor; calibration-limits claim withdrawn, overcorrections narrowed, C1 split into bias + unmeasured response |
| Dispersion-measure redshift intervals place host-identified FRBs too far away (RNAAS note) | `dmzcal/` | empirical coverage test of per-burst FRB p(z\|DM) intervals on held-out host-identified bursts (none found in a one-pass GATE-0 search), controls frozen before the data; zdm fails one-sidedly; cause left three-way (population model, host selection, or host DM) | 100% | accepted at referee round 3 ([findings](survey/dmzcal-findings.md)) | ✅ 2026-10-10 — 3 rounds (major → minor → accept); unfitted CHIME/MeerKAT models and omitted DM_G uncertainty disclosed, Galactic-σ and host ×1.5/×2 re-runs added, cause left three-way |

`make paper` builds every slice's PDF; `make papers-zip` bundles them all into one archive (the same
job runs in CI: the **`release` workflow** compiles every paper with tectonic and, on a `v*` tag,
attaches `jansky-research-papers-<tag>.zip` to a GitHub Release — and uploads it as a workflow
artifact on every manual run). `make arxiv` runs the bundled **`arxiv-submit` skill**
(`.claude/skills/arxiv-submit/`) to assemble and validate an upload package per paper
(`papers/<slice>/arxiv-submission/`: the LaTeX-source tarball with its `.bbl`, plus a `metadata.yaml`
capturing every arXiv submission property and a `CHECKLIST.md`). Orchestration is split by
cadence: static slices build through a server-less Snakemake file-DAG (`workflow/Snakefile`, run
by `make figures`, into `build/figures`, not the repo root), and the frequently-updated
e-Callisto archive is ingested by an Apache Airflow pipeline on rootless Podman (`airflow/`).

### Where to publish (and where not to)

Most of these papers are reproductions and negatives, so the venue is matched to the actual
contribution — the tooling and the reproducibility, not a novelty claim:

- **Software / citable archive:** the toolkit is meant for [JOSS](https://joss.theoj.org) (see
  `joss/paper.md`) and a [Zenodo](https://zenodo.org) DOI on release (`.zenodo.json`, `CITATION.cff`).
- **Short notes:** twelve results are condensed to
  [Research Notes of the AAS](https://journals.aas.org/research-notes/) — the frbstats validation
  (`papers/frbstats/rnaas.tex`), the WD-pulsar survey (`papers/wdpulsar/rnaas.tex`, refereed and
  ready to submit), the LPT duty-cycle constraint (`papers/lptduty/rnaas.tex`), the Voyager 2 PRA
  rotation-period null (`papers/vgpra/rnaas.tex`) and the TGSS×NVSS USS caution
  (`papers/spectra/rnaas.tex`), both accepted after confirming referee rounds on 2026-10-11, the
  3I/ATLAS reproduction (`papers/atlas3i/rnaas.tex`, accepted 2026-10-11; also on the arXiv
  track below), and the
  drift-search benchmark + Voyager-1 replication (`papers/driftsearch/rnaas.tex`, plan 94,
  refereed and ready to submit), and the LAB/VGPS terminal-velocity comparison
  (`papers/hi/rnaas.tex`, GATE-0 discharged and refereed 2026-08-31). Four are RNAAS-length as
  their `main.tex`: the FASHI DR2 environment-split caution (`papers/fashienv/`, eight referee
  rounds; round 8's minor fixes applied, not yet confirmed), the blind VLASS proper-motion limit (`papers/vlasspm/`, three rounds), the inconclusive
  two-beam blending test (`papers/hiblend/`, four referee rounds), and the FRB DM→z coverage
  audit (`papers/dmzcal/`, accepted at referee round 3).
- **arXiv:** reserved for the genuine-novelty, real-data papers. The current queue, in order, is
  `atlas3i/` (to go out alongside its RNAAS note; the two were aligned claim by claim),
  then `dr20radio/`, `lptv/`, and `innerrc/`, each of which still carries a major-revision
  verdict with no confirming round, so each needs one before submission (see the Gate column
  below). Behind them: `frblens/` (the first catalogue-level
  lensed-repeater search). Recover-a-knowns and method demos stay in the repo +
  Zenodo; software-pattern papers (`frbstats/`, `torchdsp/`, `torchfdmt/`, `ecallisto_pipeline/`)
  are JOSS candidates, not science preprints. The pure reproductions are not posted as a preprint
  batch — arXiv moderation expects a contribution, and "I reproduced a known result" belongs in
  the repo. The account-bound submission walkthrough lives in Joe's personal notes.

## Quickstart

New here? [`docs/usage.md`](docs/usage.md) is a short install-and-run guide,
[`docs/faq.md`](docs/faq.md) answers common questions about using the toolkit and the
papers in this repo, and [`CONTRIBUTING.md`](CONTRIBUTING.md) covers how to contribute,
report issues, and get support.

```bash
git clone https://github.com/joebarbere/jansky-research.git && cd jansky-research
uv sync                                   # env + jansky (installed from its pinned tag)
make test                                 # unit tests (offline, on synthetic fixtures)
make cov                                  # tests + 85% coverage floor

# Run a slice on real public data (each writes results/ + figures + macros into papers/<slice>/):
uv run python -m jansky_research.pipeline     # FRB burst statistics (CHIME catalog)
uv run python -m jansky_research.frbperiod    # FRB repeater periodicity
uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3   # USS hunt
uv run python -m jansky_research.driftsearch  # SETI injection-recovery benchmark
uv run python -m jansky_research.atlas3i --sweep  # BL 3I/ATLAS L-band reproduction (needs --extra voyager; ~13 h, 60 GB peak)
uv run python -m jansky_research.hi           # Milky Way HI rotation curve
uv run python -m jansky_research.vlass --ra 190 --dec 20 --radius 15  # VLASS variability census (needs --extra vlass)
uv run python -m jansky_research.singlepulse --benchmark --out .      # torch-fdmt benchmark (CPU;
                                          #   add --device cuda from a ROCm venv for the GPU column)
# (append --offline to run any slice on its synthetic fixture, no network)

# The papers + orchestration:
make figures                              # build every static slice via the Snakemake DAG (offline; needs --extra workflow)
make paper                                # tectonic -> all papers/<slice>/main.pdf (in a container)
make papers-zip                           # bundle every paper PDF into dist/jansky-research-papers-<TAG>.zip
make arxiv                                # assemble + validate an arXiv package per paper
make reproduce                            # fetch -> figures -> papers -> arXiv packages, end to end

# Streaming ingest (the e-Callisto archive) runs on Airflow + Podman:
make airflow-up COMPOSE="uvx podman-compose" && make dag-test DATE=2011-09-14
make ecallisto-day DATE=20110914          # the same day's scan WITHOUT Airflow (the shared worker)
```

See `REPRODUCING.md` for the full reproduction, the orchestration notes (Snakemake for the
static slices, Airflow for the streaming ingest), and offline mode.

## The rooftop station

Beyond the public archives, the [`station/`](station/) directory documents SDR-based instruments
built and operated from a Philadelphia rooftop (self-collected data, in progress) — a
[hydrogen-line receiver](station/hydrogen-line-receiver.md), a
[meteor-scatter station](station/meteor-scatter-station.md), and a planned
[two-dish interferometer](station/interferometry.md), plus [test-equipment](station/test-equipment.md)
and [long-duration operations](station/operations.md) notes. These are the build guides for the
instrument meant to feed self-collected data into future slices; the owner's working notes
(purchase log, prices, per-part rationale) live in an Obsidian vault, not this repo.

The station's **control software** is the sibling [`jansky-observe`](https://github.com/joebarbere/jansky-observe)
repo — now feature-complete across every planned milestone (capture + live view, an HI-line
classifier, calibration epochs, an unattended transit scheduler, drift-scan campaigns, a codified
observation-export bundle, and az/el rotator control, behind a read-mostly MCP surface). What's
left there is first light. Once it produces real spectra, the
[`pull-station-data`](.claude/skills/pull-station-data/SKILL.md) skill in this repo pulls its
codified observation bundles (averaged spectra + full provenance) into `data/station/` — exactly
the input format the hydrogen-line pipeline ([plan 78](plans/78-station-hline-pipeline.md))
consumes.

## Relation to `jansky`

This repo **depends on `jansky` as a library** and reuses its tested helpers (`jansky.transients`,
`jansky.rfi`, `jansky.timing`, `jansky.seti`, `jansky.sourcecounts`, `jansky.formats`,
`jansky.data`, …) rather than reimplementing them. It mirrors jansky's conventions: `uv`-managed,
ruff + mypy + pytest with an 85% coverage floor, Podman containers, and a `plans/NN-slug.md`
workflow. The `jansky` dependency installs from the pinned git tag `jansky@v0.2.0`, so a bare
`git clone && uv sync` works with no second checkout; for cross-repo development,
`eval "$(make -s dev-env)"` puts a sibling `../jansky/src` on `PYTHONPATH` ahead of the pinned tag.
See `pyproject.toml`.

## Layout

```
jansky-research/
  src/jansky_research/   # the tooling package (tested-helper pattern, 85% floor) — one module per slice
    data.py              # dataset registry + offline synthetic fallback
    frbstats.py spectra.py frbperiod.py driftsearch.py atlas3i.py hi.py vlass.py peaked.py southern.py
    offsets.py pulsarspec.py stacking.py vlbi.py solarbursts.py rmsky.py ppdot.py
    windwaves.py swaves.py triangulate.py sourcecounts.py type3synthesis.py
    ecallisto_catalog.py ecallisto_census.py stokesv.py stokesv_discovery.py lpt.py
    rmstructure.py rmdipole.py frbwait.py frblens.py junodam.py torchdsp.py wdpulsar.py fashienv.py svsbi.py lptv.py skr.py typeii.py rfitrend.py vgpra.py pte2.py glitchpop.py vlasspm.py hiblend.py dmzcal.py
    fdmt.py singlepulse.py  # torch-fdmt: pure PyTorch, device-portable (CPU or AMD GPU via ROCm)
    pipeline.py          # the FRB pipeline (shared by Make / notebook / Snakemake)
    report.py            # figure/macro emitters -> paper inputs
  survey/                # PERMANENT: literature.md, github-landscape.md, gap-analysis.md,
                         #   candidate-gaps.md + *-scan.md (backlog), and each slice's *-findings.md
  workflow/              # Snakefile: the server-less file-DAG that builds the static slices' inputs
  airflow/               # Airflow-on-Podman stack + the streaming e-Callisto ingest DAG
  papers/<slice>/        # one AASTeX paper per slice (main.tex + refs.bib tracked;
                         #   figures/, generated/, arxiv-submission/ are produced by make)
                         #   frbstats/ also has rnaas.tex (a Research Note of the AAS)
  station/               # build guides for the physical rooftop station (self-collected data, WIP)
  joss/                  # JOSS software paper (paper.md + paper.bib)
  CITATION.cff           # "Cite this repository"; .zenodo.json drives Zenodo archival
  containers/            # tectonic paper-build image
  .claude/skills/        # arxiv-submit, research-publish, traditional-style, casda-cutout-fetch,
                         #   radio-cutout, find-radio-papers, radio-source-lookup, idea-scan,
                         #   pull-station-data
  .claude/agents/        # science-reviewer, paper-presenter, paper-referee, style-editor,
                         #   dataset-analyst, pipeline-runner, results-interpreter, archive-scout,
                         #   radio-research-assistant
  plans/                 # numbered slice specs (00-99); the lasting record is each slice's
                         #   survey/*-findings.md + papers/<slice>/
  ideas.md         # current plan-ready idea list (2026-07 deep re-scan + dated updates;
                         #   supersedes the opportunity-scan shortlist)
```

## Support

The station track (plans 77–84) runs on rooftop hardware bought out of pocket. You can
help fund the buildout — rotator, second dish, coherent receiver — via
[GitHub Sponsors](https://github.com/sponsors/joebarbere) or
[Ko-fi](https://ko-fi.com/joebarbere). Honestly: that list is a general direction, not a
promise — the research plans and the observing priorities can and will change as results
come in.

## Contributing

Contributions, issues, and questions are welcome — see [`CONTRIBUTING.md`](CONTRIBUTING.md)
for how to set up, contribute a change, report a bug, and get support, and the
[Code of Conduct](CODE_OF_CONDUCT.md). Releases follow SemVer per
[`VERSIONING.md`](VERSIONING.md).

## License

Dual-licensed: the **code** (the `jansky_research` package and everything outside `papers/`) is
**MIT** — see [LICENSE](LICENSE); the **papers** in [`papers/`](papers/) are **CC BY 4.0** — see
[`papers/LICENSE`](papers/LICENSE).
