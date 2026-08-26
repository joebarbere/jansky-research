# jansky-research

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21482378.svg)](https://doi.org/10.5281/zenodo.21482378)
[![GitHub Sponsors](https://img.shields.io/github/sponsors/joebarbere?logo=githubsponsors)](https://github.com/sponsors/joebarbere)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-donate-ff5f5f?logo=kofi&logoColor=white)](https://ko-fi.com/joebarbere)

Amateur radio-astronomy research, end to end. A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky teaches radio
astronomy, this repo does original amateur research on public data.

The unit of work is a "slice": one question, one tested tool, one real dataset, one write-up.
The tools are pure NumPy/SciPy/astropy, or pure PyTorch where device portability pays (the same
code runs on CPU and on an AMD GPU via ROCm). Each slice goes through an adversarial science
review before the paper is written, and a negative result gets the same write-up a positive one
would.

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

Forty-five slices, one line each; the long version of every row is `survey/<slice>-findings.md`:

| Slice | Tool | Outcome |
|-------|------|---------|
| FRB burst-statistics | `jansky_research.frbstats` | ✅ reproduced the CHIME repeater **width** result |
| Ultra-steep-spectrum hunt | `jansky_research.spectra` | ➖ a raw TGSS×NVSS USS cut is 17% pure / 14% complete (selection on the noisy index + matched-sensitivity truncation, both measured) |
| FRB repeater periodicity | `jansky_research.frbperiod` | ✅ recovered FRB 20180916B's **16.35-day** period (transit-collapsed Z²=23.9; coincidence probability ~10⁻², the honest headline) |
| SETI drift-search benchmark | `jansky_research.driftsearch` | ➖ benchmark built; the apparent "Voyager detection" was a **DC-spike artifact** caught in development; the tool reports a null |
| HI rotation curve | `jansky_research.hi` | ✅ recovered the **flat** (non-Keplerian) inner Milky Way curve |
| VLASS multi-epoch variability | `jansky_research.vlass` | ✅/➖ 703 deg² census: catalogue variability is **artifact-dominated**, but image-confirms **FK Comae Berenices** |
| Peaked-spectrum (GPS/CSS) selection | `jansky_research.peaked` | ✅ three-frequency curvature selector; **100% recovery** of a known HFP sample, high purity vs MHz-peaked |
| Southern peaked-spectrum (GLEAM-X×RACS) | `jansky_research.southern` | ✅ multi-band curvature that **measures** the turnover ν_pk; 90 candidates over a 3° cone, two systematic fixes |
| Radio–optical offsets (ICRF3 × Gaia × MOJAVE) | `jansky_research.offsets` | ✅ reproduces the AGN offset **excess tail** (≫ Rayleigh) **and its alignment with the parsec-scale jet** (KS p=3×10⁻²²) |
| Pulsar radio spectral indices (ATNF) | `jansky_research.pulsarspec` | ✅ steep mean spectrum (−1.77 raw, −1.88 completeness-limited); MSP−normal offsets outside [−0.22, +0.26] excluded |
| Sub-threshold radio stacking (SDSS quasars × VLASS-SE) | `jansky_research.stacking` | ✅ forced-photometry stacking with an off-source control; median flux of undetected quasars 40.8 µJy at 4.9σ |
| Multi-decade VLBI variability (Astrogeo) | `jansky_research.vlbi` | ✅ **control-floor** method recovers OJ 287 & BL Lac; 12–13 of 14 non-controls above the floor across its jackknife (unselected medians 0.305 vs 0.193) |
| Solar type III exciter speed (e-Callisto) | `jansky_research.solarbursts` | ✅ drift→beam speed **~0.14 c** on a clean isolated burst (Newkirk inversion) |
| Galactic Faraday rotation sky (Taylor+2009) | `jansky_research.rmsky` | ✅ plane enhancement **5.4 ± 0.6** (sky-block jackknife), inner-Galaxy sign antisymmetry at 2.8σ/5.0σ with alias-immune sign fractions |
| Pulsar P–Ṗ diagram (ATNF) | `jansky_research.ppdot` | ✅ MSP↔normal median fields 3.6 dex apart; 98.2% above the death line (92.5–99.6 across the valley); four named anchors gated |
| Inner-heliosphere type III (Wind/WAVES) | `jansky_research.windwaves` | ✅ beam tracked to ~10 R⊙ (the Alfvén surface), with peak-time and R² caveats |
| Interplanetary type III (STEREO/WAVES HFR) | `jansky_research.swaves` | ✅ beam tracked to 0.38 AU, with an R²-inflation caveat (few independent time samples) |
| 3D type III triangulation (STEREO-A+B DF) | `jansky_research.triangulate` | ✅ geometric vs plasma-frequency distance correlate at **r=0.989**; source localized in 3D |
| Euclidean source counts (NVSS) | `jansky_research.sourcecounts` | ✅ recovers the canonical **Hopkins 2003** 1.4 GHz counts; sub-Euclidean slope −1.91 |
| RACS Stokes-V coherent emitters (ASKAP via CASDA) | `jansky_research.stokesv` | ✅/➖ leakage-floor selection validated; forced photometry recovers I, but single-epoch V is variability-limited |
| Type III occurrence census vs the solar cycle (e-Callisto × SILSO) | `jansky_research.ecallisto_census` | ✅/➖ the coverage-corrected census statistic recovers an injected solar-cycle proportionality (r=0.97); the event stream is synthetic and the paper says so — the real multi-cycle ingest is future work |
| torch-fdmt: device-portable Fast DM Transform | `jansky_research.fdmt`+`singlepulse` | ✅ oracle-validated pure-PyTorch FDMT; real Crab DM recovered to 0.3% (giant pulse S/N 14); benchmarked fairly: FDMT-on-CPU beats brute-on-GPU 3.6× |
| RACS Stokes-V discovery: two-epoch photometry | `jansky_research.stokesv_discovery` | ✅/➖ recovers GJ 65 (BL+UV Cet) with a marginal 4.2-yr V decline — 10σ against image noise alone, 3.5–5σ once a per-epoch flux-scale term is included; all else quiescent (median 5σ limit 0.83 mJy); no new detections |
| LPT population catalogue (16 objects, provenance-typed) | `jansky_research.lpt` | ✅/➖ the death-line count is structural (margins 42–10⁶×; the claim-carrying 6/6 is the no-companion subset); the ~78-min boundary not significant under any labelling (exact p, power ≈ test size at the observed offset); caught a dropped digit in the review's own data file |
| RM structure functions (SPICE-RACS DR2) | `jansky_research.rmstructure` | ✅/➖ noise-debiased SF per \|b\| bin on the deduped DR2 goodRM sample; each plateau an upper bound; disc–halo ladder measured with per-bin errors |
| Jovian DAM occurrence census (Juno/Waves) | `jansky_research.junodam` | ✅/➖ the raw ~196× proximity trend is a threshold-amplified 1/r² sensitivity effect (corrected near/far 2.2×, an upper bound); Earth-canonical Io boxes do not organise orbital-vantage occurrence |
| Cassini SKR proximity census (RPWS KEY60S) | `jansky_research.skr` | ✅/➖ recovers the published dual rotation period (10.68/10.80 h, to 0.05%); the raw 3.3× proximity trend collapses to ~1.4× under the 1/r² null, still entangled with a 28° latitude difference |
| OVRO-LWA type II burst census (F13) | `jansky_research.typeii` | ✅/➖ streamed all 768 observing days 2024–2026 in memory → a null: the 332 candidates are false-positive dominated, so a blind spectral type II census fails in this RFI-heavy band |
| PTE-II per-source giant-pulse census (F10) | `jansky_research.pte2` | ✅/➖ uniform heavy-tail test over 363 Parkes pulsars → a null: the classification is detection-power limited, the tails are too steep for classic giant pulses, and the excess does not track Ė |
| JBO glitch waiting-time classification (F11) | `jansky_research.glitchpop` | ✅ gap-robust waiting-time classification of a pinned JBO snapshot (31 classifiable: 22 exponential, 8 quasi-periodic, 1 chance-consistent clustered); J0537, Vela and the out-of-sample Howitt anchor all recovered; no classification flip survives a deep bootstrap |
| Voyager 2 PRA ice-giant rotation periods (F9) | `jansky_research.vgpra` | ✅/➖ a controlled null: synthetic injections recover, neither real ice-giant period does — the PRA total power is not a clean rotational sinusoid, so the historical beaming modelling was essential |
| e-Callisto megaconstellation RFI trend (F17) | `jansky_research.rfitrend` | ✅/➖ a burst-immune line-excess metric over 2012–2026: the two line-sampling stations disagree in sign and HUMAIN's "rise" is a reconfiguration staircase, so no Starlink attribution; the null excludes only a common rise ≳15 log-units |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `jansky_research.rmdipole` | ✅/➖ the RM sky is isotropic at dipole order in its core; the significant dipole is carried entirely by the top-1% residual tail (clip: p 0.001→0.93), read as systematics |
| Uniform Cat-2 repeater timing census (CHIME/FRB) | `jansky_research.frbwait` | ✅/➖ anchor 20180916B re-found at 16.33 d and survives a transit-group null; the census median k is explained by transit censoring (mapped injection curve), so no population clustering claim; no new periods |
| Lensed-repeater search in Cat 2 (CHIME/FRB) | `jansky_research.frblens` | ✅/➖ 0/33 detections → lensed fraction < 0.37 at 95%, dividing by the summed injection efficiency (a source-count denominator would overstate the limit 4×) |
| torch-dsp: coherent dedispersion + RFI + FFA in pure PyTorch | `jansky_research.torchdsp` | ✅/➖ CHIME baseband burst re-dedispersed to its catalogue DM (conjugate-sign control 1.4 vs 4.0); SK/SumThreshold byte-identical to the CPU oracle, cross-device checked on ROCm; FFA 10.6× on GPU (single-session benchmark); the 2.1-s Crab file gives a period-search null with a measured 2σ/pulse injection limit |
| Radio survey of 56 WD-pulsar candidates (RACS+VLASS) | `jansky_research.wdpulsar` | ✅/➖ 0/51 candidates detected (I or V) to a median 3σ V limit of 0.41 mJy; AR Sco control re-found; J1912−4410 itself undetected |
| How often are long-period transients on? (VAST archival) | `jansky_research.lptduty` | ✅/➖ 3 LPTs caught in 1.1–4.5% of snapshots, 5 more limited to 3.2–10.7%; 2 of 10 constraints voided by a phase-sampling test; 7 discovery papers publish no reference epoch |
| Environment-split HI mass function (FASHI DR1) | `jansky_research.fashienv` | ✅/➖ void HIMF knee suppressed −0.26 dex (2.9σ) vs walls, reported as an upper bound pending a density-corrected estimator; group knee survivor-biased (stated) |
| SBI population inference for RACS Stokes-V emitters | `jansky_research.svsbi` | ✅/➖ with one detection no model parameter is measured: the LF break is a lower limit (log L* ≳ 13.9), f_beam is prior-located, only the slope median is stable |
| LPT catalogue v3 + Stokes-V forced photometry | `jansky_research.lptv` | ✅/➖ catalogue extended to 16 members; the ~78-min boundary still not significant (p=0.52); 1 secure + 1 candidate circular burst; uniform V-limit table (0.47 mJy median) |
| Inner Milky Way RC replication (plan 86) | `jansky_research.innerrc` | ✅ their ρ_DM arithmetic validates, but an interior refit of their own published curve gives 0.24 GeV/cm³ (1σ 0.16–0.31) — the published 0.107 is not uniquely determined by their curve |
| BL 3I/ATLAS GBT reproduction (plan 85) | `jansky_research.atlas3i` | ✅ the BL null reproduces from the public archive (1.12M raw hits → 261 survivors → 0 confirmed) at a matching 99.2 mW EIRP limit; the survivors are a taxonomy of two-position-filter evasion |
| DR20 BHM radio-counterpart census (plan 88) | `jansky_research.dr20radio` | ✅ first radio census of SDSS-V BHM quasars: RACS south 3.95%, VLASS north 4.67%; α measured by survival analysis (−0.755 ± 0.012, systematic range −0.55 to −0.90); the contrast is quoted as a ratio (1.47) |
| Type III synthesis: corona → 0.4 AU (4 instruments) | `jansky_research.type3synthesis` | ✅ unified drift-to-distance ladder; **geometric check on the model distance** (same-event r=0.989) |

Most of these are recover-a-known validations and method work. Several are negatives, written up
as negatives: the USS candidates dissolved into selection bias on a noisy index, the SETI "Voyager
detection" was a DC-spike artifact, the blind type II census is false-positive dominated, and the
megaconstellation RFI trend is systematics-limited. The review gate has caught at least one real
physics, citation, or statistics error in every slice before write-up, which is the argument for
keeping it.

## What's next

The idea backlog is [`fable-ideas.md`](fable-ideas.md) (a 2026-07 deep re-scan of the open-data
landscape). Every entry there now has a plan (`plans/38`–`93`) and most have been executed, so
current work is depth rather than breadth:

- **Review.** Every paper passes mechanical triage (`scripts/triage_papers.py`), and a deeper
  presenter/referee round is working through the list — the Reviewed column in the papers table
  below is the record. Every round so far has come back with revisions, which is the argument
  for finishing it before submitting anything.
- **Publishing.** The submission queue and venue reasoning are in the papers section below.
- **The station.** Once the rooftop receiver produces calibrated spectra, self-collected data
  joins the public-archive slices (plan 78; `fable-ideas.md` carries a station track for it).


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
adversarial presenter/referee round; the date links to the findings file with the full verdict.
Rounds during the 2026-08 house-style conversion also re-checked each converted paper's claims
against its committed evidence, which is where several of the fixes below came from.

| Paper | `papers/…` | Framing | Tests | Gate | Reviewed |
|-------|-----------|---------|-------|------|----------|
| FRB burst statistics, validated on CHIME/FRB Cat 1 | `frbstats/` | validation (tested, reproducible tool) | 97% | **major revision** ([referee 2026-08-25](survey/findings.md)) | ✅ 2026-08-25 — the source is now the unit of analysis (width D = 0.56 at p = 1e-5; DM honestly marginal at 0.044); Figure 1 plots the waits the Weibull was fitted to (multi-source fixture added); γ carries the joint bootstrap error with Cat 1's slope as the cited comparand (0.46σ); the macro writer is guarded and the silent synthetic fallback raises |
| Recovering FRB 20180916B's 16.35-day period | `frbperiod/` | validation | 100% | **major revision** ([referee 2026-08-26](survey/period-findings.md)) | ✅ 2026-08-26 — one epoch per transit (19→14; Z²=23.9, P=16.32±0.09 by bootstrap); the abstract now prints the measured coincidence probability (0.004 uniform / 0.015 clustered null — the argument's real strength) instead of the disclaimed FAP; the null source's sensitivity is stated (twin power 0.32) and its numbers are in the JSON; per-burst MJDs committed; figure clobber hole closed |
| The flat inner Milky Way rotation curve from LAB HI | `hi/` | validation | 99% | **major revision** ([referee 2026-08-24](survey/hi-findings.md)) | ✅ 2026-08-25 — the citation became a measurement: an edge-fit estimator on the same spectra gives 243±8 (2.9% above V₀) vs 257 threshold, a measured 13.8 km/s offset; flatness quantified (slope 2.8±1.8) |
| A SETI drift-search benchmark + Voyager-1 null | `driftsearch/` | benchmark + negative result | 100% | — | ✅ 2026-08-23 — one-detector claim de-generalised; "Honest" dropped from the title |
| A TGSS×NVSS USS cut is neither pure nor complete | `spectra/` | cautionary negative | 99% | **major revision** ([referee 2026-08-25](survey/uss-findings.md)) | ✅ 2026-08-25 — whole-field reference comparison committed (456 rows): population offset +0.002 ± 0.005 (no flux-scale effect; dec-edge refuted); candidate offset −0.11 vs −0.05 ± 0.03 from the committed selection-on-noisy-α model; cut scored 17% pure / 14% complete with the matched-sensitivity truncation quantified (α* ≈ −1.01; S₁₅₀ ≥ 102 mJy at −1.65); the limit-row disagreement flagged; both papers rewritten around the measured mechanisms |
| VLASS multi-epoch variability: a 703 deg² census + FK Com | `vlass/` | methodology + validation (recovers FK Com) | 92% | **major revision** ([referee 2026-08-24](survey/vlass-findings.md)) | ✅ 2026-08-24 — completeness re-derived with the fixed error model (99.5% at 10×, was "52% saturation"); the archival rejection is committed evidence (fixed-position photometry goes negative on the artifact) |
| Three-frequency curvature selection of peaked-spectrum sources | `peaked/` | methodology + two recover-a-known validations | 98% | **major revision** ([referee 2026-08-25](survey/peaked-findings.md)) | ✅ 2026-08-25 — HFP validation re-measured with real TGSS fluxes (53% rising, 32% falling; the 100% was a flux cut); field-local limit makes the count 1+1 with the 8/6/1/0 sweep committed; candidates published with errors + in-pipeline SIMBAD/NED vetting |
| Measuring the turnover: southern peaked sources from GLEAM-X + RACS | `southern/` | methodology + measured-turnover candidate list | 97% | **major revision** ([referee 2026-08-24](survey/southern-findings.md)) | ✅ 2026-08-25 — re-run with the catalogue's own errors + a significance-tested gate: 41 candidates (was 90), density ratio vs RadioSED II 1.21 (was ~2.5×, an error-model artifact), Callingham accuracy improved to 0.091 dex, full catalogue committed |
| The AGN radio–optical offset excess and its alignment with the jet (ICRF3 × Gaia × MOJAVE) | `offsets/` | reproduction (excess + jet alignment) + reproducible catalogue | 99% | **minor revision** ([referee 2026-08-25](survey/offsets-findings.md)) | ✅ 2026-08-25 — X rebuilt as a Mahalanobis distance (both catalogues' correlations fetched), so 24.0× now stands against an exact null; directional test on the sign-symmetric null (208/83, p=1.5e-13) with the upstream/disk component reported; error model measured (core scale 1.10, inflation sweep); the 3502-row catalogue ships |
| The steep radio spectra of pulsars from the ATNF catalogue | `pulsarspec/` | reproduction + MSP/normal comparison | 100% | **major revision** ([referee 2026-08-25](survey/pulsarspec-findings.md)) | ✅ 2026-08-25 — selection restated as coverage-dominated (764/2536 have S400) with the committed completeness cut (−1.77 raw → −1.88); the exclusion is now the CI [−0.22, +0.26] with 0.24 labelled 50%-power; per-source CSV, cut sweep, permutation p, and systematics floor committed |
| Forced-photometry stacking with an off-source control: SDSS quasars in VLASS-SE | `stacking/` | methodology + population median flux | 99% | **major revision** ([referee 2026-08-25](survey/stacking-findings.md)) | ✅ 2026-08-25 — every flux is now a forced central-pixel read (40.8 µJy at 4.9σ on the grown 281-source stack); the identity "calibration" is retracted in print and replaced by a jittered injection that can fail (0.923 ± 0.024); an off-source control stack exists (−17.5 ± 8.2 µJy, no pedestal); the detected-source cut is implemented (19 excluded); the faint bin is reported as the non-detection it is; "mean" → clipped median with a skewed fixture proving the difference |
| Multi-decade parsec-scale VLBI variability from Astrogeo | `vlbi/` | control-floor method + recover-a-known (OJ 287, BL Lac) | 100% | **major revision** ([referee 2026-08-25](survey/vlbi-findings.md)) | ✅ 2026-08-25 — Fassnacht & Taylor citation fixed; floor carries its jackknife (0.186–0.201, count 13→12; 10 at the control max; 12 epoch-matched) with the median-threshold selection function stated; amplitude comparison de-circularized (unselected 0.305 vs 0.193); CTA 102 resolved at 5 decimals (above by 7×10⁻⁵, z=0.0); per-source table + m_d committed |
| A solar type III exciter speed from an e-Callisto dynamic spectrum | `solarbursts/` | method + recover-a-known | 97% | **major revision** ([referee 2026-08-24](survey/solarbursts-findings.md)) | ✅ 2026-08-25 — fit converges (0.117 c, R²=0.935, spread 0.117–0.133 committed); figure draws the fit it captions; all four candidates + the storm control committed at one parameterization, and the gate is now falsifiable |
| The Galactic Faraday rotation sky from the Taylor et al. 2009 RM catalogue | `rmsky/` | reproduction (plane enhancement + sign antisymmetry) | 98% | **major revision** ([referee 2026-08-25](survey/rmsky-findings.md)) | ✅ 2026-08-25 — every error is now a 571-block sky jackknife: ratio 5.4 ± 0.6 (bootstrap said ±0.15), inner sign pattern at 2.8σ/5.0σ backed by alias-immune sign fractions; catalogue committed; Aitoff labels fixed; variants + correlated-fixture test added |
| The pulsar P–Ṗ diagram from the ATNF catalogue | `ppdot/` | reproduction + population classes | 100% | **minor revision** ([referee 2026-08-25](survey/ppdot-findings.md)) | ✅ 2026-08-25 — magnetar median reframed as threshold-set (sweep 12.7–14.1 committed; the robust claim is the 3.64-dex MSP↔normal separation); 2017 snapshot vintage disclosed; death-valley + flux-cut sweeps, four gated anchors, and the 2052-row CSV committed |
| Tracking an inner-heliosphere type III beam with Wind/WAVES | `windwaves/` | method + recover-a-known (to the Alfvén surface) | 98% | **major revision** ([referee 2026-08-25](survey/windwaves-findings.md)) | ✅ 2026-08-25 — grid committed (0.054–0.215 c, 6–19.5 R⊙) with the f²∝n degeneracy stated; headline rebuilt as a per-sample fit, 0.102 ± 0.011 c; ridge + burst epoch committed |
| Tracking a type III beam to 0.4 AU with STEREO/WAVES | `swaves/` | method + recover-a-known (genuinely interplanetary) | 100% | **major revision** ([referee 2026-08-25](survey/swaves-findings.md)) | ✅ 2026-08-25 — per-sample fit gives 0.160 ± 0.014 c (R²=0.986; the converged per-point clip collapsed to 0.031 c — measured, and retired); grid + ridge + epoch committed; title carries the 0.2–0.4 AU range |
| 3D triangulation of a type III source with STEREO-A+B direction-finding | `triangulate/` | method + independent geometric-vs-plasma distance cross-check | 97% | **major revision** ([referee 2026-08-25](survey/triangulate-findings.md)) | ✅ 2026-08-25 — longitude wrap-fixed (179.5° ± 7.6, was 168.9); the comparison reframed additive (13.4 ± 3.3 R⊙ constant, slope 1.12 — a ~4° DF bias, measured by the committed calibration grid; density enhancement retired); Leblanc reproduced to 3.4 R⊙ rms after the offset; epochs, directions, grids committed |
| Recovering the canonical 1.4 GHz Euclidean source counts from NVSS | `sourcecounts/` | reproduction (Hopkins 2003) | 100% | **major revision** ([referee 2026-08-25](survey/sourcecounts-findings.md)) | ✅ 2026-08-25 — fixed 0.2-dex grid: ratio 1.01 ± 0.02 with a closed budget (scatter 0.054 vs computed 0.052; χ²/dof 0.7 with the reference's own residual); cosmic variance moved to the normalisation (1.7–2.1%, computed); threshold profile, cut sweep, FoF sweep, per-bin table, and a closed-form estimator test committed |
| A type III drift-to-distance framework, corona to 0.4 AU | `type3synthesis/` | synthesis + same-event geometric check on the model distance | 99% | **major revision** ([referee 2026-08-25](survey/type3synthesis-findings.md)) | ✅ 2026-08-25 — re-run against the current siblings with a staleness-guard test; ladder carries the (harmonic × density) grids inline; the cross-check reframed additive (13.4 ± 3.3 R⊙, model reproduced to 3.4 R⊙ rms after one pointing bias; log-log shape NOT validated, slope 0.653); "geometrically validated" left the title; reach = the plasma band edge 0.384 AU |
| A reproducible RACS Stokes-V coherent-emitter pipeline (and single-epoch limits) | `stokesv/` | tooling + single-epoch variability limit | 100% | **major revision** ([referee 2026-08-25](survey/stokesv-findings.md)) | ✅ 2026-08-25 — I+V pinned to one obs_id per target: the impossible ratios become interpretable I non-detections (9 valid of 15); V in 4/9 (44%, Wilson 29–61%), threshold-robust; every row committed with product identifiers |
| A streaming e-Callisto burst-ingest pipeline with cross-station coincidence QC | `ecallisto_pipeline/` | automation pattern + coincidence-vetted burst events (rejects single-station RFI) | — | — | ✅ 2026-08-23 — abstract's only result was from the wrong run; macros namespaced, both legs' evidence split |
| A coverage-corrected type III occurrence census (method + recover-a-known) | `ecallisto_census/` | census statistic + three-arm validation (one arm can fail) toward a multi-cycle census | 100% | **major revision** ([referee 2026-08-26](survey/ecallisto-census-findings.md)) | ✅ 2026-08-26 — the real leg is an audited ingest failure (45 of 168 days ingested; 123 throttled with data listed, now recorded as missing, never zero); r=0.28 retired; the validation gained the arms that can fail (growth: raw 0.707→0.975; saturating confirmation: 0.604 with cov-corr −0.568, committed); linear arm labelled correct-by-construction; macros namespaced and code-written; k=0.03 stated |
| A pure-PyTorch Fast DM Transform (device-portable dedispersion) | `torchfdmt/` | tool + oracle validation + real Crab recover-a-known + CPU/GPU benchmark | 100% | **major revision** ([referee 2026-08-24](survey/torchfdmt-findings.md)) | ✅ 2026-08-25 — figure replotted with the statistic that produced the DM; 200-rep shift null committed (p=0.03 for the peak height; the positional p=0.0035 + boxcar S/N 14 carry the detection); benchmark re-measured in one invocation (24×, not the spliced 29×) |
| Two-epoch RACS Stokes-V forced photometry of nearby M dwarfs | `stokesv_discovery/` | method + GJ 65 variability recovery + upper-limit census | 100% | **major revision** ([referee 2026-08-12](survey/stokesv-discovery-findings.md)) | ✅ 2026-08-12 — 13 findings; 10σ demoted to marginal |
| A provenance-carrying LPT population catalogue | `lpt/` | cross-checked table + regenerable P–Ṗ statistics (novelty scoped vs the review's own diagram) | 100% | **major revision** ([referee 2026-08-26](survey/lpt-findings.md)) | ✅ 2026-08-26 — the paper now contains the catalogue (generated 16-row table, all discovery papers cited, DAS with DOI); figure regenerated from the corrected CSV with a derived title; the split test is truly exact (11,440 partitions, p=0.526 published labelling) with the labelling table primary (unknowns excluded: Δ=0.413, p=0.167) and measured power (≈ size at the observed offset); 9/9 reframed as structural (margins 42–10⁶×, valley sweep, 6/6 no-companion carries the claim); "verified"→"cross-checked" with the five-defect disclosure; the CSV is cell-pinned by test |
| LPT catalogue v3 + Stokes-V forced photometry | `lptv/` | 3 verified 2026 rows (N=16) + first systematic multi-epoch forced-V limit table at all LPT positions | 98% | **major revision** ([referee 2026-08-13](survey/lptv-findings.md)) | ✅ 2026-08-14 — re-measured forced; 3 detections, veto retracted |
| Cassini SKR proximity census | `skr/` | rotation-period anchor (0.05%) + 1/r² sensitivity-null bounding of the SKR occurrence-vs-range trend to a ~1.4× near-null | 96% | **major revision** ([referee 2026-08-24](survey/skr-findings.md)) | ✅ 2026-08-24 — 1.39 retracted: three null models give 1.4/1.7/3.6, so no residual is quoted; raw trend 3.33 ± 0.93 |
| OVRO-LWA type II burst census | `typeii/` | slow-drift detector + in-memory streaming of the whole 2024–2026 archive → a null (a blind spectral census is false-positive dominated) | 96% | **major revision** ([referee 2026-08-24](survey/typeii-findings.md)) | ✅ 2026-08-24 — the confusable contaminant injected: false-positive rate ~38% against slow-drift background; match-rate deficit retracted as a coverage artifact |
| e-Callisto megaconstellation RFI trend | `rfitrend/` | notch-robust narrowband-UEM-line metric over 2012–2026 + coherence verdict with its measured power curve → a systematics-limited null | 91% | **major revision** ([referee 2026-08-26](survey/rfitrend-findings.md)) | ✅ 2026-08-26 — HUMAIN reported as a staircase (steps +7.5/−12/+15 at 2017.2/2018.9/2022.2, no regime rises; the down-step spans deployment); prieto2020 correctly attributed (DOI-verified); the null quotes its sensitivity (90% power only at ~15 log-units); equal-n test shows the pooled slope is amplitude-dominated (0.261 [0.224, 0.292]); Theil–Sen CIs, window-mirror check, shape discriminant, flank sweep with the ALMATY-matching row, and the splice-free write path all committed |
| BL 3I/ATLAS GBT L-band reproduction | `atlas3i/` | independent pipeline over the public archive → the null reproduces (261 survivors, 0 confirmed) + the filter-evasion taxonomy (satellite downlinks defeat two-position filters by design) | 99% | **minor revision** ([referee 2026-08-11](survey/atlas3i-findings.md)) | ✅ 2026-08-11 — 7 findings fixed |
| DR20 BHM radio-counterpart census | `dr20radio/` | first radio census of the first southern SDSS quasar spectra (RACS 3.95%, VLASS north 4.67%) + a north/south contrast reported as an **α-dependent range** (0.23–1.66 pp), with measured chance rates, a measured RACS-footprint bound, and carton cross-frequency controls | 99% | **major revision** ([4 referee rounds](survey/dr20radio-findings.md)) | ✅ 2026-08-12 — 4 rounds; α measured by survival analysis, 4 systematics measured |
| Inner Milky Way rotation-curve replication | `innerrc/` | two-leg replication (their published tables + raw HI4PI) of Sofue & Kohno 2025's inner RC and low halo-only local dark-matter density | 98% | **major revision** ([referee 2026-08-12](survey/innerrc-findings.md)) | ✅ 2026-08-12 — 19 findings, both blockers fixed |
| Voyager 2 PRA ice-giant rotation periods | `vgpra/` | blind Lomb-Scargle of the PRA flux + rotation-block bootstrap → a controlled null (synthetic recovers an injected period; neither real Uranus/Neptune period is recovered — historical geometric modelling was essential) | 96% | **major revision** ([referee 2026-08-13](survey/vgpra-findings.md)) | ✅ 2026-08-13 — 17 findings; null restated with FAP + selection function |
| PTE-II per-source giant-pulse census | `pte2/` | floor-robust giant-pulse excess test over 363 Parkes pulsars + ATNF Ė cross-match → a null (detection-power-limited heavy-tail fraction; no Ė trend; tails too steep for classic giant pulses) | 97% | **major revision** ([referee 2026-08-13](survey/pte2-findings.md)) | ✅ 2026-08-13 — headline re-bracketed [0, 23]; test shown non-identifiable |
| JBO glitch waiting-time classification | `glitchpop/` | monitoring-gap-robust per-pulsar waiting-time classification of the live JBO catalogue (exp/quasi-periodic) + post-2018 change table; recover-a-known passes (J0537 & Vela quasi-periodic, gap excision required) | 95% | **major revision** ([referee 2026-08-24](survey/glitchpop-findings.md)) | ✅ 2026-08-24 — the title flip was MC noise (dissolves at 2×10⁵ boots); census pinned to a committed snapshot; "~184" was vintage arithmetic (measured: 89) |
| RM structure functions from SPICE-RACS DR1 | `rmstructure/` | method + recover-a-known + bounded high-\|b\| estimate | 96% | **major revision** ([referee 2026-08-24](survey/rmstructure-findings.md)) | ✅ 2026-08-24 — sample deduped to the release's own 246,508 count; ratio now 11.0 ± 1.1 (sky-block jackknife; the bootstrap said ±0.10); quality-flag claim measured on committed runs |
| Jovian DAM occurrence from Juno/Waves | `junodam/` | census method + proximity result + reduced Io-region contrast from orbit | 100% | **major revision** ([referee 2026-08-24](survey/junodam-findings.md)) | ✅ 2026-08-24 — test committed (CI 0.53–1.59: unity and 1.6× both admitted); Io-B alone enhanced (1.67); title softened |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `rmdipole/` | method + injection validation + isotropy null (tail-carried anisotropy disclosed as systematics) | 98% | **major revision** ([referee 2026-08-24](survey/rmdipole-findings.md)) | ✅ 2026-08-24 — headline restated as a limit (|p|/m > 0.35 excluded); injection unbiased over 20 seeds; clip shown non-circular |
| Uniform Cat-2 repeater timing census | `frbwait/` | anchor recovery + population k census + no-new-periods verdict | 98% | **major revision** ([referee 2026-08-24](survey/frbwait-findings.md)) | ✅ 2026-08-24 — population claim retracted: the censoring curve explains the median k; anchor survives the grouped null, the other two peaks do not |
| Lensed-repeater search in Cat 2 | `frblens/` | first empirical lensed-fraction limit + transit selection function + null-design lesson | 98% | **major revision** ([referee 2026-08-12](survey/frblens-findings.md)) | ✅ 2026-08-12 — 16 findings; limit corrected 4× weaker |
| torch-dsp: the coherent-DSP suite in pure PyTorch | `torchdsp/` | per-kernel oracle validation + real CHIME/Crab legs (CPU) + ROCm benchmarks | 93% | **major revision** ([referee 2026-08-25](survey/torchdsp-findings.md)) | ✅ 2026-08-25 — benchmark re-measured in one `--benchmark-only` session (FFA 10.6×, derived macros, hardware introspected); portability measured (cross-device kernel checks committed); chirp checked against the dispersion law + conjugate-sign trial (1.4 vs 4.0); Crab null given its noise level (5.0) and 2σ/pulse injection limit |
| A radio survey of the WD-pulsar candidates | `wdpulsar/` | AR Sco recover-a-known + systematic RACS/VLASS non-detection limit table | 94% | **major revision** ([referee 2026-08-13](survey/wdpulsar-findings.md)) | ✅ 2026-08-14 — re-measured forced; null real, bound stated |
| How often are long-period transients on? (VAST archival; RNAAS note) | `lptduty/` | sensitivity-weighted duty-cycle limits for the LPT class | 96% | [findings](survey/lptduty-findings.md) | — |
| The environment-split FASHI HI mass function | `fashienv/` | first env split of the FASHI HIMF + injection-validated 1/Vmax + ALFALFA-void-consistent (same SDSS voids, deeper independent HI — not an independent realisation; reported as an upper bound pending a density-corrected estimator) | 96% | **major revision** ([referee 2026-08-12](survey/fashienv-findings.md)) | ✅ 2026-08-12 — 20 findings; offset reframed as an upper bound |
| SBI for the RACS Stokes-V emitter population | `svsbi/` | first calibrated beaming-fraction posterior + SBC-validated coverage + ROCm-trained NPE | 93% | **major revision** ([referee 2026-08-12](survey/svsbi-findings.md)) | ✅ 2026-08-12 — 20 findings; log L* retracted to a lower limit |

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
- **Short notes:** three results are condensed to
  [Research Notes of the AAS](https://journals.aas.org/research-notes/) — the frbstats validation
  (`papers/frbstats/rnaas.tex`), the WD-pulsar survey (`papers/wdpulsar/rnaas.tex`, refereed and
  ready to submit), and the LPT duty-cycle constraint (`papers/lptduty/rnaas.tex`).
- **arXiv:** reserved for the genuine-novelty, real-data papers. The current queue, in order, is
  `atlas3i/`, `dr20radio/`, `lptv/`, and `innerrc/` — each has been through multiple referee
  rounds (see the Reviewed column below). Behind them: `frblens/` (the first catalogue-level
  lensed-repeater search) and `fashienv/`. Recover-a-knowns and method demos stay in the repo +
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
    rmstructure.py rmdipole.py frbwait.py frblens.py junodam.py torchdsp.py wdpulsar.py fashienv.py svsbi.py lptv.py skr.py typeii.py rfitrend.py vgpra.py pte2.py glitchpop.py
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
  plans/                 # numbered slice specs (00-93); the lasting record is each slice's
                         #   survey/*-findings.md + papers/<slice>/
  fable-ideas.md         # current plan-ready idea list (2026-07 deep re-scan; supersedes the
                         #   opportunity-scan shortlist)
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
