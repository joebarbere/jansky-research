# jansky-research

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21482378.svg)](https://doi.org/10.5281/zenodo.21482378)
[![GitHub Sponsors](https://img.shields.io/github/sponsors/joebarbere?logo=githubsponsors)](https://github.com/sponsors/joebarbere)
[![Ko-fi](https://img.shields.io/badge/Ko--fi-donate-ff5f5f?logo=kofi&logoColor=white)](https://ko-fi.com/joebarbere)

**Amateur radio-astronomy research, end to end.** A sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course: where jansky *teaches* radio
astronomy, this repo *does* original amateur research — surveys the landscape, builds tested
tooling, analyzes public data, automates it, and writes it up — with honesty as the first rule.

It began as a single vertical slice (survey → gap → tool → automated pipeline → reproducible paper)
and has grown into a **deep-research survey plus a set of self-contained research slices**, each:
a tested tool — pure NumPy/SciPy/astropy, or pure PyTorch where device portability pays (the same
code runs on CPU and on the workstation's AMD GPU via ROCm) — run on real public data, put through
an adversarial science-review gate, and written up, wins *and* negatives reported plainly.

**In this README:** [Method](#method) · [Results](#results) · [What's next](#whats-next) ·
[Papers](#papers) · [Quickstart](#quickstart) · [The rooftop station](#the-rooftop-station) ·
[Relation to `jansky`](#relation-to-jansky) · [Layout](#layout)

## Method

Every result follows the same discipline, recorded as a `plans/NN-slug.md` spec: build a tested
tool (offline synthetic fixture, 85% coverage floor) → run it on real public data → an adversarial
**science-reviewer gate** (a pass that has caught real blockers every time) → honest write-up
(`survey/<slice>-findings.md` + `papers/<slice>/`). Findings are framed as *candidates /
validations / limits*, never dressed-up discoveries; every reported number is reproducible. The
non-chosen survey gaps are preserved as a backlog in `survey/candidate-gaps.md`.

## Results

Twenty-six slices plus a synthesis, honestly tallied:

| Slice | Tool | Outcome |
|-------|------|---------|
| FRB burst-statistics | `jansky_research.frbstats` | ✅ reproduced the CHIME repeater **width** result |
| Ultra-steep-spectrum hunt | `jansky_research.spectra` | ➖ USS candidates **failed** the de Gasperin cross-check (negative) |
| FRB repeater periodicity | `jansky_research.frbperiod` | ✅ recovered FRB 20180916B's **16.35-day** period |
| SETI drift-search benchmark | `jansky_research.driftsearch` | ➖ benchmark built; the apparent "Voyager detection" was a **DC-spike artifact** caught in development — the tool reports an honest null |
| HI rotation curve | `jansky_research.hi` | ✅ recovered the **flat** (non-Keplerian) inner Milky Way curve |
| VLASS multi-epoch variability | `jansky_research.vlass` | ✅/➖ 703 deg² census: catalogue variability is **artifact-dominated**, but image-confirms **FK Comae Berenices** |
| Peaked-spectrum (GPS/CSS) selection | `jansky_research.peaked` | ✅ three-frequency curvature selector; **100% recovery** of a known HFP sample, high purity vs MHz-peaked |
| Southern peaked-spectrum (GLEAM-X×RACS) | `jansky_research.southern` | ✅ multi-band curvature that **measures** the turnover ν_pk; 90 candidates over a 3° cone, two systematic fixes |
| Radio–optical offsets (ICRF3 × Gaia × MOJAVE) | `jansky_research.offsets` | ✅ reproduces the AGN offset **excess tail** (≫ Rayleigh) **and its alignment with the parsec-scale jet** (KS p=3×10⁻²²) |
| Pulsar radio spectral indices (ATNF) | `jansky_research.pulsarspec` | ✅ reproduces the **steep** mean pulsar spectrum (α≈−1.8); MSPs ≈ normal pulsars |
| Sub-threshold radio stacking (SDSS quasars × VLASS-SE) | `jansky_research.stacking` | ✅ image-plane stacking with **injection-recovery** bias calibration; mean flux of undetected quasars |
| Multi-decade VLBI variability (Astrogeo) | `jansky_research.vlbi` | ✅ **control-floor** method recovers OJ 287 & BL Lac; blazars ~1.7× more variable than steady CSO controls |
| Solar type III exciter speed (e-Callisto) | `jansky_research.solarbursts` | ✅ drift→beam speed **~0.14 c** on a clean isolated burst (Newkirk inversion) |
| Galactic Faraday rotation sky (Taylor+2009) | `jansky_research.rmsky` | ✅ plane enhancement ratio **5.4**, inner-Galaxy RM **sign antisymmetry** recovered |
| Pulsar P–Ṗ diagram (ATNF) | `jansky_research.ppdot` | ✅ three classes over ~5 orders in B, 98% above the death line; Crab validates |
| Inner-heliosphere type III (Wind/WAVES) | `jansky_research.windwaves` | ✅ beam tracked to **~10 R⊙** (Alfvén surface); honest peak-time/R² caveats |
| Interplanetary type III (STEREO/WAVES HFR) | `jansky_research.swaves` | ✅ beam tracked to **0.38 AU**; honest R²-inflation (few independent time samples) caveat |
| 3D type III triangulation (STEREO-A+B DF) | `jansky_research.triangulate` | ✅ geometric vs plasma-frequency distance correlate at **r=0.989**; source localized in 3D |
| Euclidean source counts (NVSS) | `jansky_research.sourcecounts` | ✅ recovers the canonical **Hopkins 2003** 1.4 GHz counts; sub-Euclidean slope −1.91 |
| RACS Stokes-V coherent emitters (ASKAP via CASDA) | `jansky_research.stokesv` | ✅/➖ leakage-floor selection validated; forced photometry **recovers I**, but single-epoch **V is variability-limited** (honest) |
| Type III occurrence census vs the solar cycle (e-Callisto × SILSO) | `jansky_research.ecallisto_census` | ✅/➖ coverage-corrected census statistic **recovers an injected solar-cycle proportionality** (r=0.97); SILSO is real, event stream is synthetic — real multi-cycle ingest is future work (honest) |
| torch-fdmt: device-portable Fast DM Transform | `jansky_research.fdmt`+`singlepulse` | ✅ oracle-validated pure-PyTorch FDMT; **real Crab DM recovered to 0.3%** (giant pulse S/N 14); honest benchmark: FDMT-on-CPU beats brute-on-GPU 3.6× |
| RACS Stokes-V discovery: two-epoch forced photometry | `jansky_research.stokesv_discovery` | ✅/➖ recovers **GJ 65 (BL+UV Cet)** with a **10σ 4.2-yr V change**; all else quiescent (median 5σ limit 0.83 mJy); no new detections survived the novelty bar (honest) |
| LPT population catalogue (16 objects, provenance-typed) | `jansky_research.lpt` | ✅/➖ 9/9 Ṗ-constrained members below the death line; the hinted ~78-min binary boundary **not significant at N=16** (v3; honest, with demonstrated test power); caught a typo in the review's own data file |
| RM structure functions (SPICE-RACS DR1) | `jansky_research.rmstructure` | ✅/➖ noise-debiased SF per \|b\| bin; high-\|b\| plateau is an **upper bound** (intrinsic-scatter-dominated); disc–halo contrast awaits the public DR2 file |
| Jovian DAM occurrence census (Juno/Waves) | `jansky_research.junodam` | ✅/➖ 7-month census: the raw ~196× range-quartile proximity trend is **a threshold-amplified 1/r² sensitivity effect — the sensitivity-corrected near/far is only 2.2×** (null model added 2026-07, mirroring `skr`); Earth-canonical Io boxes do NOT coherently organise orbital-vantage occurrence (per-month median 0.87); GATE-2 caught the Io-phase convention blocker (disclosed) |
| Cassini SKR proximity census (RPWS KEY60S) | `jansky_research.skr` | ✅/➖ ports the junodam census to Saturn: recovers the published dual rotation period **10.68+10.80 h to 0.05%** (validation); the raw 3.3× proximity duty-cycle trend is **a bounded near-null** — the 1/r² sensitivity null collapses it to ~1.4×, still entangled with a 28° latitude difference (honest; not a proximity law) |
| OVRO-LWA type II burst census (F13) | `jansky_research.typeii` | ✅/➖ slow-drift+harmonic detector (synthetic SNR-completeness curve 1.0→0.33); streamed **all 765 observing days 2024–2026** from AWS Open Data in memory (no disk, 0 failures) → **an honest NULL**: the 331 candidates are **false-positive dominated** (matched-CME median 478≈background 379 km/s; drift⊥CME-speed r=0.09; 83% window-saturated) — a blind spectral type II census fails in this RFI-heavy band (why the archive detector is type-III-only) |
| PTE-II per-source giant-pulse census (F10) | `jansky_research.pte2` | ✅/➖ uniform per-source heavy-tail (giant-pulse) test over all 363 Parkes PTE-II pulsars (136 fitted); synthetic recover-a-known works (completeness 0.33→1.0 vs pulse count, FP ~0.05) → **an honest NULL**: 19% flagged heavy-tailed but the classification is **detection-power limited** (`count_limited`; heavy sources 2.4× more pulses, flag rate 0.09→0.29 by count), the tails are steep (index ~12, not the ~2–3 of true giant pulses), and the excess **does not correlate with Ė** (Spearman −0.03) — no giant-pulse-vs-Ė signal recoverable from S/N-only archival data |
| JBO glitch waiting-time classification (F11) | `jansky_research.glitchpop` | ✅ a monitoring-gap-robust per-pulsar waiting-time classification of the live Jodrell Bank catalogue (727 glitches/223 pulsars; 33 with ≥5 glitches → 23 exponential, 10 quasi-periodic) + the post-2018 change table (2 flips, 3 newly classifiable vs the end-2018 Basu subset); **recover-a-known passes** — J0537−6910 & Vela recovered as quasi-periodic, which *required* excising monitoring gaps (a 2264-day RXTE→NICER gap otherwise mis-classifies J0537 as exponential — the methodological point) |
| Voyager 2 PRA ice-giant rotation periods (F9) | `jansky_research.vgpra` | ✅/➖ blind Lomb-Scargle of the Voyager PRA total-power flux (Uranus + Neptune encounter volumes, PDS-PPI) → a **controlled NULL**: recovers a clean injected rotation in synthetic tests (same wide 14–20 h window) but **recovers neither real period** (Uranus peak 18.4 h wanders across sub-bands; Neptune rails the search bound) — the auroral total-power flux isn't a clean rotational sinusoid, so the historical beaming/magnetic-longitude modelling was essential; the ~2 h flyby precision is hundreds× coarser than the 28-s HST shift (Lamy+2025) regardless |
| e-Callisto megaconstellation RFI trend (F17) | `jansky_research.rfitrend` | ✅/➖ a burst-immune, gain-cancelling narrowband-UEM-line metric over the continuous e-Callisto archive 2012–2026 (intrinsic Starlink lines 125/135/150/175 MHz, Di Vruno+2023; GRAVES 143.05 excluded); primary + secondary metrics recover an injected trend offline → **a systematics-limited NULL**: the two line-sampling stations trend in **opposite signs** (HUMAIN +0.45/yr rising, ALMATY falling), so the pipeline's cross-station coherence test returns `coherent_rise=False` — no Starlink attribution; HUMAIN flagged only as a satellite-pass-gated follow-up candidate |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `jansky_research.rmdipole` | ✅/➖ RM sky **isotropic at dipole order in its core**: the significant power dipole is carried entirely by the top-1% \|residual\| tail (clip → p 0.001→0.93); no apex within 80° of the CMB; real-footprint injection recovered (honest null) |
| Uniform Cat-2 repeater timing census (CHIME/FRB) | `jansky_research.frbwait` | ✅/➖ first one-statistic census of all 83 repeaters: **anchor 20180916B re-found at 16.33 d** (107 cycles, p=0.001, duty 0.21); median k=0.83, 3 clustered; the 2 other p≤0.01 peaks have ≤5 cycles — labelled epoch degeneracies, **no new period claims** (honest) |
| Lensed-repeater search in Cat 2 (CHIME/FRB) | `jansky_research.frblens` | ✅/➖ first catalogue-level fixed-delay search: **0/33 detections → lensed fraction < 0.091** (95%, injection-scoped); documents the day-scramble false-positive mode the real data exposed (honest methods lesson) |
| torch-dsp: coherent dedispersion + RFI + FFA in pure PyTorch | `jansky_research.torchdsp` | ✅/➖ three empty niches filled; **CHIME baseband burst re-dedispersed, S/N peaks at its catalogue DM on a ROCm GPU**; SK/SumThreshold byte-identical to the CPU oracle (sequential); FFA 10.5× on GPU; Crab period re-find = honest null (2.1-s file) |
| Radio survey of 56 WD-pulsar candidates (RACS+VLASS) | `jansky_research.wdpulsar` | ✅/➖ first systematic radio search of the Pelisoli+2025 list: **0/51 candidates detected** (I or V) to a median 3σ V limit **0.41 mJy**; AR Sco control re-found (4.2 mJy, circular); J1912−4410 itself undetected — the duty-cycle caveat made concrete (honest null) |
| Environment-split HI mass function (FASHI DR1) | `jansky_research.fashienv` | ✅/➖ first environment split of the FASHI HIMF (41,741 sources): **void knee suppressed −0.26 dex (2.9σ) vs walls**, an independent FAST measurement of the ALFALFA void HIMF (Moorman+2014); group knee survivor-biased +0.19 dex (stated); DR1 first leg (DR2 embargoed to ~Aug 2026), absolute slope caveated |
| SBI population inference for RACS Stokes-V emitters | `jansky_research.svsbi` | ✅/➖ first SBI of the M-dwarf coherent-emitter population: **first calibrated beaming-fraction posterior f_beam = 0.30 (90% CI 0.07–0.53)**, SBC-validated coverage; LF weakly constrained (beaming–luminosity degeneracy, honest); NPE trained on the ROCm GPU |
| LPT catalogue v3 + Stokes-V forced photometry | `jansky_research.lptv` | ✅/➖ extends the LPT catalogue to **16 members** (3 verified 2026 rows); the ~78-min binary boundary **still not significant at N=16** (p=0.52); first systematic multi-epoch forced Stokes-V at all LPT positions → **1 secure + 1 candidate single-epoch circular burst detection** + a uniform V-limit table (0.47 mJy median; persistent circular pol not a class property; confusion-vetoed, honest) |
| Type III synthesis: corona → 0.4 AU (4 instruments) | `jansky_research.type3synthesis` | ✅ unified drift-to-distance ladder; **geometric check on the model distance** (same-event r=0.989) |

A long run of recover-a-known validations and methodology contributions, two honest negatives (the USS
candidates and the SETI "Voyager detection"), one mixed result (VLASS catalogue variability is
artifact-dominated but image-confirms a real variable star), and **zero overclaims that survived
review** — the science-reviewer caught the USS candidates evaporating against the authoritative
catalog, the SETI "detection" being an instrument artifact, and (every slice) at least one real
physics/citation/statistics fix before write-up. The negatives are arguably the most instructive part.
Each slice's honest assessment is in `survey/*-findings.md`.

## What's next

New slice ideas live in **[`fable-ideas.md`](fable-ideas.md)** (2026-07-05, a 12-agent deep
re-scan of the open-data landscape) — it supersedes the shortlist in
`survey/opportunity-scan-2026-07.md`, and both retire the earlier "reliable no-auth sources are
largely used up" assumption. Suggested first moves there: the SPICE-RACS RM-dipole test (F1 —
data already on disk), the CHIME/FRB Catalog 2 bundle (F2+F5 — gated on the archive recovering;
most time-sensitive), and the `torch-dsp` GPU suite (F6). `plans/29-lotss-deep-144mhz-counts.md`
remains scoped (it should migrate to LoTSS DR3 — see F33); plan 28's single-pulse science was
absorbed into the merged `torchfdmt` slice (`plans/34`). Once the rooftop station (below) produces
calibrated spectra, self-collected data joins the public-archive slices — `fable-ideas.md`
carries a station track (S1–S8) for that too.

## Papers

Every slice is written up as its own **AASTeX paper** under `papers/<slice>/` (authored by Joseph
Barbere, with Claude credited via an AI-use disclosure + a `\software{}` citation — an AI/LLM is not
an eligible author). Each is built reproducibly with containerized tectonic, takes every headline
number from a pipeline-generated `generated/macros.tex` (no figure typed by hand), and is honest
about what it is — mostly recover-a-known validations and methodology, with two honest negatives:

| Paper | `papers/…` | Framing |
|-------|-----------|---------|
| FRB burst statistics, validated on CHIME/FRB Cat 1 | `frbstats/` | validation (tested, reproducible tool) |
| Recovering FRB 20180916B's 16.35-day period | `frbperiod/` | validation |
| The flat inner Milky Way rotation curve from LAB HI | `hi/` | validation |
| A SETI drift-search benchmark + Voyager-1 null | `driftsearch/` | benchmark + honest negative |
| TGSS×NVSS USS selection is dominated by the flux scale | `spectra/` | cautionary negative |
| VLASS multi-epoch variability: a 703 deg² census + FK Com | `vlass/` | methodology + validation (recovers FK Com) |
| Three-frequency curvature selection of peaked-spectrum sources | `peaked/` | methodology + two recover-a-known validations |
| Measuring the turnover: southern peaked sources from GLEAM-X + RACS | `southern/` | methodology + measured-turnover candidate list |
| The AGN radio–optical offset excess and its alignment with the jet (ICRF3 × Gaia × MOJAVE) | `offsets/` | reproduction (excess + jet alignment) + reproducible catalogue |
| The steep radio spectra of pulsars from the ATNF catalogue | `pulsarspec/` | reproduction + MSP/normal comparison |
| Image-plane stacking with injection-recovery: SDSS quasars in VLASS-SE | `stacking/` | methodology + calibrated population-mean flux |
| Multi-decade parsec-scale VLBI variability from Astrogeo | `vlbi/` | control-floor method + recover-a-known (OJ 287, BL Lac) |
| A solar type III exciter speed from an e-Callisto dynamic spectrum | `solarbursts/` | method + recover-a-known |
| The Galactic Faraday rotation sky from the Taylor et al. 2009 RM catalogue | `rmsky/` | reproduction (plane enhancement + sign antisymmetry) |
| The pulsar P–Ṗ diagram from the ATNF catalogue | `ppdot/` | reproduction + population classes |
| Tracking an inner-heliosphere type III beam with Wind/WAVES | `windwaves/` | method + recover-a-known (to the Alfvén surface) |
| Tracking a type III beam to 0.4 AU with STEREO/WAVES | `swaves/` | method + recover-a-known (genuinely interplanetary) |
| 3D triangulation of a type III source with STEREO-A+B direction-finding | `triangulate/` | method + independent geometric-vs-plasma distance cross-check |
| Recovering the canonical 1.4 GHz Euclidean source counts from NVSS | `sourcecounts/` | reproduction (Hopkins 2003) |
| A type III beam from the corona to 0.4 AU, geometrically validated | `type3synthesis/` | synthesis + same-event geometric check on the model distance |
| A reproducible RACS Stokes-V coherent-emitter pipeline (and single-epoch limits) | `stokesv/` | tooling + honest single-epoch/variability result |
| A streaming e-Callisto burst-ingest pipeline with cross-station coincidence QC | `ecallisto_pipeline/` | automation pattern + coincidence-vetted burst events (rejects single-station RFI) |
| A coverage-corrected type III occurrence census (method + recover-a-known) | `ecallisto_census/` | census statistic + recover-a-known validation toward a multi-cycle census |
| A pure-PyTorch Fast DM Transform (device-portable dedispersion) | `torchfdmt/` | tool + oracle validation + real Crab recover-a-known + honest CPU/GPU benchmark |
| Two-epoch RACS Stokes-V forced photometry of nearby M dwarfs | `stokesv_discovery/` | method + GJ 65 variability recovery + upper-limit census |
| A provenance-carrying LPT population catalogue | `lpt/` | verified table + regenerable P–Ṗ statistics (novelty scoped vs the review's own diagram) |
| LPT catalogue v3 + Stokes-V forced photometry | `lptv/` | 3 verified 2026 rows (N=16) + first systematic multi-epoch forced-V limit table at all LPT positions |
| Cassini SKR proximity census | `skr/` | rotation-period anchor (0.05%) + 1/r² sensitivity-null bounding of the SKR occurrence-vs-range trend to a ~1.4× near-null |
| OVRO-LWA type II burst census | `typeii/` | slow-drift detector + in-memory streaming of the whole 2024–2026 archive → an honest null (blind spectral census is false-positive dominated) |
| e-Callisto megaconstellation RFI trend | `rfitrend/` | notch-robust narrowband-UEM-line metric over 2012–2026 + pipeline cross-station coherence verdict → a systematics-limited null (stations trend in opposite signs; no Starlink attribution) |
| Voyager 2 PRA ice-giant rotation periods | `vgpra/` | blind Lomb-Scargle of the PRA flux + rotation-block bootstrap → a controlled null (synthetic recovers an injected period; neither real Uranus/Neptune period is recovered — historical geometric modelling was essential) |
| PTE-II per-source giant-pulse census | `pte2/` | floor-robust giant-pulse excess test over 363 Parkes pulsars + ATNF Ė cross-match → an honest null (detection-power-limited heavy-tail fraction; no Ė trend; tails too steep for classic giant pulses) |
| JBO glitch waiting-time classification | `glitchpop/` | monitoring-gap-robust per-pulsar waiting-time classification of the live JBO catalogue (exp/quasi-periodic) + post-2018 change table; recover-a-known passes (J0537 & Vela quasi-periodic, gap excision required) |
| RM structure functions from SPICE-RACS DR1 | `rmstructure/` | method + recover-a-known + bounded high-\|b\| estimate |
| Jovian DAM occurrence from Juno/Waves | `junodam/` | census method + proximity result + reduced Io-region contrast from orbit |
| The first RM dipole/isotropy test (SPICE-RACS DR2) | `rmdipole/` | method + injection validation + honest isotropy result (tail-carried anisotropy disclosed as systematics) |
| Uniform Cat-2 repeater timing census | `frbwait/` | anchor recovery + population k census + honest no-new-periods verdict |
| Lensed-repeater search in Cat 2 | `frblens/` | first empirical lensed-fraction limit + transit selection function + null-design lesson |
| torch-dsp: the coherent-DSP suite in pure PyTorch | `torchdsp/` | per-kernel oracle validation + real CHIME/Crab legs on ROCm + honest benchmarks |
| A radio survey of the WD-pulsar candidates | `wdpulsar/` | AR Sco recover-a-known + systematic RACS/VLASS non-detection limit table |
| The environment-split FASHI HI mass function | `fashienv/` | first env split of the FASHI HIMF + injection-validated 1/Vmax + ALFALFA void confirmation |
| SBI for the RACS Stokes-V emitter population | `svsbi/` | first calibrated beaming-fraction posterior + SBC-validated coverage + ROCm-trained NPE |

`make paper` builds every slice's PDF; `make papers-zip` bundles them all into one archive (the same
job runs in CI: the **`release` workflow** compiles every paper with tectonic and, on a `v*` tag,
attaches `jansky-research-papers-<tag>.zip` to a GitHub Release — and uploads it as a workflow
artifact on every manual run). `make arxiv` runs the bundled **`arxiv-submit` skill**
(`.claude/skills/arxiv-submit/`) to assemble and validate an upload package per paper
(`papers/<slice>/arxiv-submission/`: the LaTeX-source tarball with its `.bbl`, plus a `metadata.yaml`
capturing every arXiv submission property and a `CHECKLIST.md`). The orchestration is right-sized: the
static slices build through a server-less **Snakemake** file-DAG (`workflow/Snakefile`, run by
`make figures`), while a **streaming Apache Airflow pipeline on rootless Podman** (`airflow/`) ingests
the frequently-updated e-Callisto archive.

### Where to publish (and where not to)

These papers are mostly **reproductions and honest negatives**, so the venue is matched to the actual
contribution — the *tooling and reproducibility*, not a novelty claim:

- **Software / citable archive:** the toolkit is meant for [JOSS](https://joss.theoj.org) (see
  `joss/paper.md`) and a [Zenodo](https://zenodo.org) DOI on release (`.zenodo.json`, `CITATION.cff`).
- **A short note in the literature:** the frbstats validation is condensed to a
  [Research Note of the AAS](https://journals.aas.org/research-notes/) (`papers/frbstats/rnaas.tex`,
  built by `make paper`).
- **arXiv:** reserved for the genuine-novelty, real-data papers (full-repo triage 2026-07-11, PR #115) —
  `lptv/` (the lead: the first uniform multi-epoch forced Stokes-V survey of the long-period-transient
  class, with a secure single-epoch circular detection), `junodam/` (a from-orbit Juno/Waves DAM census
  showing the apparent range-occurrence "law" is a 1/r² detection effect), and `fashienv/` (the first
  environment-split FASHI DR1 HI mass function) — all three packaged and validated; next up are
  `frblens/` (first catalogue-level lensed-repeater search) and `rmstructure/` (first structure function
  of SPICE-RACS DR2). The pre-triage shortlist (`type3synthesis/`, `vlass/`, `peaked/`, `triangulate/`)
  is demoted to repo + Zenodo — recover-a-knowns and method demos, not discoveries. Software-pattern
  papers (`frbstats/`, `torchdsp/`, `torchfdmt/`, `ecallisto_pipeline/`) are JOSS candidates, not
  science preprints. The pure reproductions/negatives are **not** posted as a preprint batch — arXiv
  moderation expects a contribution, and "I reproduced a known result" belongs in the repo + Zenodo.
  The account-bound submission walkthrough lives in Joe's personal notes
  (Obsidian vault: `efforts/radio_astronomy/research_paper_todo.md`).

## Quickstart

New here? [`docs/usage.md`](docs/usage.md) is a short install-and-run guide, and
[`CONTRIBUTING.md`](CONTRIBUTING.md) covers how to contribute, report issues, and get
support.

```bash
# Requires the jansky repo checked out next to this one (../jansky) for local dev.
uv sync                                   # env + jansky (from ../jansky)
make test                                 # unit tests (offline, on synthetic fixtures)
make cov                                  # tests + 85% coverage floor

# Run a slice on real public data (each writes results/ + figures + macros into papers/<slice>/):
uv run python -m jansky_research.pipeline     # FRB burst statistics (CHIME catalog)
uv run python -m jansky_research.frbperiod    # FRB repeater periodicity
uv run python -m jansky_research.spectra --ra 180 --dec 30 --radius 3   # USS hunt
uv run python -m jansky_research.driftsearch  # SETI injection-recovery benchmark
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

See `REPRODUCING.md` for the full reproduction, the right-sized-orchestration notes
(Snakemake static / Airflow streaming), and offline mode.

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
workflow. The `jansky` dependency is a local path source (`../jansky`) for cross-repo dev, switching
to the pinned git tag `jansky@v0.1.0` for clean-checkout CI. See `pyproject.toml`.

## Layout

```
jansky-research/
  src/jansky_research/   # the tooling package (tested-helper pattern, 85% floor) — one module per slice
    data.py              # dataset registry + offline synthetic fallback
    frbstats.py spectra.py frbperiod.py driftsearch.py hi.py vlass.py peaked.py southern.py
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
  .claude/skills/        # arxiv-submit, research-publish, casda-cutout-fetch, radio-cutout,
                         #   find-radio-papers, radio-source-lookup, idea-scan, pull-station-data
  .claude/agents/        # dataset-analyst, pipeline-runner, results-interpreter (+ reused jansky)
  plans/                 # numbered slice specs (00-37); the lasting record is each slice's
                         #   survey/*-findings.md + papers/<slice>/. #29 is planned, not started
                         #   (#28 was absorbed into #34/torchfdmt)
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

MIT — see [LICENSE](LICENSE).
