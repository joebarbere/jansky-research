# Findings — e-CALLISTO as an accidental 15-year RFI observatory (plan 54, fable-ideas F17)

`jansky_research.rfitrend`: a burst-immune, gain-cancelling occupancy metric applied to the
continuous e-CALLISTO FITS archive 2012–2026, trending the Starlink unintended-emission (UEM) band
and attributing to the public Starlink on-orbit count. The only published RFI-trend study on the
archive (Prieto/Pérez+2020, SoPh 295:11) is a two-epoch 2012-vs-2019 campaign; this extends it into
a continuous trend across the megaconstellation era.

## GATE 0 (2026-07-10)

- **Novelty PASS**: Pérez+2020 is a two-snapshot (2012, 2019) RFI-occupancy comparison at Spanish
  sites finding a ~2× rise, ending *before* Starlink scaled. No post-2019 e-CALLISTO RFI-trend /
  megaconstellation census exists. Dedicated-instrument Starlink-UEM studies (LOFAR Di Vruno+2023;
  gen-2 Bassa+2024) are calibrated-interferometer campaigns that motivate, not pre-empt, an
  archival-spectrograph trend.
- **Starlink UEM bands in 45–870 MHz**: broadband 110–188 MHz; **intrinsic** narrowband lines at
  **125/135/150/175 MHz** (Di Vruno+2023). The 143.05 MHz feature is reflected GRAVES radar (NOT
  intrinsic — excluded). [Corrected 2026-07-10: an initial 137.05 MHz line was a transcription
  error for the GRAVES 143.05, which is external anyway — GATE-2 catch.]
- **Control**: FM 87.5–108 MHz is the cleanest flat terrestrial control; DAB/TV carry
  analog-switchover step changes.
- **Attribution regressor**: the public planet4589 Starlink cumulative-count time series — no
  Space-Track login needed.
- Reuse: `solarbursts.fetch_ecallisto` (one 15-min FITS), `ecallisto_catalog`, `ecallisto_census`.

## The metric (and a real-data design pivot)

**The systematics problem.** e-CALLISTO is uncalibrated log-power: station gain drifts over years
(would fake a trend) and solar bursts are broadband transients (would mask one). Both are
**common-mode** — they move the whole spectrum together.

- **Burst immunity**: per-channel *time median* over the 15-min sweep = the persistent occupied
  level; a burst/satellite pass occupies a small fraction of the sweep and does not move a median.
  Verified: turning bursts on in up to 90% of synthetic months does not change the recovered slope.
- **Planned gain-cancelling metric**: UEM-band occupancy minus an FM-control band at the same
  station (common-mode gain cancels in the difference).

**The pivot (found on first real contact).** Smoke-testing on real HUMAIN spectra: the station
(focus code 59) **notches out the FM band entirely** — a hardware gap 84→112 MHz, zero sampled FM
channels — and notches other sub-bands too. e-CALLISTO operators configure station-specific
RFI-avoidance notches, so a fixed FM control is **not universally observed**. Response:
- **Primary metric ⇒ the narrowband line-vs-adjacent excess**: the level at each sampled intrinsic
  UEM line (125/135/150/175 MHz, Di Vruno+2023) over its immediately adjacent channels.
  Self-normalizing *within* the UEM band, so it cancels *common-mode* station gain by construction
  and is robust to per-station notches (a line simply drops out if unsampled). It does **not**
  cancel differential local flank contamination — which turns out to be the residual systematic
  that limits the study. Most Starlink-specific spectral feature available.
- **Control ⇒ station-adaptive** (`pick_control_band`): FM if sampled, else a sampled clean band
  below/above the UEM window; reported per station.
- **Config-stability enforced**: a station must keep the *same* sampled UEM lines across its
  retained months, else the differing months are dropped (an instrument-reconfiguration guard).

## Synthetic recover-a-known (offline, in CI)

168 synthetic monthly spectra with a known Starlink-shaped UEM-level rise + per-month common-mode
gain drift + broadband bursts:

- Differential **recovers the injected trend at 0.18/yr** (`diff_trend_p` ≈ 0), correlating with
  the injected Starlink shape at **r = 0.999**.
- A **null control** (two bands both outside the UEM window) does **not** trend (slope 0.001,
  `control_flat` True) — the metric does not manufacture a trend from gain drift.
- **Burst immunity explicit**: recovered slope unchanged at 0% vs 90% burst months.

## The real archive (ran 2026-07-10): a SYSTEMATICS-LIMITED NULL, with one flagged candidate

Monthly-sampled 2012–2026 at HUMAIN, ALMATY, GLASGOW (**286 usable station-months** — finite
line-excess, after config-stability screening). Streamed one small gzipped FITS at a time, in memory.

| station | UEM lines (MHz) | months | line-excess slope/yr | p | r(Starlink) | Pérez 2012→2019 |
|---|---|---|---|---|---|---|
| HUMAIN  | 150, 175       | 161 | **+0.448** | 2e-5  | +0.48 | +7.0 log-units |
| ALMATY  | 125, 135, 150  | 125 | **−0.128** | 0.012 | +0.04 | — |
| GLASGOW | none stable    | 162 | — | — | — | — |

**The two line-sampling stations disagree in SIGN.** HUMAIN's UEM-line excess rises strongly and
significantly (+0.45/yr), positively correlated with the Starlink count (r=0.48), and its raw UEM
occupancy rose ~7 log-units 2012–2013 → 2018–2019 — qualitatively consistent with Pérez+2020's
pre-Starlink rise (but raw gain-contaminated level at a different station, so a consistency note,
not a reproduction). Taken alone HUMAIN looks textbook. **But ALMATY — overhead the same global
constellation — falls significantly (−0.13/yr) AND with essentially zero Starlink correlation
(r=0.04)**, so its trend is plainly not a Starlink signal.

- **A real global megaconstellation signal should raise every station together.** The pipeline's
  cross-station coherence test (`summarize_stations`) finds the significant slopes do **not** agree
  in sign (**1 rising, 1 falling; `coherent_rise` = False**).
- The weakly positive **pooled** slope (+0.31/yr, p<1e-5) with only **r=0.35** Starlink correlation
  is a **sample-size artifact** of the more populous rising station (HUMAIN 161 vs ALMATY 125
  months), not a coherent trend.
- **Consistent with per-station effects dominating**: local RFI-environment change, receiver/antenna
  reconfiguration (uncalibrated, operator-tuned instruments), and station-specific occupancy of the
  flanking channels that normalize each line. This is exactly the systematics limit the
  differential/multi-station design was built to expose.
- **Coherence-test caveat**: with only 2 line-sampling stations (GLASGOW drops out) the test is
  underpowered, and opposite signs do not *strictly* disprove a global signal (the two stations
  sample partly different lines; the self-normalized excess can read as falling if a station's flank
  channels accrue local RFI faster than the line). Since the conclusion is a null, this asymmetry is
  conservative — it can only make us *less* likely to claim a detection.

**Honest bottom line: we claim no Starlink attribution and no megaconstellation detection.** The
archive can measure per-station UEM-band trends but cannot, from spectral data alone, separate a
global constellation signal from co-located local-RFI growth — and the sign incoherence indicates
the latter dominates here. **HUMAIN is flagged as a candidate**: its rise is real, significant, and
Starlink-correlated in the right band; the follow-up that would break the local-RFI degeneracy is
satellite-pass-gated occupancy (does the excess appear only when Starlink is above HUMAIN's
horizon?), which the archive cannot provide. We flag it and claim nothing more. The verdict is
pipeline-generated (`cross_station_signs_agree`, `coherent_rise`), so the artifact carries the
conclusion, not just the prose.

- Data-quality caveat: a handful of archived FITS were server-truncated (astropy warns); the
  median-based metrics are robust to partial time coverage.

## GATE-2 (2026-07-10) — PASS on honest framing, with required fixes (all applied)

The reviewer confirmed the core is an honestly-stated systematics-limited null (no Starlink
attribution; HUMAIN properly hedged as a candidate; pooled slope correctly called an artifact) and
that the `coherent_rise=False` verdict is load-bearing and correctly computed. Required fixes, all
applied:

- **Line-frequency correction (the key catch)**: the initial `137.05 MHz` UEM line was a
  transcription error. Di Vruno+2023's *intrinsic* Starlink narrowband lines are **125/135/150/175
  MHz**; 143.05 MHz is reflected **GRAVES** space-surveillance radar (external), now excluded via
  `GRAVES_MHZ`. Verified against the primary source; the census was re-run — HUMAIN (150/175) is
  unchanged, ALMATY moved to 125/135/150 (and its Starlink correlation dropped to ~0, sharpening
  the null).
- **Overclaim on Pérez "reproduction"** → downgraded to a hedged consistency note (raw
  gain-contaminated level, different station/metric).
- **Coherence-test caveat added**: n=2, different lines sampled, flank contamination can flip a
  sign — the test can fail to detect but not strictly disprove a global signal (conservative for a
  null). Softened "shows … dominates" → "indicates / consistent with".
- **Sample-count fix**: the reported months excluded GLASGOW's all-NaN line series (448 → 286
  usable).
- **Primary-metric validation**: the synthetic now injects narrowband lines, so `line_vs_adjacent`
  (the primary metric) is validated end-to-end, not only the differential cross-check.
- **Citation fix**: Bassa+2024 → A&A 689, L10, doi 10.1051/0004-6361/202451856.
- Cleanups: p-values formatted as upper bounds (no bare `0.0`); removed the mislabeled pooled
  `diff_*` duplication; dead `CONFOUND_MHZ`/`DAB_CONTROL` removed/rewired.

## Reproduce

Offline (metric + synthetic recover-a-known + tests): `uv run python -m jansky_research.rfitrend
--offline --out .`
Real (streamed, in memory): `uv run python scripts/rfitrend_real.py --stations HUMAIN ALMATY
GLASGOW --start 2012 --end 2026` (writes `results/rfitrend_metrics.json`, `is_real=True`).

## Referee round on the style conversion (2026-08-24)

The null, the HUMAIN/ALMATY sign disagreement, the underpowered-discriminator caveat and the "we flag
it, and claim nothing more" close all survived at full strength, and the abstract trim (340 -> 316
words) removed only rhetoric. Fixed:

1. **MAJOR.** The bold signpost "The cross-station coherence test is the verdict" had become "**the
   decisive check**" --- an evaluative claim about the test's evidential power that the paper
   withdraws twelve lines later ("with only two line-sampling stations the test is an underpowered
   global-vs-local discriminator", "it does not prove the absence of any Starlink contribution").
   A test cannot be both decisive and underpowered. Now "the paper's primary discriminator".
2. The abstract's "A genuine global megaconstellation signal **would** raise every station's UEM
   lines together" (body keeps "should"): the abstract carries no version of the caveat that the
   implication can fail, so "would" presented as deductive a premise the paper shows is
   probabilistic. Restored.
3. A sentence split left "This is the gain-contaminated raw level at a different station" with
   *Perez's campaign* as its nearest antecedent, i.e. saying their result was contaminated.

## Three pre-existing defects, fixed here

**1. The abstract attributed the trend three sentences before disclaiming attribution.** "We
attribute the trend to the public Starlink on-orbit count" vs "We therefore claim no Starlink
attribution". Now "We test the trend against".

**2. `\rfRealNMonths` = 286 was attached to three stations but counts two.** 286 = 161 (HUMAIN) +
125 (ALMATY); the three-station total is 448, since GLASGOW's 162 months sample no UEM line. Scoped
to the line analysis.

**3. The selection function read as a screen, and its attrition was unreported.** The three stations
are a hard-coded list (`ECALLISTO_STATIONS`), not the survivors of an archive-wide stability screen
over ~150 stations, and "configuration-stable" described them post hoc. The screen's attrition is
also badly uneven and was invisible: HUMAIN 174 -> 161 (7%), GLASGOW 163 -> 162 (1%), **ALMATY
174 -> 125 (28%)** --- and ALMATY is the station whose falling slope carries the sign disagreement the
null rests on. `n_months_raw` was in the committed metrics all along; it is now emitted as
`\rfReal<ST>NMonthsRaw` and stated in the paper, so the attrition is macro-backed rather than absent.

**Still open:** the synthetic recover-a-known injects only common-mode gain drift and broadband
bursts, both of which the line-vs-adjacent difference cancels **algebraically** --- hence `\rfSynCorr`
= 1.0 and `\rfSynLineCorr` = 1.0. The systematic the paper itself calls decisive, differential local
contamination of a line's *flanking* channels, is never injected: ask what would make this test fail,
and the answer is nothing in the regime that matters. The paper now says what the validation does and
does not cover; adding a flank-contamination arm and reporting the recovered slope bias is the real
fix.

### Resolved 2026-08-24: the flank-contamination arm exists and the validation can now fail

`synthetic_month_stack(flank_rise=...)` injects local RFI into the lines' *flanking* channels only,
as a linear ramp uncorrelated with the Starlink curve --- the systematic the paper names as decisive
and the line-vs-adjacent difference cannot cancel. Measured (committed as `flank_*` in the metrics
and `\rfSynFlankSlope`/`\rfSynFlankBias`):

- clean arm: +0.2404 /yr recovered (the injected rise)
- flank arm: **-0.1907 /yr** --- a bias of -0.4311 /yr that **flips the recovered sign**

The ALMATY mechanism is therefore reproduced end-to-end: a station whose flanks fill in faster than
its lines registers a falling excess under a rising true signal. The paper's validation paragraph
now reports this instead of disclaiming it, and states the consequence --- a station's falling line
excess cannot, by itself, be read as a falling UEM signal. Two tests pin the arm (bias is negative;
sign flips; the committed bias equals the difference of the two committed slopes).

## Full referee round (2026-08-26): MAJOR REVISION, 17 findings, two BLOCKERs

The structure is right and the null is honest: every macro reproduces exactly from the
committed JSON, the figure matches the data, the coherence verdict is pipeline-computed, and
divruno2023/bassa2024/benz2009 all Crossref-verify. What blocks it: the one positive claim is
contradicted by its own committed time series, and the reference carrying the novelty argument
is mis-attributed.

**BLOCKER 1: HUMAIN's "rise" is two single-month instrumental steps, not a trend.** From the
committed line_excess array: a −12.0-unit step at 2018.875 (held 3.3 yr) and a +15.0-unit step
at 2022.208 (held to the end); the paper's fitted total change over 13.4 yr is 6.0 units —
less than half of either step. Within each of the four regimes the Theil–Sen slope is zero or
negative (none rises significantly), and the DOWNWARD step spans exactly the interval
(2019.0–2021.9) when Starlink went 0→~1900 — a signal that cannot behave like this. All
retained months share the same stable_lines, so the config screen cannot catch it; the paper
concedes within-set reconfigurations survive the screen, then exempts HUMAIN. "Whose rise is
real, significant, and Starlink-correlated" must become a description of a piecewise-constant
series with two discontinuities; the regime table belongs in the JSON.

**BLOCKER 2: perez2020 is mis-attributed** — the real paper (Sol. Phys. 295, article 11, DOI
10.1007/s11207-019-1577-5, Crossref-verified) is "Increase in Interference Levels in the
45–870 MHz Band at the Spanish e-CALLISTO Sites over the Years 2012 and 2019" by Prieto,
Bussons Gordo, Rodríguez-Pacheco, et al. — no Pérez-Torres exists on it, the title in refs.bib
is wrong, and there is no DOI. The lawrance2024 pattern, on the single reference carrying "the
only published RFI-trend study on the archive" and a results macro. The hand-typed "~2× rise"
needs checking against the paper's own numbers.

**MAJORs (all measured):**
- NO SENSITIVITY STATEMENT ON THE NULL (the frblens lesson): injecting a common Starlink-shaped
  rise into both real series (circular-shifted residuals, real cadence) and applying the
  paper's own coherent_rise criterion: ~60–70% power against a global signal the size of the
  largest excursion in the data (A≈6–8 units), 28% at half that. "Indicates per-station
  effects dominate" is not earned at 60% power — "is consistent with" is; commit the power
  curve and quote the exclusion amplitude.
- "SAMPLE-SIZE ARTIFACT" IS REFUTED BY THE OBVIOUS TEST: subsampling HUMAIN to ALMATY's 125
  months leaves the pooled slope at 0.261 [0.224, 0.292] vs published 0.306. The pooled slope
  is amplitude-dominated (HUMAIN ±13 units vs ALMATY ±2), not count-dominated; the abstract's
  stated reason for dismissing it is the wrong reason.
- THE REAL RESULTS JSON IS A SPLICE: three flank_* keys only _synthetic_metrics produces were
  appended to the real file (git show 016471d), no _merge marker, and scripts/rfitrend_real.py
  writes with a bare write_text that bypasses preserve_live_results entirely. Values verified
  correct; provenance is the defect (the torchfdmt pattern).
- PER-MONTH FILE PROVENANCE COMPUTED AND DISCARDED: sample_month_metrics returns the file
  name; _real_trend stores only years+excess. Which day/UT each month sampled is unrecorded —
  and that is now the single most important question about the two HUMAIN steps.
- THE ARXIV.YAML OVERRIDE REINSTATES A RETRACTED CLAIM: "configuration-stable stations" (the
  post-hoc description the findings file records as previously fixed) replaces main.tex's
  "long-running", and the override drops "We therefore claim no Starlink attribution." The
  yaml's mechanical check preserved numbers but not claim strength.

**MINOR/NIT:** Theil–Sen intervals discarded (HUMAIN +0.45 [0.27, 0.57]; ALMATY −0.13
[−0.215, −0.023] — "falls significantly" rests on an upper bound of −0.023); ALMATY's window
truncation checked and BOUNDED (HUMAIN on ALMATY's window: +0.546 — the confound does not
drive the null; disclose the spans); the recover-a-known never states the known (truth slopes
0.2382/0.1787, recovery ratios 1.009/1.015 — genuinely good, uncredited); r=1.0 is 1.0 by
construction on the clean arms (say which arm can fail); the flank-arm sweep shows sign
crossover at flank_rise ≈ 0.8× the true line rise and flank_rise=5.0 reproduces ALMATY's
measured pair almost exactly (−0.119/−0.036 vs −0.128/+0.041) — a quantitative match,
unreported; the Starlink regressor is ten hand-typed rounded anchors, nonzero from 2019.1
(paper says zero before 2019.4), clamped after 2026.0, and possibly year-shifted (needs one
planet4589 lookup; +1 yr shift moves r 0.480→0.516, nothing turns on it); the
shape-discriminant (Starlink-shape r=0.480 vs linear-ramp 0.294 vs step 0.005) is the only
quantitative evidence for the HUMAIN flag and is unreported; FM-notch 84/112/focus-59 numbers
uncommitted (control_name="low" does back the usable half); stale 137 MHz prose in module +
driver docstrings and one vacuous test assertion; GLASGOW's "162 of 163 retained" retains an
empty line set (contributes nothing); burst-immunity claim true (0.2382 vs 0.2406, measured)
but in no artifact.

**Verdict: MAJOR REVISION.** The single change: report HUMAIN as a step function, not a trend
— the null gets stronger, the systematics framing gets its sharpest example, and the flagged
candidate is replaced by something the data support.
