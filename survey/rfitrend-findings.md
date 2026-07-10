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
