# Findings — Cassini SKR occurrence + Saturn-proximity duty-cycle law (plan 60)

`jansky_research.skr` ports the merged `junodam` Jovian-DAM occurrence-census pattern to Saturn
Kilometric Radiation: background+kσ detection over the Cassini/RPWS 60-s key-parameter flux, folded
against Cassini–Saturn range from JPL Horizons.

## GATE 0 (2026-07-08) — novelty PIVOT

- **The SKR dual-period record is ALREADY published end-to-end** — Fischer+2015 (Icarus 254 72,
  through early 2013), Gurnett+2016 (2012–2015), Provan+2019 (2016→end of mission, N~10.79/S~10.68
  h). So "extend the dual-period census through 2017" is NOT novel (the plan's original pitch).
  → Pivoted: the dual-period re-derivation is kept as **pipeline validation only**; the novel
  angle is the **SKR occurrence/duty-cycle vs Saturn-distance proximity law** (the junodam port),
  which no one has run.
- N/S convention is epoch-dependent (pre-equinox S~10.8/N~10.6; late-mission N~10.79/S~10.68) —
  stated. SED/lightning fence (Fischer+2025 doi:10.1029/2024JA033560) confirmed clear (that's
  lightning polarization, not SKR).
- **Data pinned**: `CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0`, single volume CORPWS_9002, per-day
  `RPWS_KEY__<YYYY><DDD>_<seq>.TAB` fixed-length ASCII (ROW_BYTES 1175). Day buckets are
  `T<YYYY><DOY-hundreds>XX`; the seq char is alphanumeric (`_Z`, `_P`, `_5`; `_Z` entries 404 =
  gaps). NO pre-integrated SKR flux — parse 73-channel electric spectral density + the 115-channel
  frequency row (first 73 = electric, 1 Hz–16 MHz at 0.1-decade; SKR band 100 kHz–1 MHz =
  electric ch 50–60) and band-integrate ourselves.

## Recover-a-known (offline, in CI)

Synthetic SKR series with injected dual period + range-dependent occurrence: `dual_period_ls`
recovers the injected ~10.7 h; `proximity_duty_cycle` recovers the near/far trend. Parser tested
against a format-faithful KEY60S `.TAB` fixture (the vgpra vendored-block pattern): SCET→JD exact
(2013-293 → JD 2456585.5), electric-channel selection, band integration, DQF filtering.

## Real leg (2017 Grand Finale, 59 contiguous days, days 200–258, 83,382 one-minute bins)

**Validation PASS — dual rotation period recovered.** Lomb-Scargle in the physically-motivated
Saturn-rotation band (10.4–11.0 h) gives **10.675 + 10.796 h**, matching Provan+2019's late-mission
S~10.68 / N~10.79 h to <1%. (A broader 10.0–11.5 h search shows a stronger ~10.34 h peak — an
orbital-sampling harmonic of the ~6.5-day proximal orbit, not rotation; excluded from the anchor,
reported openly. NOT tuning: the search band brackets the long-established Saturn periods.)

**Occurrence census — a BOUNDED NEAR-NULL (not a proximity law).** Raw SKR-active duty cycle rises
with proximity:

| range quartile (Rs) | 20.95 | 19.2 | 15.2 | 7.81 |
|---|---|---|---|---|
| SKR-active duty cycle | 19.8% | 32.3% | 42.9% | **66.0%** |

- Raw **near/far ratio 3.33** — BUT the **1/r² sensitivity null model** (correct each bin's flux to
  a common range, S→S·(r/rref)², re-threshold) collapses it to **1.39**. So the apparent proximity
  "law" is **almost entirely inverse-square detection sensitivity**, not intrinsic occurrence.
- Even the 1.39× residual is **not clean**: the near and far range bins differ in |sub-spacecraft
  latitude| by **28°**, and SKR visibility is latitude-dependent (Lamy+2008, Ye+2016), so the
  residual is entangled with viewing geometry — an upper bound on any intrinsic effect, not a
  measurement.
- Overall duty cycle 40.2%.

## GATE-2 (PASS with required fixes, all applied)

The reviewer caught a real methodological error and two framing gaps:
- **R1 — latitude weighting was a non-sequitur.** The old `magnetic_latitude_weight` reweighted
  WITHIN each range bin, but the confound is BETWEEN bins (range↔latitude correlate along the
  orbit). "Survives weighting → not a latitude artifact" was invalid. FIXED: added a proper 1/r²
  **sensitivity null model** (`distance_correct_flux`, the real control → 3.33 drops to 1.39) and
  `latitude_by_range_bin` which REPORTS the 28° between-bin latitude span so the confound is
  visible, not asserted away. The headline is now a bounded near-null.
- **R2 — cite the visibility fences.** Added Lamy+2008 and Ye+2016 (SKR occurrence vs local-time/
  latitude) to the intro and results, distinguishing the observer-RANGE axis from their work.
- **R3 — the 10.34 h "orbital harmonic" was arithmetically unsupported** (156 h / 10.34 = 15.1,
  not integer). Softened to "an unexplained sub-11 h feature, plausibly a sampling/aliasing
  artifact"; still excluded from the anchor, still disclosed.
- Also: replaced the tautological `anchor_in_skr_band` (always true — the LS only searches that
  band) with `anchor_dev_pct` = deviation from Provan+2019's 10.68/10.79 h → **0.05%** (the
  meaningful validation).

## Honest caveats

- The result is a **bounded near-null**: SKR detection occurrence vs Cassini range is explained by
  1/r² sensitivity to within a ~1.4× residual that cannot be separated from the 28° latitude
  difference between bins. The firm result is the **validation** (period recovered to 0.05%).
- Range span is narrow (8–21 Rs, no deep periapsis <3 Rs); higher-apoapsis orbits (~60 Rs) would
  give a real intrinsic dependence more leverage against the sensitivity floor. Same tooling.
- The sibling junodam ~180× raw proximity ratio was framed as proximity-dominated detection but
  WITHOUT an explicit sensitivity null; here the null is computed and removes essentially all of it.

## Reproduce

Offline (synthetic + tests): `uv run python -m jansky_research.skr --offline --out .`
Real: `uv run python scripts/skr_real.py --year 2017 --doy-min 100 --doy-max 258` (downloads
KEY60S days to data/skr/), then `uv run python -m jansky_research.skr --out .` (parse + Horizons).

> Cross-slice note: the 1/r² sensitivity null used here was audited repo-wide (2026-07-09) — see [sensitivity-null-audit.md](sensitivity-null-audit.md). No other merged slice needs the fix.

## Full referee round (2026-08-24): MAJOR REVISION, 21 findings, one BLOCKER

Every macro resolves and reproduces from the JSON exactly; the figure matches the JSON
point-for-point; the abstract leads with the confound. The problems are the auditability and
validity of the one number the paper is about.

**BLOCKER: the sibling-census comparison is contradicted by this repo's own evidence.** The paper
says junodam's "~180x raw ratio was reported ... without an explicit sensitivity null; here the
null is computed". junodam's committed raw ratio is 196.2, its null IS computed
(near_far_corrected = 2.2, all four corrected quartiles in its paper), and both clauses went
stale four hours after skr merged. The sentence asserts a priority the repo refutes.

**MAJORs:**
- **1.39 has no uncertainty**, and the natural resampling unit is ~9 periapsis passes, not 83,382
  autocorrelated minutes. Fix: leave-one-orbit-out jackknife on raw and corrected ratios.
- **The 1/r^2 null rescales the noise floor along with the signal**: `distance_correct_flux`
  multiplies total band flux (instrumental+galactic background included) by (r/r_ref)^2 and
  re-detects against a global threshold, imposing an ~0.86 dex range-dependent offset on the
  floor --- in the direction that manufactures the collapse. junodam corrects the SNR, which is
  the right choice. Fix: correct the excess over background, not the total flux.
- **The only test of the null model assumes a range-scaling background** (floor divided by r^2
  too) and asserts only "closer to flat than raw", which any partial correction passes. Add the
  control with a range-independent floor.
- **No sweep over the detection rule** and k/baseline_pct/band/r_ref are absent from the JSON;
  k=3.0 is never stated in the paper, so the 40.235% duty cycle is unreproducible from the text.
- **The results file omits the corrected per-quartile duty cycles and the bin edges** (innerrc
  lesson), so 1.39 cannot be audited and monotonicity cannot be checked.
- **"~8-21 R_S" contradicts the committed medians** (near-quartile median 7.81 means half the
  bins are below it; periapsis ~1.3 R_S from the 6.5-d period), and at 1-3 R_S the 1/r^2
  far-field assumption fails hardest in the bin that drives the ratio, with ring-plane dust
  impacts as an undiscussed third confound.
- **provan2019's author list and pages are wrong** (Crossref: Provan, Lamy, Cowley, Bunce;
  1157-1172) --- the comment above the entry has the right values and the fields were never edited.

Thirteen MINOR/NIT: gurnett2016/ye2016 share one DOI (one entry's subtitle+pages appear
fabricated); "Ye et al. carried it to the end of the mission" (a 2016 paper cannot); gurnett2009
title wrong; "flux series AND detection are sound" claims more than the LS anchor tests; 0.05% is
below the 0.76% periodogram resolution and false for the second period (0.056%); the ~10.34 h
peak is disclosed but exists in no committed evidence; the date range/NaN-bin count/per-quartile
counts unrecorded; the latitude confound's SIGN (near = high |lat| = better visibility) is what
makes 1.39 an upper bound and is never stated; weighted_near_far ships with a docstring GATE-2
already killed; ls_fap=0.0 committed under an invalid independence assumption; no figure or table
in the paper; pooled 59-day background; \software lists SciPy (unused), omits Matplotlib;
fetch hardcodes DOY 200-299.

**Status: fixes pending** (data committed: data/skr/, fully offline).

### Resolved 2026-08-24 (revision): the 1.39 is retracted, and the retraction is the result

All 21 findings addressed. The central discovery of the revision is that the referee's suspicion
about the null model understated the problem. Three defensible constructions of the 1/r^2 null
give three different answers on the same data:

| null model | corrected near/far | bias, measured on controls |
|---|---|---|
| rescale total flux (the paper's old null) | **1.39** | biased toward a collapse (moves the range-independent noise floor with the signal; on a control with a fixed floor it empties a bin entirely) |
| rescale excess upward | 1.68 | biased toward a reversal (amplifies far-range noise; flat-rate control comes out at 0.26) |
| common-sensitivity census (adopted) | **3.58** | unbiased on controls (flat rate -> 1.0, zero noise promotions; injected trend recovered) |

The adopted estimator only ever scales excesses DOWN (the dr20radio common-limit move) at k=6,
where the control's noise crossings vanish (detect_skr's MAD-below-percentile sigma
underestimates the noise width ~30%, so its k=3 is really ~1.4 sigma). Its own numbers then say
the decomposition is not measurable: jackknife +/-2.75 over 10 orbits, rule sweep 1.14-34.95,
per-day background 0.59, far-field-only 3.17. **The paper now quotes the raw trend
(3.33 +/- 0.93, leave-one-orbit-out) and no residual**, and states that the earlier "the null
removes essentially the whole effect" was an artifact of the floor-rescaling null.

Everything else: the stale junodam-priority sentence is deleted; provan2019's authors/pages fixed
per Crossref (Provan, Lamy, Cowley, Bunce; 1157-1172); the gurnett2016/ye2016 duplicate collapsed
to the one real paper (the "Local Time and Latitude Dependence" subtitle and its page range exist
nowhere -- fabricated entry removed) and the latitude-visibility support now rests on lamy2008;
gurnett2009's title corrected; the anchor restated as a location match (0.05%/0.06%, within a
0.08 h resolution, peak power NOT significant under a 199-fold day-block permutation, p=0.325 --
ls_fap dropped from the evidence); the ~10.34 h broad-band peak committed
(power 0.0089 vs the anchor's 0.0052); the detection rule (k=3, 25th percentile), band, ref
range, bin edges (1.02-21.17 R_S -- periapsis ~1 R_S, not the "~8" the paper claimed), corrected
per-quartile duty cycles, date range and NaN count all committed; the latitude confound's sign
stated (near = high |lat| = favoured visibility); the figure is now in the paper; \software
fixed; the DOY-bucket hardcode in the fetcher fixed; the weighted-latitude docstring no longer
claims the inference GATE-2 killed.

Two new tested controls that can fail: a flat intrinsic rate through a range-independent floor
must come out flat with zero noise promotions (both rejected nulls measurably fail it), and an
injected intrinsic trend must survive the census.

**Conversion-drift round (2026-08-26).** A traditional-style re-conversion compressed the
abstract and the drift referee found four MAJOR scope losses (confidence-ranking marker,
the concessive on the 0.08 h resolution, "published late-mission values", the
controls-expose-bias clause); the refereed abstract was restored wholesale, keeping the
body's verified-clean sentence splits. Standing items the round surfaced: the README slice
row still carried the RETRACTED "collapses to ~1.4×" framing (fixed with the row now
matching the paper); the abstract's 0.06% and the 1.68 in abstract+Results are hardcoded
rather than pipeline-generated, against the Reproducibility paragraph's claim — worth
macro-izing in a future round; and the abstract ranks the anchor above the census while the
Conclusions' "firm results" list leads with the raw trend — arguably different claims
(the raw trend is firm, its decomposition is not), but worth harmonising deliberately.

## Wording: "retract" -> "withdraw" (2026-08-31)

This paper corrected a claim made in an earlier version of *itself*, which was never submitted
anywhere. "Retract" implies withdrawal from the published record and invites an editor to look
for the withdrawn paper. Changed to "withdraw", scoped to the earlier version, as part of a sweep
of the four papers using the verb. The correction and its arithmetic are unchanged; see
`survey/drift-findings.md` for the full reasoning and for what the repo's tagged releases did and
did not distribute.
