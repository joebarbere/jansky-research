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
