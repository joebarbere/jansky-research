# 26 — 3D triangulation of a type III source with STEREO-A + STEREO-B

Status: ✅ done (tooling + synthetic recover-a-known + real run + GATE-2 + paper)

## Context

The `swaves`/`windwaves` slices map a type III burst's emission frequency to a heliocentric distance
*through a density model* (observed frequency = harmonic × local plasma frequency; a density model
gives the radius). STEREO/WAVES also ships a Level-3 **direction-finding** product (goniopolarimetry;
Cecconi et al. 2008, Krupar et al. 2012): per time and frequency, the direction of arrival of the
emission in the heliocentric HEEQ frame, plus the spacecraft position. With **two** spacecraft (A
ahead, B behind) each giving a line of sight, the source is located in 3D as the intersection of the
two rays — **no density model**. That makes the geometric distance *independent* of the
plasma-frequency distance, so the two can be cross-checked, and the geometry adds the source longitude
and latitude the 1-D drift cannot. The honest catch is direction-finding noise: a single type III has
a ~60° apparent source size, so per-sample directions scatter by tens of degrees.

## Deliverables

- `src/jansky_research/triangulate.py` — pure-NumPy tooling:
  - `direction_unit`, `mean_direction` (intensity-weighted vector mean — the correct circular mean),
    `triangulate_rays` (least-squares closest point of two rays, with the forward/`t>0` gate and miss
    distance), `triangulate_track` (per-frequency triangulation over a burst window + Leblanc
    cross-check), `synthetic_event` (offline two-spacecraft fixture with injected angular scatter),
    `fetch_stereo_df` (SPDF L3 CDF, `# pragma: no cover`), `run`/`_figure`/`_write_macros`/`_main`.
  - Reuses `windwaves.emission_radius`/`R_AU_RSUN`/`C_KMS` (and so `jansky.solar`) for the plasma model.
- `tests/test_triangulate.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry entry `stereo-l3-df`; `papers/triangulate/` (AASTeX); `survey/triangulate-findings.md`.

## Approach

1. **Tooling + synthetic recover-a-known.** Place a radially outflowing source at a known
   (longitude, latitude); give each spacecraft the true arrival direction plus Gaussian scatter; the
   offline `run` recovers the injected longitude/latitude (≈1–2° at 9° scatter) and the geometric
   distance correlates with the Leblanc distance (r > 0.9).
2. **Real run (2013-05-15, STEREO-A + STEREO-B, 82° baseline).** 38 channels (0.125–1.975 MHz)
   triangulate to HEEQ lon 169° / lat +5° (near ecliptic), r 15–106 R⊙ (0.07–0.49 AU), median miss
   17 R⊙. Geometric vs plasma distance correlate at **r = 0.989**; absolute scale runs ×2.18 high
   (outward DF-noise bias + likely active-region density enhancement).
3. **GATE-2 science review** (science-reviewer): PASS-WITH-FIXES — corrected the DF-noise attribution
   to the real 82° baseline (ratio 1.08 at 9°, 1.58 at 25°), added the trivial-monotonic-correlation
   caveat (r ≈ 0.75 baseline), clarified the arrival-direction convention, added Cecconi 2008 / Krupar
   2014 citations.
4. **Paper** `papers/triangulate/` — a method/reproducibility contribution: the robust outputs are the
   distance *correlation* and the source *direction*; the absolute radial scale is reported as
   upper-biased.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the injected longitude/latitude within the DF noise and gives corr > 0.8.
- (Real-data) `python -m jansky_research.triangulate --date 20130515` reproduces the metrics, figure,
  and macros from the public SPDF L3 CDFs; GATE-2 sign-off before the write-up.
