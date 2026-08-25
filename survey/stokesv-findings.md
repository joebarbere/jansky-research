# Findings — RACS Stokes-V coherent-emitter selection (in progress)

`jansky_research.stokesv` selects coherent radio emitters (radio stars / ultracool & brown dwarfs /
pulsars) by their circular polarization $|V|/I$ in the ASKAP RACS survey. Circular polarization is a
near-unambiguous flag of a coherent emission process, so a high $|V|/I$ at a stellar position is a
clean finder; the dominant false positive is instrumental Stokes-I$\to$V leakage, defeated by a
per-region leakage floor ($7\times$ the median field $|V/I|$) plus a Gaia proper-motion confirmation.
This file records the data-access reality and the first (credential-free) real-data result.

## Data access: VizieR for catalogues, CASDA (authenticated) for V cutouts

- **CASDA catalogue TAP is unexpectedly unusable from a script here.** Every ADQL query against the
  `AS110.*` RACS catalogue tables (and even `ivoa.obscore` / `TAP_SCHEMA`) on
  `https://casda.csiro.au/casda_vo_tools/tap` returns `relation "..." does not exist` (PostgreSQL
  lowercases the case-sensitive schema), via both the sync and async endpoints, authenticated or not.
  The documented `FROM AS110.racs_mid_components_v01 WHERE 1=CONTAINS(...)` form fails identically.
  CASDA's async results listing also omits the `xlink:href`, so `pyvo` cannot fetch results. The
  supported `astroquery.casda` path only returns ObsCore *image* products, not catalogue rows.
- **Resolution:** use **VizieR** for the catalogue layer (it mirrors RACS and is what the other slices
  already use) and reserve **authenticated CASDA SODA cutouts** for the Stokes-V *images* the forced
  photometry needs (that login path is verified working — `Authentication successful`).
- VizieR catalogues wired: RACS-low DR1 Stokes-I (Hale+2021, `J/other/PASA/38.58`, `Fpk`/`Noise`) and
  the Sydney Radio Star Catalogue (Driessen+2024, `J/other/PASA/41.84`) — whose `radio` sub-table
  carries per-detection Stokes I **and V** peak fluxes (`SpeakI`, `SpeakV`, `localrmsV`) across
  RACS-low/mid/high and VLASS.

## Recover-a-known (credential-free): the SRSC V-detected stars classify as coherent emitters

Running the selection helpers over the SRSC radio table, on the **176** RACS detections that have both
a Stokes I and a Stokes V peak flux:

| quantity | value |
|---|---|
| median $|V|/I$ | **0.67** (deeply circularly polarized — coherent emission) |
| 90th-percentile $|V|/I$ | 0.96 |
| classified `circular` or `highly_circular` by `classify_emitter` | **176/176 = 100%** |
| of those, `highly_circular` ($|V|/I\ge0.3$) | 160 |

This confirms `fractional_circular_pol` + `classify_emitter` behave correctly on real known coherent
emitters (`validate_srsc`). It is a *positive-set* check — the SRSC radio table only lists detections,
so the complementary purity test (rejecting the unpolarised + leakage population) needs the
forced-photometry leg below.

## Forced Stokes-V/I photometry (the science core, tested offline)

`measure_circular_pol` performs the forced measurement at a locked target position: it finds the
Stokes-I peak within a small search box, then reads Stokes V **at that same pixel** (V is signed) —
the physically correct measurement for a point-like coherent emitter, where the circular-polarization
peak coincides with the total-intensity peak. It returns $I$, $V$ (signed), per-Stokes annulus RMS,
$|V|/I$, and the I-peak offset, so a non-detection becomes an honest upper limit rather than a miss.
Tested on a synthetic SIN-projection image: recovers an injected $|V|/I=0.4$ and the LCP sign. It is
archive-agnostic — it runs on a cutout array + WCS from whichever service serves the RACS images.

## CASDA VO-service outage (blocks the live cutout fetch right now)

The authenticated login works (`Authentication successful`), but **every CASDA programmatic query
service is erroring** as of this run, so the standard `query_region` → `cutout` discovery flow cannot
run:

| CASDA service | endpoint | result |
|---|---|---|
| catalogue TAP | `casda_vo_tools/tap` | `relation "..." does not exist` (sync + async, authed + not) |
| ObsCore TAP | `ivoa.obscore` | same `relation does not exist` |
| SIA2 image query | `casda_vo_tools/sia2/query` | HTTP 500 NPE: `Cannot invoke "java.util.Map.size()" because "m" is null` |
| `*/availability` | — | report `available=true` (the front ends are up; the query backend is not) |

So the live image fetch is blocked by **CASDA infrastructure, not credentials or our code**.

### Data Central is not a Stokes-V alternative (investigated)

`datacentral.org.au` hosts RACS and has a working POST cutout API
(`/api/services/cutout/`, FITS output via `fits:true`, band selection by integer PK), but its RACS
holdings are **RACS Low1 / Low3 Data Release 1 — Stokes I total intensity only**. There is **no
Stokes V** in any of its 186 imaging bands (the only radio bands are GLEAM/GLEAM-X; RACS data releases
73/74 carry no Stokes-V band, and RACS bands aren't even enumerated in the cutout schema). So Data
Central cannot supply the V images this slice needs. (It *is* a clean GLEAM/GLEAM-X Stokes-I cutout
source — relevant to the queued southern-curvature runner-up, not here.) **CASDA remains the only
public RACS Stokes-V image source**, and it must recover before the forced-photometry leg can run.

### CADC mirror also checked — RACS Stokes I only (no V)

CADC *does* host a `RACS` collection (799 planes), but a CAOM2 TAP query shows it is **`polarization_states = /I/` only — Stokes I, no V**. So neither public mirror (CADC, Data Central)
carries RACS Stokes V; CASDA is the sole source, confirmed from three independent archives. (Bonus:
CADC RACS-I *is* usable via the `radio-cutout` skill's CADC SODA path if a slice ever needs RACS total
intensity.)

### Browser-automation skill (`.claude/skills/casda-cutout-fetch/`): works, but the outage reaches the web service too

As a route around the broken VO APIs, a Playwright skill drives the **web** CASDA Cutout Service with
the OPAL login. Verified live: the **OPAL login, modal handling, navigation, results-URL construction,
results parsing, and all the guards work** — the script authenticates and reaches the results page
correctly. But it **could not download a FITS**, because the CASDA Cutout Service returns
**`RACS-low DR1 0  RACS-mid DR1 0  RACS-high DR1 0` for every position — including the service's own
example target, PSR B1919+21**. So the discovery-backend outage that takes down TAP/SIA2 also takes
down the web cutout service; no automation can download what the backend won't surface. The skill is
kept (it exits with a distinct code 5 for "outage" vs 8 for "broken automation") and should succeed
once CASDA recovers — the honest "best-effort: kept the working skill, recorded the blocking finding".

## CASDA RECOVERED (2026-06) — the forced-photometry leg is now unblocked

Re-tested live: the CASDA query backend is up (`query_region` returns RACS products), the OPAL login
works, and — the thing that was ERROR-ing — **SODA cutout staging now succeeds** (returns real
`cutout-*.fits` download URLs). So `fetch_racs_cutout` is wired and the forced-photometry leg runs.

### Two real complications, and the honest finish

1. **The `noiseMap`/`meanMap` trap.** The CASDA image query returns, per field, a `noiseMap.image.i.…`
   and `meanMap.image.i.…` alongside the science `image.i.…` (all carry `.i.`/`.v.` + `restored`).
   Selecting the science `image.{i,v}.` product is essential — the noise map otherwise reads as a flat
   ~0.2 mJy field and manufactures a null (`_racs_science_mask` handles this).
2. **Single-epoch V is variability-limited.** Coherent stellar emission is *bursting*; the catalogued
   RACS-LOW V detections come from whichever epoch caught each star flaring. Forced photometry on a
   single RACS-low DR1 snapshot **recovers Stokes I well at the known position** (validating the CASDA
   cutout + forced-photometry pipeline) but recovers significant **V only for the subset caught in a
   polarised state** — an honest lower bound set by the duty cycle, not a pipeline failure. (Example:
   a catalogued |V/I|=0.90 emitter with I=19 mJy → image I=10.5 mJy recovered, image |V/I|=0.03.)
3. **CASDA auth is intermittently flaky** — the datalink step occasionally returns HTTP 401; the fetch
   retries with a fresh login.

**Finish (chosen: honest single-epoch):** `fetch_racs_cutout` (CASDA SODA, science-image filter, retry)
+ `forced_photometry_recover` over the brightest RACS-LOW emitters → report I recovered (median
image/catalogue ratio) and the variability-limited V fraction. Framed as methods + tooling + honest
limits; the multi-epoch blind survey (leakage floor over the field, VAST+RACS epochs) is the natural
next step the tooling is ready for. Paper at `papers/stokesv/`.

## Honest caveats so far

- The leakage floor and the negative (purity) test require the field V measurements from the images;
  the credential-free result validates only the positive set.
- The SRSC `e_SpeakV`/`localrmsV` columns are sparsely populated, so a uniform V-SNR cut is not yet
  applied in the validation (the $|V|/I$ and classification are robust regardless).
- ASKAP's absolute $V$ sign convention varies by pipeline/epoch — handedness is recorded but not
  physically interpreted without the per-epoch convention.

## Full referee round (2026-08-25): MAJOR REVISION, 14 findings, three BLOCKERs

The referee reconstructed the fifteen uncommitted rows from the committed figure's vector
coordinates (calibration verified against the committed medians) -- and the reconstruction is
what exposed the blockers.

**BLOCKER 1: two of the nine "circular detections" have image |V|/I of 568% and 135%** --
physically impossible for a single source -- and `classify_emitter` has no upper bound, so both
count as `highly_circular`. Excluding them: 7/15, not 9/15. The claim "the I recovery confirms
the source and the measurement are correct" is falsified for at least two targets.

**BLOCKER 2: the I and V cutouts are selected by two INDEPENDENT unordered CASDA queries** --
no obs_collection/band/obs_id constraint, first row of whatever matches
`image.<stokes>.*restored*conv` (RACS-low/mid/high all match) -- and V is then indexed with I's
WCS with no grid check. The "single-epoch RACS-low DR1" provenance is enforced by no line of
code; a cross-band/epoch pairing economically explains both |V|/I>1 and the I scatter. The
correct pattern already exists in scripts/stokesv_discovery_real.py (obs_id-grouped TAP, same
observation required).

**BLOCKER 3: the entire real result exists only inside a PDF.** No per-target CSV;
forced_photometry_recover's rows (cat/img fluxes, offsets) are reduced to five scalars and
discarded -- the innerrc lesson verbatim.

**MAJORs:** the 0.92 median I ratio hides min 0.068 / max 2.75 (0.48 dex scatter, 8/15 outside
2x) -- "confirms the pipeline works" cannot be carried; the variability interpretation ("often
the VAST monitoring caught each star flaring") is FALSE as written -- every target is selected
`Survey == RACS-LOW`, and if catalogue row and image are the same observation the variability
explanation is unavailable (epochs never recorded); none of the three advertised gates
(leakage floor, V-SNR, proper motion) runs on the real leg -- "significant" means |V|/I >= 0.06,
a constant that appears only in a figure label; the synthetic validation cannot fail (injected
at 0.2-0.8 vs a 0.006 leakage population; PM gate "confirmed" by zero exercised cases; all 500
targets bright); the results file claims real provenance ("RACS-low DR1 (CASDA)") for a
half-synthetic file, against the repo's own mixed-marker rules; 9/15 quoted with no interval
and \svFracVcirc typesets as "0.6 of them".

**MINOR/NIT:** the 12" peak-search mode is used while the prose says "the physically correct
forced measurement" (the sibling's flagged phrasing) and the per-row offsets that would test
the mitigation are computed and thrown away; the 15-target cap is code, not data, and the
parent population count is unrecorded; the 6% threshold is duplicated as magic numbers in
classifier and figure and is uncalibrated against off-axis leakage (4 of 9 "detections" sit at
9-23%, inside the documented near-edge leakage range); ~/.casda_pw undocumented and the query
non-deterministic; pritchard2021 title wording; caption says "panels" for one panel; the 7x
convention uncited.

**Status: fixes pending** (needs one CASDA re-run with obs_id-pinned pairs; credentials:
CASDA_USERNAME explicit + ~/.casda_pw).

**Status: RESOLVED (2026-08-25).** The same-observation re-run replaced the paper's central
claim with what the data actually say, and every blocker's mechanism was confirmed.

**Blockers 1+2, closed together:** `fetch_racs_low_pair` pins each target's I and V cutouts to
one ObsCore `obs_id` (RACS collection, RACS-low band, taylor.0 restored conv; identifiers and
`t_min` committed per target) and asserts the grids match. Re-measured, the physically
impossible ratios are gone from the statistics: 9 of 15 targets are valid, and the 6 invalid
are Stokes-I non-detections (<3 sigma) in that observation -- single-epoch states, reported as
such and excluded from every number. The decisive surprise: the I ratio scatter SURVIVES the
pinned pairing (median 0.56, range 0.10-2.56), so it is real -- single-epoch stellar
variability and flux-convention differences, not the query bug -- and the paper no longer
claims "Stokes I is recovered" beyond the median.

**Blocker 3:** results/stokesv_targets.csv commits every row: fluxes, both local rms values,
V-SNR, offsets (median 1.85", max 7.9" -- the search stayed on-source), validity, obs_id,
t_min, and both product filenames.

**The headline:** significant V ( |V|/sigma_V >= 5 AND |V|/I >= 6% ) in **4 of 9 valid
targets (44.4%, Wilson 29-61%)**, identical at a 5% or 10% threshold (the four detections
carry V-SNR 12-47). The parent population (19) is committed, so 15 measured is nearly the
whole sample. The "variability-limited duty cycle" interpretation is withdrawn to what the
evidence supports: a single-observation occurrence rate, with the committed t_min enabling the
epoch comparison a duty-cycle claim needs.

**The rest:** the V-SNR gate now runs on real data ("significant" means something); the
source string names both legs (allowlisted, per the singlepulse precedent); \svFracVcirc (the
"0.6 of them" typo carrier) is retired for \svFracVcircPct with a Wilson interval; the figure
excludes invalid targets and requires the same gates as the statistics, and was redrawn from
the committed CSV; "the physically correct forced measurement" is replaced by a plain
description of the 12" search with committed offsets; the synthetic validation's scope is
stated (injections far from the floor; zero PM rejections exercised).
