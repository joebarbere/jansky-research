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
