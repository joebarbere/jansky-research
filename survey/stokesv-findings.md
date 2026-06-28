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

## Status and what's next

- **Done (this step, credential-free):** tested helpers (leakage floor, $|V|/I$ selection, PM
  confirmation, classification); VizieR fetchers (`fetch_radio_stars`, `fetch_racs_i`,
  `fetch_radio_star_measurements`); the pure `match_targets_to_radio` cross-match; the `validate_srsc`
  recover-a-known above.
- **Next (authenticated):** forced photometry of a curated late-type-star / UCD target list in RACS
  Stokes-V CASDA SODA cutouts — measure $|V|/I$ (or a $3$–$4\sigma$ upper limit) *below* the blind
  extraction threshold, estimate the per-beam leakage floor from local field sources, apply the V-SNR
  and proper-motion gates, and vet with SIMBAD. This is the new-findings arm: candidates the blind
  RACS V catalogues miss. Then GATE-2 and `papers/stokesv/`.

## Honest caveats so far

- The leakage floor and the negative (purity) test require the field V measurements from the images;
  the credential-free result validates only the positive set.
- The SRSC `e_SpeakV`/`localrmsV` columns are sparsely populated, so a uniform V-SNR cut is not yet
  applied in the validation (the $|V|/I$ and classification are robust regardless).
- ASKAP's absolute $V$ sign convention varies by pipeline/epoch — handedness is recorded but not
  physically interpreted without the per-epoch convention.
