# 88 — `dr20radio`: the first radio-counterpart census of the SDSS-V DR20 Black Hole Mapper

Status: ✅ complete 2026-08-07 (all increments + GATE-2; see survey/dr20radio-findings.md) — GATE 0 done 2026-08-07 (catalog pinned, novelty verified, rights
landmine scoped out; details below). Standing remainder: same-week ADS re-search before the
first commit. (Plan 87 = lineconf, on its local branch.)

## Context

SDSS-V DR20 (released 2026-08-03; Almeida et al., arXiv:2607.26149) delivers the Black Hole
Mapper's ~500k spectroscopically confirmed quasars/galaxies — including **the first optical
SDSS spectra ever taken from the southern hemisphere** (BOSS on the du Pont at Las Campanas).
No radio cross-match of any SDSS-V BHM catalog exists; the crowded SDSS×radio genre is built
entirely on legacy catalogs (Milliquas/DR16Q×VLASS arXiv:2105.12985; RHzQCat DR16Q×RACS at
z>3 arXiv:2512.03415; red-quasar×VLASS arXiv:2603.24739). The clean first: **VLASS cannot see
south of δ = −40°, and southern SDSS quasar spectra did not exist before DR20 — so "SDSS
spectroscopic quasars × RACS south of −40°" was categorically impossible until this release.**
The slice: a uniform radio-counterpart census of the DR20 quasar table against VLASS (north,
already on local disk from the `vlass` slice) and RACS (south), with radio-detection fraction
vs redshift/luminosity and the north/south selection contrast as the headline products.

**GATE 0 (done 2026-08-07).** *Catalog:* `spAll-lite-v6_2_1-allepoch.fits.gz`, **177 MiB**,
at `data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/summary/allepoch/` (verified listing);
datamodel confirms `RACAT/DECCAT`, `Z`, `ZWARNING`, `CLASS` (select `QSO`), `SUBCLASS`,
`FIRSTCARTON`, `SDSS_ID`, and **`OBS` (APO/LCO)** — the north/south flag is a first-class
column, no inference needed. Quasar-property VACs exist (Wu et al. spectral properties;
Dwelly et al. visual inspection) for later enrichment. *Novelty:* no DR20 or BHM×radio work
found (ADS/arXiv, Aug 2026). One caveat to handle, not a blocker: DR20 carries BHM
open-fiber cartons `openfibertargets_bhm_racsradio_boss` / `..._lofarradio_boss` — a small
set of fibers targeted BECAUSE they are radio sources; flag via `FIRSTCARTON` and exclude
from (or split out of) the detection-fraction statistics, else the census is circular for
those objects. *Rights:* the eROSITA/SPIDERS X-ray leg is limited to the German-consortium
hemisphere (l = 180–360°; Merloni et al. arXiv:2401.17274, unchanged at DR2 2026-07-31) —
**deferred entirely**; this slice is optical×radio only. *Radio:* VLASS δ > −40° (2–4 GHz,
2.5″, ~120 μJy/beam/epoch; three-epoch catalogs local); RACS-low −80°≤δ≤+30° (887.5 MHz,
25″) + RACS-mid to δ≤+49° (1367.5 MHz, ~11″); overlap band −40°→+30° enables a
cross-survey consistency check. Literature match radii cluster 1.5–6″; adopt ~2.5″ primary
with a wider secondary pass for extended sources.

## Deliverables

- `src/jansky_research/dr20radio.py`:
  - `fetch_spall` (SAS download, resumable, `# pragma`) + `load_quasars` (FITS → table:
    `CLASS=='QSO'`, `ZWARNING==0` handling, `OBS` split, `FIRSTCARTON` radio-carton flag).
  - `match_vlass` (local QL epoch catalogs, all three epochs; 2.5″ primary / wider secondary
    with the false-match rate measured by shifted-position trials, not assumed).
  - `match_racs` (RACS-low/mid via CASDA TAP or Data Central; same false-match calibration).
  - `detection_fraction` (vs z, vs magnitude/luminosity proxy, vs OBS hemisphere; radio-carton
    objects excluded/split), `northsouth_contrast` (the selection-function-aware comparison —
    VLASS 3 GHz vs RACS 0.9/1.4 GHz depths differ; convert to common luminosity limits
    honestly rather than comparing raw fractions), `multi_epoch_flags` (VLASS E1/E2/E3
    variability for matched objects — flagged as catalog-level with the `vlass` slice's
    documented artifact caveat).
  - `synthetic_catalog` (offline fixture with known injected counterpart fractions and a
    known radio-targeted subset → recover-a-known incl. the circularity exclusion).
  - `run/_figure/_write_macros/_main`; committed-real-results from day one (evidence JSONs
    force-tracked; macros only from committed evidence).
- Tests to the 85% floor; `papers/dr20radio/`; `survey/dr20radio-findings.md`; wiring.

## Approach

1. Tooling + synthetic recover-a-known (fractions recovered; radio-carton circularity
   exclusion verified to matter in the synthetic case).
2. Real leg A (local, no bandwidth to speak of): spAll-lite fetch (177 MB) → quasar table →
   VLASS match against the on-disk epoch catalogs → northern census.
3. Real leg B: RACS match (TAP-scale queries, modest) → the southern-first census — the
   headline: radio properties of the first southern SDSS quasars.
4. False-match calibration by position-shift trials; detection fractions with binomial
   errors; the north/south contrast at matched luminosity limits.
5. GATE-2: selection-function honesty (BHM targeting ≠ uniform quasar sample — state what
   the parent sample is and is not), frequency-mismatch caveats, carton circularity, no
   dichotomy overclaims (report fractions, not verdicts on the radio-loudness debate).
6. Paper (`papers/dr20radio/`): census tables + fractions + the categorical-first framing,
   tightly scoped to DR20-new objects and the southern pairing.

## Verification

Synthetic round-trip recovers injected fractions and the circularity exclusion; the
VLASS-RACS overlap band (−40°→+30°) gives an internal consistency check on match rates;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **BHM targeting is heterogeneous** (eROSITA-selected, variability-selected, filler
  cartons) → the census is per-carton-class aware from the start; the paper reports by
  selection class, never as "the quasar population".
- **Frequency/depth mismatch north vs south** → common-luminosity-limit comparison, both
  raw and matched numbers shown.
- **RACS service availability** (CASDA moods) → Data Central fallback per the
  `radio-cutout` skill's routing; TAP queries are resumable per-chunk.
- **The racsradio carton** → excluded from fractions, reported separately (it is itself a
  small validation set: radio-selected objects should match at ~100%, a nice pipeline check).
