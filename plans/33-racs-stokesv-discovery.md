# 33 — RACS Stokes-V discovery slice: forced target-list photometry + two-epoch V variability

Status: ✅ done — GATE 0s passed (mid1/mid2 4.2-yr pair replaced nonexistent low1 V); tooling PR #75; real leg + GJ 65 10σ two-epoch V recovery + upper-limit census + GATE-2-fixed paper merged PR #78

## Context

The `stokesv` slice (#15, merged) delivered the *methods*: CASDA auth + SODA cutouts
(`_casda_session`, `fetch_racs_cutout` with the noiseMap-excluding `_racs_science_mask`), the
off-axis-leakage floor (`leakage_floor`, 7× median |V/I|), forced photometry
(`forced_photometry_recover`), proper-motion vetting (`proper_motion_confirm`), and the honest
finding that **single-epoch V is variability-limited** (catalogue |V/I|=0.90 → image 0.03 for a
flaring emitter). This plan builds the *discovery* half those methods were built for.

RACS-low2 **Paper VIII** (arXiv:2606.16182, 2026-06) published the survey team's **blind** V
catalogue (221 V detections: 61 stars, 85 pulsars, 43 AGN) — and, verified by the 2026-07
literature scan, explicitly left two angles unclaimed:

1. **Forced V photometry on a curated late-type-star input list** — measuring |V|/I (or a 3σ upper
   limit) at every known nearby-M-dwarf/UCD position pushes *below* the blind 5σ extraction
   threshold, where coherent bursts caught at partial amplitude live. Each genuinely new
   V-detected star is a citable find.
2. **RACS-low1 (2020) vs RACS-low2 (2025) two-epoch V variability** — no published two-epoch V
   comparison exists, *not even for Paper VIII's own 61 stars*. A ~5-yr baseline turns our
   documented single-epoch variability limitation into the measurement: coherent emitters are
   expected to appear/disappear between epochs.

Both angles complement (not duplicate) Pritchard+2021/2024, Driessen+2024 SRSC, and Paper VIII.
CASDA is verified working; the TAP endpoint moved to `https://casda.csiro.au/casda_vo_tools/tap`.

## Deliverables

- `src/jansky_research/stokesv_discovery.py` — composing `stokesv` (no re-implementation):
  - `build_target_list()` (network, `# pragma: no cover`) — curated nearby M-dwarf/UCD list with
    Gaia DR3 positions, parallaxes, and **proper motions** (GATE 0b picks the source), plus the
    Paper VIII star list for the two-epoch leg.
  - `epoch_position()` — propagate each target to the RACS-low1 / low2 observation epochs with
    its Gaia PM (nearby M dwarfs move arcsec/yr; RACS pixels are ~2.5″ — mandatory, and
    `proper_motion_confirm` already holds the machinery).
  - `forced_epoch_photometry()` — per target per epoch: I+V forced peak photometry on SODA
    cutouts (reuse `fetch_racs_cutout` + `forced_photometry_recover`), per-field leakage floor.
  - `two_epoch_variability()` — pairwise metrics (ΔV/σ, fractional change, appear/disappear
    flags) for the target list and the Paper VIII 61 stars; adapt `vlass.variability_metrics`
    where it fits two epochs.
  - `select_candidates()` — V above the per-epoch leakage floor + PM-consistent + I/V handedness
    sanity (`handedness`, `classify_emitter` reuse).
  - `synthetic_epoch_pair()` — offline fixture: inject steady / flaring / leakage-only V sources
    into two synthetic epochs (extend `synthetic_field`); selection + variability round-trip with
    known completeness/purity.
  - `run(offline=...)`, `_figure` (V/I vs I with floors + the two-epoch ΔV plane), `_write_macros`
    (macro union so the paper builds offline and real), `_main`.
- `tests/test_stokesv_discovery.py` — synthetic-fixture tests to the 85% floor; mocked-CASDA test
  for the real-path integration (the `stokesv` pattern); no network.
- `data.py` registry entries for the target-list source(s); `papers/stokesv_discovery/` (AASTeX);
  `survey/stokesv-discovery-findings.md`; Makefile SLICES/figures/reproduce wiring
  (`CASDA_USERNAME`-gated like `stokesv`).

## Approach

0. **GATE 0a — RACS-low1 V image availability (do FIRST).** `stokesv` pulled RACS-low V cutouts —
   confirm specifically that **both** low1 *and* low2 epochs' `image.v.*` products are
   SODA-servable per field, and record each epoch's observation MJD (needed for PM propagation
   and for the "5-yr baseline" claim). If low1 V is not servable, the slice degrades to the
   forced-photometry angle alone — still publishable; say so and continue.
1. **GATE 0b — target list.** Pick one machine-readable nearby-M-dwarf source and record it in
   `data.py`: candidates are the CNS5 (Golovin+2023, VizieR `J/A+A/670/A19`), Reylé+2021 10 pc
   sample (`J/A+A/650/A201`), or a Gaia DR3 parallax+colour cut (ϖ > 40 mas, M-dwarf locus).
   Cross-check against the Sydney Radio Star Catalogue so *known* radio stars are labelled, not
   claimed as new. GATE 0c: a machine-readable Paper VIII star table (CDS/VizieR or the paper's
   own supplementary); if unavailable, transcribe the 61-star table with provenance noted.
2. **Tooling + synthetic recover-a-known.** Two synthetic epochs with injected (i) steady V
   emitters, (ii) single-epoch flares, (iii) leakage-only contaminants at the floor: the
   selection must recover (i)+(ii) and reject (iii); the variability step must flag exactly (ii).
   Report completeness/purity as the offline metrics.
3. **Real run (bounded).** Forced I+V photometry for the N≈50–100 nearest targets (Dec < +45°,
   inside RACS-low coverage) in both epochs + the Paper VIII 61 stars; per-field leakage floors;
   candidate table with PM confirmation. Bound the CASDA load (the `stokesv` slice's timeout
   lesson): batch cutouts, retry-with-relogin (already implemented), and run the long leg as a
   background job writing incremental CSV (the census-ingest pattern).
4. **GATE-2 science review** — the candidate list must survive the leakage floor and PM vetting;
   every "new" detection must be checked against SRSC/Paper VIII/SIMBAD before the word "new" is
   used; upper limits reported for non-detections; the variability claim framed against the
   two-epoch sampling caveat (two epochs cannot distinguish flaring from secular change).
5. **Write-up** `papers/stokesv_discovery/` — a candidate list + two-epoch variability census,
   with the honest framing: sub-threshold forced photometry finds *candidates* for follow-up,
   not confirmed coherent emitters.

## Verification

- `make test` / `cov` green on synthetic fixtures (85% floor); `ruff` + `mypy` clean; paper
  builds offline via the macro union in CI.
- Offline: selection completeness/purity on the injected epoch pair within tolerance; the
  variability step flags exactly the injected flares.
- Real: every candidate above the leakage floor, PM-confirmed, and cross-checked against
  SRSC/Paper VIII/SIMBAD; `make reproduce` (with `CASDA_USERNAME`) regenerates the tables,
  figures, and macros. GATE-2 sign-off.

## Risks & mitigations

- **Scoop risk (highest) →** this is the time-sensitive opener; build the bounded version first
  (50–100 targets + the 61 Paper VIII stars), publish honestly, extend later.
- **CASDA flakiness / long real leg →** retry-with-relogin exists; incremental background CSV;
  never run the full-list fetch inside a test or a single foreground session.
- **RACS-low1 V not servable →** GATE 0a degrades the slice to forced photometry only (still
  unclaimed); document.
- **Leakage masquerading as detections →** per-field floor (not survey-global), reject anything
  within a beam-edge annulus, and keep the leakage-only synthetic contaminants in the fixture.
- **Positional errors at 5-yr baseline →** PM propagation is mandatory; verify on a
  high-PM control (e.g. Kapteyn's star region if in coverage).
