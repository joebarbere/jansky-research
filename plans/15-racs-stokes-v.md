# 15 — RACS Stokes-V coherent radio emitters (radio stars / UCDs / pulsars)

Status: 🚧 in progress (tooling + credential-free catalogue path + SRSC recover-a-known done; forced-photometry leg next)

## Context

Circular polarization is a near-unambiguous flag of *coherent* radio emission (electron-cyclotron
maser, coherent plasma emission): extragalactic AGN are at most a fraction of a percent circularly
polarized, but flaring M-dwarfs, ultracool/brown dwarfs, magnetic chemically-peculiar stars, RS CVn
binaries, and pulsars reach tens of percent. So a Stokes-V (|V|/I) selection in a wide-field survey is
a clean coherent-emitter finder — and each genuinely new V-detected star is a real, citable find, not
a reproduction. This is the slice with the highest *new-findings* potential in the
`survey/new-findings-scan.md` scout.

The ASKAP RACS survey covers the southern sky (Dec ≲ +48°) in three bands and now releases Stokes-V
images/catalogues publicly on CASDA. The blind single-epoch V search is already done by the survey
team (Pritchard et al. 2021, MNRAS 502, 5438; the RACS-low2 Paper VIII V catalogue, 2026), so we do
**not** repeat it. The open, complementary angle for an independent researcher is a **targeted
forced-photometry V measurement of a curated late-type-star / UCD input list** — measuring |V|/I (or a
3σ upper limit) at each known stellar position, which reaches *below* the blind 5σ extraction
threshold and turns non-detections into burst-rate constraints. Optionally, a RACS-low1 vs RACS-low2
two-epoch V comparison flags transient circular-polarization events (M-dwarf radio duty cycle ~8.5%,
so single-epoch non-detections are expected).

## Reuse (this slice composes existing tooling)

- `vlass.forced_photometry` / `vlass.fetch_*` cutout pattern — forced peak photometry at a locked
  position, repointed from CADC SODA to CASDA SODA (near-direct).
- `spectra.crossmatch` — positional matching (radio ↔ stellar input list, and candidate ↔ leakage
  reference AGN in the same beam).
- `spectra.spectral_index` — the Stokes-I two-point index across RACS bands (gyrosynchrotron vs ECME
  discriminant) for any multi-band detection.
- `vlass.vet_candidates` (SIMBAD/NED) — identify/annotate each candidate.

## Deliverables

- `src/jansky_research/stokesv.py` — the new, tested logic (pure NumPy/astropy):
  - `fractional_circular_pol` — |V|/I with error propagation.
  - `leakage_floor` — robust per-region Stokes-I→V leakage estimate from unpolarised reference sources
    (median |V/I|), returning the credible threshold `n_sigma_floor × median` (default 7×, the RACS
    convention) so candidates must clear *leakage*, not just image noise.
  - `select_circular_pol` — candidates with |V|/I above both the leakage floor and an SNR cut, sign of
    V (handedness) recorded.
  - `proper_motion_confirm` — Gaia-DR3-style PM check: is the radio position consistent with the
    star's position propagated to the radio epoch (genuine), or a chance coincidence?
  - `classify_emitter` — coarse class from |V|/I, I-band spectral index, and stellar type.
  - `synthetic_field` — offline fixture: an unpolarised AGN population carrying realistic
    position-dependent leakage, injected circularly-polarised stars (high |V/I|) with a matched
    proper-motion stellar catalogue, so tests recover the injected stars and reject leakage.
  - `run(offline=...)` — synthetic offline; real CASDA forced-photometry online; writes
    results + figure + macros.
  - `validate_known` (later) — recover the Pritchard 2021 / Sydney Radio Star Catalogue stars.
- `tests/test_stokesv.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** Leakage-floor estimator + |V|/I selection + PM confirmation, validated on
   a synthetic field that recovers injected polarised stars at low contamination and rejects the
   leakage population.
2. **Real data (next).** Forced-photometry V at a curated M-dwarf/UCD target list (LSPM / CARMENES /
   Reylé 10-pc) in RACS-low2 + RACS-mid V cutouts from CASDA SODA; per-beam leakage floor from local
   AGN; SIMBAD/Gaia vetting + PM confirmation. Honest depth + leakage caveats.
3. **Validation.** Recover the known Pritchard 2021 / SRSC radio stars from the public V images as a
   recover-a-known; quantify the false-positive (leakage/sidelobe) rejection.
4. **GATE-2** before write-up — candidates survive the leakage, sidelobe, sign-convention, and
   chance-match caveats.
5. **Write-up** as `papers/stokesv/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers ≳80% of injected polarised stars with low leakage contamination, and rejects
  the unpolarised+leakage population (≈0 false positives above the 7× floor).
- `leakage_floor`, `fractional_circular_pol`, `proper_motion_confirm` match hand-computed values.
- (Real-data, later) recovers known RACS V radio stars; GATE-2 sign-off before write-up.
