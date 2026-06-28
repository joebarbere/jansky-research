# 13 — VLASS three-epoch radio variability catalogue

Status: 🚧 in progress (tooling)

## Context

The VLA Sky Survey (VLASS; Lacy et al. 2020) imaged the sky north of Dec −40° at 2–4 GHz, 2.5″
resolution, in **three epochs** (2017–2024). The CIRADA Quick-Look component catalogues
(Gordon et al. 2021) give ~3.4M components per epoch, free over HTTP/VO from CADC/CIRADA.
Cross-matching a source across epochs by position yields a 2–3-point radio light curve, from which
the standard transient-survey variability diagnostics select variable/transient candidates.

The published "census of variable radio sources at 3 GHz" uses only Epochs 1–2. **Epoch 3 is public
and not yet incorporated** into a three-epoch population study — the priority window this slice aims
at. The deliverable is a reproducible tool + a vetted variable/transient **candidate catalogue**,
not a discovery claim (see `survey/data-source-scan.md` for why this was chosen over the other
shortlisted directions).

## Deliverables

- `src/jansky_research/vlass.py` — tested, CPU-only variability tooling:
  - `eta_metric` / `v_metric` — the two canonical transient-survey statistics (de Vries 2004;
    Scheers 2011; Swinbank 2015): η = weighted reduced χ² vs. a constant flux (significance), and
    V = coefficient of variation (amplitude).
  - `debiased_modulation_index`, `variability_metrics` (η, V, m, m_debiased, χ², p, mean flux).
  - `crossmatch_epochs` — positional N-epoch cross-match (astropy) → per-source light curves.
  - `select_candidates` — the 2-D (log η, log V) outlier selection (Rowlinson et al. 2019).
  - `synthetic_epochs` — offline fixture: a steady population + an injected variable subset with
    known truth labels, so the tests recover the injected variables with low contamination.
  - `fetch_vlass_epoch` (network) + `run(offline=...)` writing `results/vlass_metrics.json` and an
    η–V candidate figure.
- `tests/test_vlass.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for the CIRADA VLASS Quick-Look catalogues.

## Approach

1. **Tooling (this step).** Pure-NumPy/SciPy metrics + astropy cross-match + the 2-D selection,
   validated on a synthetic 3-epoch population (steady sources at η≈1, V≈measurement error;
   injected variables at high η and high V). Offline `run` recovers the injected variables.
2. **Real-data fetch (done).** `run(center=(ra,dec), radius_deg=...)` fetches each epoch and
   cross-matches within 2.5″. There is **no single TAP** for all epochs: Epoch 1 is queried by cone
   from VizieR TAP (`J/ApJS/255/30/comp`, Gordon et al. 2021); Epochs 2–3 are the bulk NRAO
   catalogues (CSV.gz / FITS), cached and region-filtered locally, with the per-version schema and
   quality cuts (`Duplicate_flag<2, Quality_flag∈{0,4}, S_Code≠E` for E1/E2; `Flag=0, S_Code≠E` for
   E3). Each epoch's peak flux is put on the common Perley-Butler scale by `apply_flux_scale`
   (VLASS Memos 13/22) with the ~7% residual scale scatter folded into the errors — without this the
   epoch offset alone manufactures variables. Needs the `vlass` extra (`pyvo`).
   **Still to run:** execute on a real field and cross-check survivors against known variable classes
   (AGN, stars, pulsars) and artefacts (sidelobes near bright sources, ghosts).
3. **GATE-2 science review** before any write-up — the candidate list must survive the known VLASS
   QL caveats (the ~10–15% QL peak-flux underestimate in early epochs, CLEAN/snapshot bias,
   component-vs-source blending) before a single source is called variable.
4. **Write-up** as `papers/vlass/` once the real run + review are done.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers ≳80% of injected variables with low false-positive contamination.
- η and V match hand-computed values on a constant and a known-variable light curve.
- (Real-data, later) candidates survive the VLASS QL caveats and a literature cross-check; GATE-2
  sign-off before write-up.

The GPU is not needed for the metrics/cross-match (I/O- and CPU-bound); it is reserved for an
optional later ML classification of the candidate list.
