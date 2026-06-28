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
   **Validated on real data.** A live run (RA 190°, Dec +20°, b≈+80°; Epochs 1–3) exposed and fixed
   two real bugs that live verification caught (the agent's notes had them wrong): the VizieR Epoch-1
   `Fpeak` is **µJy** mislabelled "mJy/beam" (×1000 off → would flag every source variable), and the
   CIRADA flags are stored as floats (`"0.0"`) so `int()` silently dropped all 3.4M rows. After the
   fixes, all three epochs read consistent ~mJy medians (1.8/1.9/1.7) and cross-match cleanly (1022
   of 1148 sources in ≥2 epochs over a 2.5° cone). Loud guards now reject an empty epoch or an
   implausible per-epoch median flux.
3. **GATE-2 vetting (built + run).** `vet_candidates` cross-checks each survivor against SIMBAD + NED;
   `run` writes `results/vlass_candidates.csv`. On this field the conservative 3σ-in-both-metrics cut
   yields a **single candidate** (1.3→5.3→0.7 mJy, a one-epoch ×4 spike) with **no SIMBAD/NED
   counterpart** — which for a single-epoch brightening points to a Quick-Look imaging artefact, not a
   transient. Honest first-field outcome: the pipeline works; this candidate does not survive vetting.
4. **Image vetting (done).** `fetch_vlass_cutout` + `image_lightcurve` read the actual QL image peak
   per epoch (CADC SODA, collection VLASS) as the ground truth. On the first field, **all** candidates
   are catalogue artefacts: image light curves are flat (cand 1 catalogue 1.3/5.3/0.7 → image
   4.6/4.7/4.6; cand 2 catalogue 5.5/0.9/– → image 5.3/5.0/4.6). Failure modes: component deblending
   (fixed by `isolated_mask`), cross-match incompleteness, and PyBDSF-vs-image flux inconsistency. See
   `survey/vlass-findings.md` — an honest negative: catalogue-only QL variability is artefact-dominated
   and image vetting is mandatory.
5. **Forced-photometry confirmation (done).** `run` now auto-vets each candidate: `confirm_candidates`
   does forced peak photometry (`measure_image_flux`) at the locked position in every epoch's image,
   recomputes η/V on the image light curve, and flags `image_confirmed` only for genuinely variable
   ones (p<0.01, V>0.3) — reported as `n_image_confirmed` and written to the candidate CSV. Both
   first-field candidates are auto-rejected (forced V = 0.11, 0.16). The candidate list is now
   trustworthy without manual inspection.
6. **Census + completeness (done).** `injection_recovery` measures data-driven completeness vs flare
   amplitude (inject a single-epoch flare into the real steady light curves, re-run the cut); `run`
   reports the curve, the variable fraction, and `n_image_confirmed`, and writes a completeness figure.
   A **703 deg²** census (42 721 sources, 36 956 usable, ~10 min) gives 40 candidates and **2
   image-confirmed variables** after a forced-peak **centring gate** (`center_arcsec`, which removed
   three bright-neighbour box-edge artefacts, 5→2). One survivor is the known radio-active star **FK
   Comae Berenices** (clean ~7→1 mJy decline, 0.9″ match) — a genuine recovery; the other a faint
   centred fader. Completeness still saturates ~50% (3-epoch ceiling), so these are a floor. Efficiency
   fixes: batched `fetch_vlass_cutouts` (one CADC query/candidate), bounded NED + `use_ned`,
   `max_confirm` cap. See `survey/vlass-findings.md`.
   **Next:** write `papers/vlass/` (tool + QL-systematic-aware methodology + the census: a validation
   (FK Com) plus a completeness-bounded limit). A looser cut / more sky / per-epoch forced photometry
   would tighten it — future work.
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
