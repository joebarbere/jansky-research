# 59 — RACS Stokes-V two-epoch census of the Gaia UCD sample

Status: 🔨 in progress (built in a Claude Code session, not by scout; started 2026-10-15) — GATE 0 run 2026-10-15 (scout, ideas.md "Update 2026-10-15 (2)", PR #343):
novelty PASS at the full-text level for the two closest papers; data PASS for the CASDA listing,
**not** for cutout throughput (S0 below). Controls frozen 2026-10-15, before any real V cutout
of a target is opened. Slice name: `ucdv`.

## Context

The blind VLASS Stokes-I UCD search is closed: Perger, Seli & Vida (arXiv:2506.21169, full text,
Sect. 2 and 3.1) matched M7.5–L4 candidates against VLASS in **Stokes I only** and found no
brown-dwarf counterparts; they use no Stokes V and no RACS. RACS-low2 (arXiv:2606.16182, full
text) measures V only at Stokes-I catalogue positions (Sect. 3.3), reports 61 stars (Sect. 5.2.3),
names no UCDs and gives no per-target limits; its Sect. 6 announces a Stokes-V release for all
RACS epochs, which is the scoop clock for anything catalogue-wide. ADS keyword queries (PR #343)
found no targeted forced-V UCD census. UCD radio emission is coherent and often highly circularly
polarized, so forced V at ~0.9 GHz probes a different mechanism and band than the closed 3 GHz
Stokes-I search. Expected yield is near zero: **the two-epoch V limit table is the paper**, and
that framing is committed now.

## Corrections to the 2026-07 version of this plan (from GATE 0)

- The UCD table is **Reylé et al. 2018** (A&A 619, L8; VizieR `J/A+A/619/L8`), Gaia **DR2**
  based, not a table from arXiv:2506.21169. Row count to be fixed from the table itself at step 1
  (14,915 in the abstract of 2506.21169; the VizieR TSV returned 14,918 lines including headers).
- **RACS-low1 has no Stokes-V products on CASDA** (TAP listing, PR #343). The epoch pair is
  **RACS-low2 (887.5 MHz) + RACS-low3 (943.5 MHz)**; RACS-mid1 + mid2 (1367.5 MHz, same band) is
  the fallback if low3 V coverage of the sample is too thin. The 6% frequency gap is a stated
  limitation of the low2/low3 pair, not a variability signal.
- Plan 53 (EIRP sweep) is **not** bundled: scout's arithmetic puts its limits orders above
  narrowband SETI limits, a different axis. It can reuse this slice's I/O later.
- "~a week wall-clock" was not measured. The last CASDA run (`stokesv_discovery`) took 115 min for
  60 targets × 2 epochs and failed on ~52% of target-epoch slots. Linear scaling to ~15k targets
  is ~500 h. Hence S0.

## Frozen controls (written before real target data is opened)

Provenance rule (CLAUDE.md "an oracle that fed a fit may never certify"): anything that sets a
threshold is production-side. Where a control calibrates, a **disjoint** set certifies.
Split: targets and null positions are assigned to set A/B by Gaia `source_id` parity (even = A).

- **S0 — throughput and scope gate (first real action; touches no target V pixel values).**
  Fetch I+V cutout pairs for 40 targets drawn at random (seed 59) from set B, both epochs, with
  (a) the existing one-cutout-per-job path and (b) a batched SODA job (several CIRCLE positions per
  image) if CASDA accepts it. Record wall-clock per pair and failure rate; the pixel data from
  these 40 is **not inspected** until the controls below have run. Rule: the census sample is the
  full RACS-footprint UCD set if projected wall-clock ≤ 120 h; otherwise a distance-limited
  subsample `d < D`, with D the largest distance whose projected wall-clock is ≤ 120 h. D is
  chosen from parallax and the throughput number **only** — never from any V or I value.
- **C1 — blank-sky null (calibration on A, certification on B).** Forced photometry
  (`measure_circular_pol(search_arcsec=0)`) at **5 positions per target-epoch**, each 2′ from the
  target, at bearings 72° apart with one random rotation (seed 59), same epoch and cutout; nulls
  are then > 2 beams apart. **n ≥ 400 null measurements** (40 S0 targets × 2 epochs × 5). Pass, on
  B only: fraction of V > 0 within the binomial **99%** interval of 0.5; std of V/σ_V in
  [0.85, 1.20]; count with |V|/σ_V ≥ 5 at most the Gaussian expectation + 3. A positive-only null
  means the measurement is searching, and stops the run.
  *False-failure rate of this rule on pure Gaussian noise (20,000 simulations, 2026-10-15): 0.8% at
  n = 400.* The first frozen version (95% sign interval, one null per target-epoch, n ≈ 80) failed
  a perfect pipeline 6.6% of the time (17% at n = 40); the synthetic offline run tripped it on
  pure noise (std 0.83). Amended before any real data was opened (CLAUDE.md, the `hiblend` gate
  lesson).
- **C2 — injection–recovery (real cutouts, set B).** Inject beam-shaped (header BMAJ/BMIN/BPA)
  V point sources at the target pixel, at 3, 5, 7, 10, 20 σ_V, 200 per level. Pass: completeness
  ≥ 90% at 7 σ_V for the 5 σ_V threshold; median recovered/injected flux in [0.95, 1.05] at
  ≥ 10 σ_V. The completeness curve is reported, not tuned.
- **C3 — recover-a-known.** (a) Every star of the RACS-low2 Paper VIII V catalogue
  (arXiv:2606.16182, Sect. 5.2.3) inside a fetched cutout: forced low2 V flux vs catalogue V flux,
  median ratio in [0.85, 1.15]. (b) WISE J0623−0456 (arXiv:2306.15219) if it lies in a low2/low3
  V image: report its forced V in both epochs. (b) is **reported, not a pass/fail**, because
  its detectability at these epochs is not established (scout had abstract/snippet only; read the
  full text before quoting any flux from it).
- **Detection rule (frozen).** A target-epoch is a candidate only if: |V|/σ_V ≥ 5 at the forced
  pixel, AND |V|/I above the local leakage floor (`stokesv.leakage_floor`, `LEAKAGE_NSIGMA` = 7,
  from bright field sources only) when I/σ_I ≥ 5, AND no I source brighter than 100 mJy within
  3′ (sidelobe/leakage guard). With N ≈ 3×10⁴ target-epochs the Gaussian false-alarm expectation
  at 5σ two-sided is ~0.02. Candidates are named only after `classify_emitter` vetting and image
  inspection; nothing is called a detection on the photometry alone.
- **Cutout size.** 0.06° radius (216″): a 2′ null plus its 30–90″ RMS annulus reaches 210″, and the
  bright-source guard needs I out to 3′. The existing `stokesv` default (0.03°) is too small.
- **Positions.** Gaia DR2 positions propagated to each epoch's observation MJD with
  `stokesv_discovery.epoch_position`: Gaia DR2 epoch 2015.5 to RACS-low2/low3 (2022–2024) is 7–9 yr,
  and nearby UCDs move up to ~1″/yr, i.e. several arcsec — a sizeable fraction of the beam.
  Test: a high-PM injected source at its propagated position is recovered; at the catalogue
  position it is not required to be.
- **Limits.** Per target-epoch 3σ upper limit on |V| = |V_meas| + 3 σ_V (conservative, stated in
  the paper), and on L_ν,V with the Gaia distance. No limit is quoted for a target-epoch that
  failed C1's per-cutout sanity check (annulus RMS undefined or > 5× the epoch median).
- **Two-epoch comparison.** Any inter-epoch V change must survive a 10% per-epoch flux-scale
  term (the `stokesv_discovery` V/I-constant lesson), and the scale is bounded from compact
  I sources common to both cutouts. With low2/low3 the band differs: a change is reported as
  "differs between epochs", never as variability, unless it also exceeds what a spectral index
  of ±3 could produce over 887→943 MHz.

## Deliverables

- `src/jansky_research/ucdv.py`: `load_reyle_ucds` (`# pragma`), `null_positions`,
  `inject_point_source`, `forced_census` (wraps `measure_circular_pol`, forced mode only),
  `apply_detection_rule`, `limit_table`, `two_epoch_compare`, plus a batched-SODA fetcher if S0
  shows it works; `run/_figure/_write_macros/_main` with namespaced `ucdvSyn*`/`ucdvReal*` macros
  and `preserve_live_macros` / `preserve_live_results`.
- Synthetic offline fixture (noise + injected sources) for the tests, 85% floor.
- `results/ucdv_*.json` with `source` set; `survey/ucdv-findings.md`; `papers/ucdv/` RNAAS note
  (software-cites `jansky-research` per CLAUDE.md).

## Order of work

1. Tooling + synthetic tests (no network). 2. S0 on set B (throughput only). 3. C1–C3 on set B.
4. Freeze D (if needed) and fetch the census. 5. Detection rule, limits, comparison.
6. GATE-2 science review, then the paper. Any control added after step 3 is labelled post hoc in
the findings.

## Risks

- **Null result plus leakage near bright sources** → the limit table is the paper; the 100 mJy/3′
  guard and the leakage floor are fixed above.
- **Throughput** → S0 decides the scope before the census starts; a distance-limited census is
  still citable and is the honest version if the full one is not feasible.
- **Duty cycle** → two 15-min snapshots miss bursty emitters; limits apply to those epochs only.
- **Scoop** → the announced all-epoch RACS V release (2606.16182 Sect. 6); speed matters less than
  the per-target limit table, which a blind catalogue does not provide.
