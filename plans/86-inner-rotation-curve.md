# 86 — `innerrc`: independent replication of the Sofue & Kohno inner rotation curve

Status: ✅ done 2026-08-05 — **both replication legs land, GATE-2 PASS with fixes** (PRs
#160–#165). Anchor: their ρ_DM = 0.107 arithmetic reproduces exactly (4π convention, tested);
an unconstrained refit of their own table prefers ρ_DM = 0.24 (rms 11.8 vs 13.6), all 8
sensitivity variants in 0.19–0.32 — the low density is degeneracy width, not tension. Raw
HI4PI (1,113 sightlines): Table 2 reproduces at ~4% (8.8 km/s); threshold estimator measured
+17.9 km/s high of Gaussian decomposition (σ calibration 26.6 vs 9.5, bracketing their 15) —
the `hi` caveat closed with a number; E/W asymmetry replicates in period/phase. Five
paper-style figures + `papers/innerrc/main.tex` (compiles; macros from committed evidence
only). Remaining: publication sequencing (Obsidian todo). CO surveys + their-exact-fit
replication remain stretch legs. — GATE 0 done 2026-07-31/08-01: full-text read (all 13 pp), novelty PASS,
data-access verified (details below). Standing remainder: same-week ADS re-search before the
first commit. (Note: the local `slice/lineconf` branch also carries a "plan 85" — renumber it
to 87 when it revives; 85 is atlas3i on main.)

## Context

Sofue & Kohno 2025 (PASJ 77, 1335; arXiv:2509.23581) derive the definitive modern inner
rotation curve (RC) of the Milky Way: terminal-velocity method (TVM) on longitude–velocity
diagrams from **five public surveys** (HI4PI; CfA-Chile, FUGIN/Nobeyama, Nobeyama GC, and
Mopra ¹²CO), with terminal velocities from **Gaussian decomposition** (highest-velocity
component), an iterated velocity-dispersion correction (±15 km/s HI, ±5 km/s CO), Gaussian
running-average binning (δR = 2 pc), and modern constants (R₀ = 8.178 kpc, Θ₀ = 235.1 km/s).
They publish the RC as ASCII tables (3 tables in the arXiv source), fit a BH + Plummer
bulge/disc + NFW halo decomposition, extract an E/W asymmetry curve (sinusoid, amplitude
decaying with R — weak-bar interpretation, plus an unexplained "sawtooth" fine structure),
and quote a local dark-matter density **ρ_DM = 0.107 ± 0.003 GeV/cm³** — which the paper
itself flags as ~3× below the ~0.3 GeV/cm³ consensus, explains via its monotonically-declining
outer fit, and frames as a **lower limit** (NFW halo component only; disk DM not separated).

No independent replication or critique exists (checked 2026-08-01; the paper is cited but
only used, not audited). The slice: replicate the HI leg end-to-end from public data, anchor
on their own published tables, and map the **sensitivity of ρ_DM to the decomposition
choices** — the analysis their own lower-limit framing invites. This also retires our own
`hi` slice's documented weakness: the 2 K threshold estimator reads ~7 km/s high
(McClure-Griffiths & Dickey 2016); implementing the proper Gaussian-decomposition TVM and
measuring threshold-vs-decomposition on the same sightlines closes that caveat with a number.
Public hook: this is the professional reference behind the IEEE Spectrum backyard-DM article
and the live HN thread on it.

**GATE 0 (done).** *Data:* HI4PI at CDS `J/A+A/594/A116` (verified reachable; `CUBES/GAL/`
CAR-projection ~20° tiles per `cubes_gal.dat` — the |b|≈0 row across the inner Galaxy is
~9–11 tiles, ~15–20 GB total, fetch-search-delete friendly); CO surveys public (CfA/Dame
2001; FUGIN via JVO; Mopra DOI 10.25919/9z4p-mj92) but **out of minimum scope**; the paper's
RC ASCII tables ship in the arXiv e-print source (anchor data). VERA/VLBA/Gaia comparison
RCs are published tables. *Novelty:* no replication/critique found; ρ_DM meaning pinned from
§5.4 full text (lower limit, halo-only — the slice frames sensitivity, not refutation).

## Deliverables

- `src/jansky_research/innerrc.py`:
  - `fetch_hi4pi_tile` (CDS, `# pragma`, resumable, delete-after option) + `lv_diagram`
    (tile → |b|<3° averaged (ℓ, v) diagram).
  - `gaussian_tvm` — per-spectrum Gaussian decomposition, terminal velocity = highest-velocity
    component above threshold; `dispersion_correction` (±σ_v, iterated to meet Θ₀ at R₀ per
    their §3.3); `threshold_tvm` (our old 2 K estimator, for the head-to-head).
  - `rotation_curve_weighted` — their eq. 13–15 Gaussian-weighted binning with errors.
  - `parse_paper_tables` (arXiv-source ASCII tables → anchor arrays).
  - `decompose_rc` — BH point mass + Plummer bulge/disc + NFW halo least-squares (their
    eq. 19–27); `rho_dm_local` + `sensitivity_scan` (halo profile choice, disc model, outer-RC
    data in/out, fit range → ρ_DM spread).
  - `ew_asymmetry` — E vs W curves + their eq. 29 sinusoid fit + residual (sawtooth) check.
  - `synthetic_lv` — synthetic LV diagram with known injected RC + dispersion for the offline
    round-trip; `run/_figure/_write_macros/_main`.
- **Committed-real-results pattern from day one** (post-2026-07-31 integrity rule): the real
  HI4PI-leg outputs land in force-tracked `results/innerrc_*.json`, and paper macros/figures
  are generated from the committed real results — the offline synthetic leg feeds tests only
  and never writes the paper inputs.
- Tests to the 85% floor; `papers/innerrc/`; `survey/innerrc-findings.md`; wiring
  (Snakefile via the committed-results path, README, CHANGELOG).

## Approach

1. Tooling + synthetic recover-a-known: inject a known RC into `synthetic_lv`; the full
   Gaussian-TVM + dispersion-correction + weighted-binning chain must recover it; the
   threshold estimator must show its known high bias on the same synthetic data.
2. **Anchor leg (mandatory before variants):** HI4PI tiles → our TVM RC vs their Table 2
   (δR = 50 pc). Agreement within their quoted errors licenses everything after; a failure is
   itself the finding and narrows scope to the methods comparison.
3. Estimator head-to-head on the real sightlines: Gaussian-decomposition vs 2 K threshold —
   the number that retires the `hi` slice caveat.
4. E/W asymmetry: reproduce their sinusoid + check the sawtooth residual structure they flag
   as unexplained (fresh eyes on their own residuals; report, no over-interpretation).
5. Decomposition + ρ_DM sensitivity scan: reproduce 0.107 GeV/cm³ from their tables/our RC,
   then vary halo profile (NFW vs Burkert vs isothermal), disc treatment, and outer-RC
   choices; report the ρ_DM range against the ~0.3 consensus — quantifying their lower-limit
   caveat rather than contradicting it.
6. GATE-2 science review: axisymmetry/bar caveats (their §5.6–5.7 and Chemin et al. 2015
   TVM-in-bar critique must be carried honestly), estimator-comparison fairness, ρ_DM framing.
7. Paper (`papers/innerrc/`, full-length; RNAAS-able core if the anchor + sensitivity story
   is compact). CO surveys and the unified-25 kpc fit are stretch legs, not dependencies.

## Verification

Synthetic round-trip recovers the injected RC; anchor reproduces their published table within
errors before any variant runs; the guard (`make guard-real`) plus committed-real-results
wiring keeps synthetic output out of the paper by construction; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Gaussian decomposition is fiddly in crowded spectra** (their Fig. 2 shows ~10-component
  fits) → deterministic initialization from peak-finding, seeded, tested on synthetic crowded
  spectra; report decomposition-failure rate per longitude rather than hand-tuning.
- **Their dispersion correction is calibrated circularly** (iterated to meet Θ₀ at R₀ — their
  own description) → implement exactly as specified for the anchor, then show the RC with and
  without it; the difference is part of the sensitivity story.
- **Bar-region non-circularity** (their §5.6, Chemin+2015): inner ~2 kpc is not a mass tracer
  → same exclusion honesty as the `hi` slice; the sensitivity scan fits R > 2 kpc variants.
- **HI4PI tile volume** (~15–20 GB): fetch-search-delete, one tile at a time; bandwidth-aware
  scheduling per the release-day protocol.
