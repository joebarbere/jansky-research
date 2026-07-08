# Findings — SBI population inference for RACS Stokes-V coherent emitters (plan 40)

`jansky_research.svsbi` builds the first neural simulation-based inference (SBI) of the
population parameters of M-dwarf coherent radio emitters — the radio luminosity function (slope,
break) and the **beaming fraction** — conditioned on the merged `stokesv_discovery` census.

## GATE 0 (2026-07-08)

- **Novelty PASS**: SBI has hit pulsar populations (arXiv:2312.14848, 2412.04070), magnetars
  (arXiv:2503.11875), FRB selection functions (arXiv:2606.26334), and per-source QU-fitting
  (VROOM-SBI, arXiv:2605.27538) — but no SBI *population* inference for coherent stellar
  emitters. All five arXiv IDs verified against their abstracts. Classical prior art is
  qualitative (Callingham+2021, Pritchard+2021, Driessen+2024) — no calibrated beaming-fraction
  posterior exists.
- **`sbi`/ROCm smoke test PASS** (the critical technical gate): `sbi` 0.26.1 is pure PyTorch and
  trains NPE natively on the RX 7600 XT (gfx1102), no code changes — a 2-param toy recovered its
  posterior in 4.7s. Installed into `~/.venvs/rocm-test`.
- **Inputs confirmed**: the `stokesv_discovery` census CSV — 60 targets × 2 RACS-mid epochs, 39
  with usable two-epoch V measurements, 2 confident detections (GJ 65 / CNS5 425 at V~9 mJy),
  median V rms 0.167 mJy, leakage floor ~5.7% (the selection function).

## Design

- **theta = (lf_slope, log10 L_break, f_beam)**; forward model draws each real star a coherent
  luminosity from a power-law-with-exponential-break LF, beams it with prob f_beam, converts to
  V flux at the star's distance (S = L/4πd², 1 erg/s/Hz/cm² = 1e26 mJy), adds per-target V noise
  over 2 epochs (keep the brighter), detects at 5σ AND above the leakage floor. Summary =
  (n_det, log brightest |V|, log summed |V|, detected fraction).
- **Physics is pure NumPy** (draw_population, forward_model, summary_stats) → tested in core CI;
  only train_npe/sbc_ranks touch `sbi` (new `[sbi]` extra, kept out of core CI like `fdmt`).
- Verified the flux scale: GJ 65 (V~10 mJy at ~2.7 pc) ↔ logL ≈ 13.9, so the informative break
  regime is logL* ~ 13.5–15.

## Recover-a-knowns

- **Physics (offline, CI)**: detection count rises monotonically with both f_beam (0.24 → 1.97
  as f_beam 0.05 → 0.50) and log L* (0.69 → 2.36 as 12.5 → 14.5) on a 400-star synthetic parent —
  the pipeline is sensitive to what it infers. Run on a synthetic parent (the real 39-star
  census has too few detections to average the trend out of Poisson noise — a machinery check,
  not a data check).
- **SBC (real leg, GPU)**: 150 injected true θ ranked within their own 150-sample posteriors →
  KS distance from uniform ≤ 0.09 on all three parameters (< 0.15 threshold) → **the posterior's
  credible intervals have nominal coverage**. This is the gate: only a calibrated posterior is
  reported.

## Result (real census, NPE on 12,000 sims, ROCm GPU; post-GATE-2)

- **Beaming fraction f_beam = 0.35, 90% CI [0.11, 0.55]** — the FIRST *calibrated* constraint on
  the fraction of M-dwarf coherent emitters that beam detectably toward us. BUT the interval is
  **0.86 of the prior width** and the median sits near the prior midpoint: one detection barely
  updates the beaming prior. Framed honestly as a calibrated framework, NOT a measurement.
- **LF break is the better-constrained parameter**: log L* = 14.6 (90% CI [13.96, 14.94], ~0.36
  of prior width) — the one bright detection's flux at GJ 65's true 2.7-pc distance pins the
  bright-end scale. Slope 1.29 (CI [1.03, 2.10], ~0.59 of prior). The beaming–luminosity
  degeneracy leaves f_beam (not the luminosity scale) near its prior.

## GATE-2 (PASS with required fixes, all applied)

- **f_beam near-prior-width disclosed (R1)**: was framed as a "measured interval" — corrected to
  a calibrated-but-near-vacuous constraint; the 0.86 posterior/prior width ratio is now a
  pipeline-generated macro in the paper.
- **LF sampler Jacobian fixed (R2)**: the log-grid inverse-CDF was missing the dL∝L measure, so
  the realized slope was slope+1. Added `np.gradient(grid)`; verified input 2.0 → realized 2.03.
- **Unresolved binary deduplicated (R3)**: the "2 detections" were GJ 65's two Gaia components
  (CNS5 424/425, byte-identical forced photometry — the RACS beam can't resolve the ~2″ pair).
  Now n_det = 1 among 38 physical targets.
- **SBC threshold tightened (R4)**: was KS<0.15 (looser than the formal 95% KS critical value
  0.111 at n=150); now uses 1.358/√n. Result still passes: max KS 0.08 < 0.111.
- **Real Gaia distances (S1)**: replaced the fixed 10-pc default with Gaia DR3 parallaxes by
  gaia_id (GJ 65 at 2.7 pc, not 10 — the 10-pc assumption had biased log L* bright by ~1 dex).

## Honest caveats

- **The deliverable is the calibrated posterior + pipeline, not a number.** A one-detection
  census yields a beaming posterior barely narrower than the prior; the value is that its
  intervals are SBC-validated and it tightens directly with a larger census (RACS-low2 V,
  Sydney Radio Star Catalogue).

## Reproduce

Offline (physics + tests): `uv run python -m jansky_research.svsbi --offline --out .`
Real (NPE, ROCm GPU): `PYTHONPATH=src:../jansky/src ~/.venvs/rocm-test/bin/python -m
jansky_research.svsbi --device cuda --n-sims 12000 --out .` (sbi installed in the venv; ~min).
CPU also works (`--device cpu`, slower). New dep: the `[sbi]` extra (`uv sync --extra sbi`).
