# 40 — SBI population inference for the RACS Stokes-V coherent-emitter class

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) + `sbi` v0.26.x
ROCm smoke test on the pinned torch wheel

## Context

Neural simulation-based inference has hit pulsars (arXiv:2312.14848, 2412.04070), magnetars
(arXiv:2503.11875), FRB selection functions (arXiv:2606.26334), and per-source QU-fitting
(VROOM-SBI, arXiv:2605.27538) — but no SBI population inference exists for circularly-polarized
stellar/coherent emitters. This repo already owns the two hard inputs: a validated forced-V
measurement (`stokesv.py`, `measure_circular_pol`) and a per-field leakage-floor model — which is
effectively the selection function — plus the real detection/non-detection target list from the
merged `stokesv_discovery` census. Data: RACS-low/mid V via CASDA (verified access pattern),
SRSC (VizieR `J/other/PASA/41.84`) as labelled positives, Gaia CNS5/DR3 M-dwarfs as the parent
sample. The deliverable is the first calibrated posterior on the population parameters, not a
tight number — posterior widths are reported honestly. GPU: NPE training is exactly what the
16 GiB card is for; pure-PyTorch `sbi` only (avoid LtU-ILI's TF backend; no CUDA kernels).

## Deliverables

- `src/jansky_research/svsbi.py`: `draw_population` (luminosity-function slope/break, beaming
  fraction; distances from the Gaia parent sample), `forward_model` (fold through RACS
  sensitivity + the per-beam leakage floor), `fetch_parent_sample` + `fetch_racs_cutouts`
  (# pragma), `train_npe` (`sbi` NPE, ROCm), `posterior_coverage` (simulation-based calibration),
  `inject_recover_v` (synthetic V populations into real cutouts via `measure_circular_pol`),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/svsbi/`; `survey/svsbi-findings.md`; wiring (new dep: `sbi`).

## Approach

0. GATE 0: full-text pass on arXiv:2312.14848, 2412.04070, 2503.11875, 2605.27538, 2606.26334
   to confirm no coherent-emitter population SBI exists; re-verify CASDA V access, the SRSC
   VizieR table, and that `sbi` v0.26.x trains on the pinned ROCm torch wheel (gfx1102).
1. Tooling + synthetic recover-a-known: inject synthetic V populations with known (slope,
   beaming fraction, break) into real RACS cutouts via `measure_circular_pol`; verify NPE
   posterior coverage via simulation-based calibration before any real inference.
2. Real leg (GPU training, days): NPE conditioned on the `stokesv_discovery`
   detections/non-detections; report the posterior with SBC-validated coverage.
3. GATE-2 science review: model misspecification, the beaming/duty-cycle degeneracy at small N,
   leakage-floor fidelity as the selection function, prior sensitivity.
4. Paper: first SBI population inference for coherent circularly-polarized emitters. Do NOT
   also start the two sibling SBI slices (RM-structure-function turbulence; LPT population via
   the arXiv:2509.06315 completeness) — pick one later, per fable-ideas.

## Verification

Simulation-based calibration shows nominal posterior coverage on injected populations; injected
(slope, beaming, break) recovered within stated credible intervals from real-cutout injections;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Model misspecification / beaming–duty-cycle degeneracy at small N** → report posterior
  widths honestly; the deliverable is the first calibrated posterior, not a tight number.
- **ROCm/`sbi` compatibility** → GATE-0 smoke test on the pinned wheel; fall back to CPU
  training (slower but viable) if NPE hits a ROCm wall.
- **Leakage floor mismodelled** → reuse the `stokesv_discovery` per-field vetting; propagate the
  floor uncertainty through the forward model rather than fixing it.
